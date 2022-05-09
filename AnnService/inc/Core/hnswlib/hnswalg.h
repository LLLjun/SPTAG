#pragma once

#include "visited_list_pool.h"
#include "hnswlib.h"
#include <atomic>
#include <random>
#include <stdlib.h>
#include <assert.h>
#include <unordered_set>
#include <list>
#include <map>
#include <stack>
#include "config.h"
#include "omp.h"

#include "inc/Core/VectorIndex.h"

class Timer {
    typedef std::chrono::high_resolution_clock _clock;
    std::chrono::time_point<_clock> check_point;

public:
    Timer() : check_point(_clock::now()) {
    }

    void reset() {
        check_point = _clock::now();
    }

    double getElapsedTimeus() const {
        return (double) std::chrono::duration_cast<std::chrono::microseconds>(
                   _clock::now() - check_point).count();
    }

    double getElapsedTimes() const {
        return getElapsedTimeus() * 1e-6;
    }
};

// load file. store format: (uint32_t)num, (uint32_t)dim, (data_T)num * dim.
template<typename data_T>
void LoadBinToArray(std::string& file_path, data_T *data_m, uint32_t nums, uint32_t dims, bool non_header = false){
    std::ifstream file_reader(file_path.c_str(), std::ios::binary);
    if (!non_header){
        uint32_t nums_r, dims_r;
        file_reader.read((char *) &nums_r, sizeof(uint32_t));
        file_reader.read((char *) &dims_r, sizeof(uint32_t));
        if ((nums != nums_r) || (dims != dims_r)){
            printf("Error, file %s is error, nums_r: %u, dims_r: %u\n", file_path.c_str(), nums_r, dims_r);
            exit(1);
        }
    }

    for (int i = 0; i < nums; i++)
        file_reader.read((char *) (data_m + dims * i), dims * sizeof(data_T));
    file_reader.close();
    printf("Load %u * %u Data from %s done.\n", nums, dims, file_path.c_str());
};

template<typename data_T>
uint32_t compArrayCenter(const data_T *data_m, uint32_t nums, uint32_t dims){
    std::cout << "Comput the center point: \n";
    float *sum_m = new float[dims]();
    float *avg_m = new float[dims]();
    for (size_t i = 0; i < nums; i++){
        for (size_t j = 0; j < dims; j++){
            sum_m[j] += data_m[i * dims + j];
        }
    }
    for (size_t j = 0; j < dims; j++){
        avg_m[j] = sum_m[j] / nums;
    }

    float cur_max = std::numeric_limits<float>::max();
    uint32_t center_pt_id = 0;
// #pragma omp parallel for
    for (size_t i = 0; i < nums; i++){
        float tmp_sum = 0;
        for (size_t j = 0; j < dims; j++){
            tmp_sum += powf((data_m[i * dims + j] - avg_m[j]), 2);
        }
// #pragma omp cratical
        {
            if (tmp_sum < cur_max){
                cur_max = tmp_sum;
                center_pt_id = i;
            }
        }
    }
    std::cout << center_pt_id << "\n";
    return center_pt_id;
};

namespace hnswlib {
    using namespace SPTAG;

    typedef int32_t tableint;
    typedef uint32_t linklistsizeint;

    template<typename T=float>
    class HierarchicalNSW {
    public:
        static const tableint max_update_element_locks = 65536;

        HierarchicalNSW(const std::string &location, bool nmslib = false, size_t max_elements=0) {
            loadIndex(location, max_elements);
        }

        HierarchicalNSW(size_t max_elements, size_t dim, size_t M = 16, size_t ef_construction = 200, size_t random_seed = 100) :
                link_list_locks_(max_elements), link_list_update_locks_(max_update_element_locks), element_levels_(max_elements) {
            max_elements_ = max_elements;

            has_deletions_=false;

            vecdim = dim;
            data_size_ = vecdim * sizeof(T);

            M_ = M;
            maxM_ = M_;
            maxM0_ = M_ * 2;
            ef_construction_ = std::max(ef_construction,M_);
            ef_ = 10;

            level_generator_.seed(random_seed);
            update_probability_generator_.seed(random_seed + 1);

            size_links_level0_ = maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);
            size_data_per_element_ = size_links_level0_ + data_size_ + sizeof(labeltype);
            offsetData_ = size_links_level0_;
            label_offset_ = size_links_level0_ + data_size_;
            offsetLevel0_ = 0;

            data_level0_memory_ = (char *) malloc(max_elements_ * size_data_per_element_);
            if (data_level0_memory_ == nullptr)
                throw std::runtime_error("Not enough memory");

            cur_element_count = 0;

            visited_list_pool_ = new VisitedListPool(1, max_elements);

            //initializations for special treatment of the first node
            enterpoint_node_ = -1;
            maxlevel_ = -1;

            linkLists_ = (char **) malloc(sizeof(void *) * max_elements_);
            if (linkLists_ == nullptr)
                throw std::runtime_error("Not enough memory: HierarchicalNSW failed to allocate linklists");
            size_links_per_element_ = maxM_ * sizeof(tableint) + sizeof(linklistsizeint);
            mult_ = 1 / log(1.0 * M_);
            revSize_ = 1.0 / mult_;

        }

        struct CompareByFirst {
            constexpr bool operator()(std::pair<float, tableint> const &a,
                                      std::pair<float, tableint> const &b) const noexcept {
                return a.first < b.first;
            }
        };

        ~HierarchicalNSW() {

            free(data_level0_memory_);
            for (tableint i = 0; i < cur_element_count; i++) {
                if (element_levels_[i] > 0)
                    free(linkLists_[i]);
            }
            free(linkLists_);
            delete visited_list_pool_;
        }

        size_t max_elements_;
        size_t cur_element_count;
        size_t size_data_per_element_;
        size_t size_links_per_element_;

        size_t M_;
        size_t maxM_;
        size_t maxM0_;
        size_t ef_construction_;

        double mult_, revSize_;
        int maxlevel_;


        VisitedListPool *visited_list_pool_;
        std::mutex cur_element_count_guard_;

        std::vector<std::mutex> link_list_locks_;

        // Locks to prevent race condition during update/insert of an element at same time.
        // Note: Locks for additions can also be used to prevent this race condition if the querying of KNN is not exposed along with update/inserts i.e multithread insert/update/query in parallel.
        std::vector<std::mutex> link_list_update_locks_;
        tableint enterpoint_node_;


        size_t size_links_level0_;
        size_t offsetData_, offsetLevel0_;


        char *data_level0_memory_;
        char **linkLists_;
        std::vector<int> element_levels_;

        size_t data_size_;

        bool has_deletions_;


        size_t label_offset_;
        size_t vecdim;

        std::unordered_map<labeltype, tableint> label_lookup_;

        std::default_random_engine level_generator_;
        std::default_random_engine update_probability_generator_;


        inline labeltype getExternalLabel(tableint internal_id) const {
            labeltype return_label;
            memcpy(&return_label,(data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_), sizeof(labeltype));
            return return_label;
        }

        inline void setExternalLabel(tableint internal_id, labeltype label) const {
            memcpy((data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_), &label, sizeof(labeltype));
        }

        inline labeltype *getExternalLabeLp(tableint internal_id) const {
            return (labeltype *) (data_level0_memory_ + internal_id * size_data_per_element_ + label_offset_);
        }

        inline char *getDataByInternalId(tableint internal_id) const {
            return (data_level0_memory_ + internal_id * size_data_per_element_ + offsetData_);
        }

        int getRandomLevel(double reverse_size) {
            std::uniform_real_distribution<double> distribution(0.0, 1.0);
            double r = -log(distribution(level_generator_)) * reverse_size;
            return (int) r;
        }

        static const unsigned char DELETE_MARK = 0x01;
//        static const unsigned char REUSE_MARK = 0x10;
        /**
         * Marks an element with the given label deleted, does NOT really change the current graph.
         * @param label
         */
        void markDelete(labeltype label)
        {
            has_deletions_=true;
            auto search = label_lookup_.find(label);
            if (search == label_lookup_.end()) {
                throw std::runtime_error("Label not found");
            }
            markDeletedInternal(search->second);
        }

        /**
         * Uses the first 8 bits of the memory for the linked list to store the mark,
         * whereas maxM0_ has to be limited to the lower 24 bits, however, still large enough in almost all cases.
         * @param internalId
         */
        void markDeletedInternal(tableint internalId) {
            unsigned char *ll_cur = ((unsigned char *)get_linklist0(internalId))+2;
            *ll_cur |= DELETE_MARK;
        }

        /**
         * Remove the deleted mark of the node.
         * @param internalId
         */
        void unmarkDeletedInternal(tableint internalId) {
            unsigned char *ll_cur = ((unsigned char *)get_linklist0(internalId))+2;
            *ll_cur &= ~DELETE_MARK;
        }

        /**
         * Checks the first 8 bits of the memory to see if the element is marked deleted.
         * @param internalId
         * @return
         */
        bool isMarkedDeleted(tableint internalId) const {
            unsigned char *ll_cur = ((unsigned char*)get_linklist0(internalId))+2;
            return *ll_cur & DELETE_MARK;
        }

        std::priority_queue<std::pair<float, tableint>, std::vector<std::pair<float, tableint>>, CompareByFirst>
        searchBaseLayer(tableint ep_id, std::shared_ptr<VectorIndex> p_index, const void *data_point, int layer) {
            VisitedList *vl = visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;

            std::priority_queue<std::pair<float, tableint>, std::vector<std::pair<float, tableint>>, CompareByFirst> top_candidates;
            std::priority_queue<std::pair<float, tableint>, std::vector<std::pair<float, tableint>>, CompareByFirst> candidateSet;

            float lowerBound;
            if (!isMarkedDeleted(ep_id)) {
                float dist = p_index->ComputeDistance(data_point, getDataByInternalId(ep_id));
                top_candidates.emplace(dist, ep_id);
                lowerBound = dist;
                candidateSet.emplace(-dist, ep_id);
            } else {
                lowerBound = std::numeric_limits<float>::max();
                candidateSet.emplace(-lowerBound, ep_id);
            }
            visited_array[ep_id] = visited_array_tag;

            while (!candidateSet.empty()) {
                std::pair<float, tableint> curr_el_pair = candidateSet.top();
                if ((-curr_el_pair.first) > lowerBound) {
                    break;
                }
                candidateSet.pop();

                tableint curNodeNum = curr_el_pair.second;

                std::unique_lock <std::mutex> lock(link_list_locks_[curNodeNum]);

                int *data;// = (int *)(linkList0_ + curNodeNum * size_links_per_element0_);
                if (layer == 0) {
                    data = (int*)get_linklist0(curNodeNum);
                } else {
                    data = (int*)get_linklist(curNodeNum, layer);
//                    data = (int *) (linkLists_[curNodeNum] + (layer - 1) * size_links_per_element_);
                }
                size_t size = getListCount((linklistsizeint*)data);
                tableint *datal = (tableint *) (data + 1);
#ifdef USE_SSE
                _mm_prefetch((char *) (visited_array + *(data + 1)), _MM_HINT_T0);
                _mm_prefetch((char *) (visited_array + *(data + 1) + 64), _MM_HINT_T0);
                _mm_prefetch(getDataByInternalId(*datal), _MM_HINT_T0);
                _mm_prefetch(getDataByInternalId(*(datal + 1)), _MM_HINT_T0);
#endif

                for (size_t j = 0; j < size; j++) {
                    tableint candidate_id = *(datal + j);
//                    if (candidate_id == 0) continue;
#ifdef USE_SSE
                    _mm_prefetch((char *) (visited_array + *(datal + j + 1)), _MM_HINT_T0);
                    _mm_prefetch(getDataByInternalId(*(datal + j + 1)), _MM_HINT_T0);
#endif
                    if (visited_array[candidate_id] == visited_array_tag) continue;
                    visited_array[candidate_id] = visited_array_tag;
                    char *currObj1 = (getDataByInternalId(candidate_id));

                    float dist1 = p_index->ComputeDistance(data_point, currObj1);
                    if (top_candidates.size() < ef_construction_ || lowerBound > dist1) {
                        candidateSet.emplace(-dist1, candidate_id);
#ifdef USE_SSE
                        _mm_prefetch(getDataByInternalId(candidateSet.top().second), _MM_HINT_T0);
#endif

                        if (!isMarkedDeleted(candidate_id))
                            top_candidates.emplace(dist1, candidate_id);

                        if (top_candidates.size() > ef_construction_)
                            top_candidates.pop();

                        if (!top_candidates.empty())
                            lowerBound = top_candidates.top().first;
                    }
                }
            }
            visited_list_pool_->releaseVisitedList(vl);

            return top_candidates;
        }

        mutable std::atomic<long> metric_distance_computations;
        mutable std::atomic<long> metric_hops;
        mutable std::atomic<long> metric_hops_L;

        void getNeighborsByHeuristic2(std::shared_ptr<VectorIndex> p_index,
                std::priority_queue<std::pair<float, tableint>, std::vector<std::pair<float, tableint>>, CompareByFirst> &top_candidates,
        const size_t M) {
            if (top_candidates.size() < M) {
                return;
            }

            std::priority_queue<std::pair<float, tableint>> queue_closest;
            std::vector<std::pair<float, tableint>> return_list;
            while (top_candidates.size() > 0) {
                queue_closest.emplace(-top_candidates.top().first, top_candidates.top().second);
                top_candidates.pop();
            }

            while (queue_closest.size()) {
                if (return_list.size() >= M)
                    break;
                std::pair<float, tableint> curent_pair = queue_closest.top();
                float dist_to_query = -curent_pair.first;
                queue_closest.pop();
                bool good = true;

                for (std::pair<float, tableint> second_pair : return_list) {
                    float curdist =
                            p_index->ComputeDistance(getDataByInternalId(second_pair.second),
                                         getDataByInternalId(curent_pair.second));
                    if (curdist < dist_to_query) {
                        good = false;
                        break;
                    }
                }
                if (good) {
                    return_list.push_back(curent_pair);
                }
            }

            for (std::pair<float, tableint> curent_pair : return_list) {
                top_candidates.emplace(-curent_pair.first, curent_pair.second);
            }
        }


        linklistsizeint *get_linklist0(tableint internal_id) const {
            return (linklistsizeint *) (data_level0_memory_ + internal_id * size_data_per_element_ + offsetLevel0_);
        };

        linklistsizeint *get_linklist0(tableint internal_id, char *data_level0_memory_) const {
            return (linklistsizeint *) (data_level0_memory_ + internal_id * size_data_per_element_ + offsetLevel0_);
        };

        linklistsizeint *get_linklist(tableint internal_id, int level) const {
            return (linklistsizeint *) (linkLists_[internal_id] + (level - 1) * size_links_per_element_);
        };

        linklistsizeint *get_linklist_at_level(tableint internal_id, int level) const {
            return level == 0 ? get_linklist0(internal_id) : get_linklist(internal_id, level);
        };

        tableint mutuallyConnectNewElement(const void *data_point, tableint cur_c, std::shared_ptr<VectorIndex> p_index,
                                       std::priority_queue<std::pair<float, tableint>, std::vector<std::pair<float, tableint>>, CompareByFirst> &top_candidates,
        int level, bool isUpdate) {
            size_t Mcurmax = level ? maxM_ : maxM0_;
            getNeighborsByHeuristic2(p_index, top_candidates, M_);
            if (top_candidates.size() > M_)
                throw std::runtime_error("Should be not be more than M_ candidates returned by the heuristic");

            std::vector<tableint> selectedNeighbors;
            selectedNeighbors.reserve(M_);
            while (top_candidates.size() > 0) {
                selectedNeighbors.push_back(top_candidates.top().second);
                top_candidates.pop();
            }

            tableint next_closest_entry_point = selectedNeighbors.back();

            {
                linklistsizeint *ll_cur;
                if (level == 0)
                    ll_cur = get_linklist0(cur_c);
                else
                    ll_cur = get_linklist(cur_c, level);

                if (*ll_cur && !isUpdate) {
                    throw std::runtime_error("The newly inserted element should have blank link list");
                }
                setListCount(ll_cur,selectedNeighbors.size());
                tableint *data = (tableint *) (ll_cur + 1);
                for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {
                    if (data[idx] && !isUpdate)
                        throw std::runtime_error("Possible memory corruption");
                    if (level > element_levels_[selectedNeighbors[idx]])
                        throw std::runtime_error("Trying to make a link on a non-existent level");

                    data[idx] = selectedNeighbors[idx];

                }
            }

            for (size_t idx = 0; idx < selectedNeighbors.size(); idx++) {

                std::unique_lock <std::mutex> lock(link_list_locks_[selectedNeighbors[idx]]);

                linklistsizeint *ll_other;
                if (level == 0)
                    ll_other = get_linklist0(selectedNeighbors[idx]);
                else
                    ll_other = get_linklist(selectedNeighbors[idx], level);

                size_t sz_link_list_other = getListCount(ll_other);

                if (sz_link_list_other > Mcurmax)
                    throw std::runtime_error("Bad value of sz_link_list_other");
                if (selectedNeighbors[idx] == cur_c)
                    throw std::runtime_error("Trying to connect an element to itself");
                if (level > element_levels_[selectedNeighbors[idx]])
                    throw std::runtime_error("Trying to make a link on a non-existent level");

                tableint *data = (tableint *) (ll_other + 1);

                bool is_cur_c_present = false;
                if (isUpdate) {
                    for (size_t j = 0; j < sz_link_list_other; j++) {
                        if (data[j] == cur_c) {
                            is_cur_c_present = true;
                            break;
                        }
                    }
                }

                // If cur_c is already present in the neighboring connections of `selectedNeighbors[idx]` then no need to modify any connections or run the heuristics.
                if (!is_cur_c_present) {
                    if (sz_link_list_other < Mcurmax) {
                        data[sz_link_list_other] = cur_c;
                        setListCount(ll_other, sz_link_list_other + 1);
                    } else {
                        // finding the "weakest" element to replace it with the new one
                        float d_max = p_index->ComputeDistance(getDataByInternalId(cur_c), getDataByInternalId(selectedNeighbors[idx]));
                        // Heuristic:
                        std::priority_queue<std::pair<float, tableint>, std::vector<std::pair<float, tableint>>, CompareByFirst> candidates;
                        candidates.emplace(d_max, cur_c);

                        for (size_t j = 0; j < sz_link_list_other; j++) {
                            candidates.emplace(
                                    p_index->ComputeDistance(getDataByInternalId(data[j]), getDataByInternalId(selectedNeighbors[idx])), data[j]);
                        }

                        getNeighborsByHeuristic2(p_index, candidates, Mcurmax);

                        int indx = 0;
                        while (candidates.size() > 0) {
                            data[indx] = candidates.top().second;
                            candidates.pop();
                            indx++;
                        }

                        setListCount(ll_other, indx);
                        // Nearest K:
                        /*int indx = -1;
                        for (int j = 0; j < sz_link_list_other; j++) {
                            float d = p_index->ComputeDistance(getDataByInternalId(data[j]), getDataByInternalId(rez[idx]));
                            if (d > d_max) {
                                indx = j;
                                d_max = d;
                            }
                        }
                        if (indx >= 0) {
                            data[indx] = cur_c;
                        } */
                    }
                }
            }

            return next_closest_entry_point;
        }

        std::mutex global;
        size_t ef_ = 1024;

        void setEf(size_t ef) {
            ef_ = ef;
        }

        void resizeIndex(size_t new_max_elements){
            if (new_max_elements<cur_element_count)
                throw std::runtime_error("Cannot resize, max element is less than the current number of elements");


            delete visited_list_pool_;
            visited_list_pool_ = new VisitedListPool(1, new_max_elements);


            element_levels_.resize(new_max_elements);

            std::vector<std::mutex>(new_max_elements).swap(link_list_locks_);

            // Reallocate base layer
            char * data_level0_memory_new = (char *) realloc(data_level0_memory_, new_max_elements * size_data_per_element_);
            if (data_level0_memory_new == nullptr)
                throw std::runtime_error("Not enough memory: resizeIndex failed to allocate base layer");
            data_level0_memory_ = data_level0_memory_new;

            // Reallocate all other layers
            char ** linkLists_new = (char **) realloc(linkLists_, sizeof(void *) * new_max_elements);
            if (linkLists_new == nullptr)
                throw std::runtime_error("Not enough memory: resizeIndex failed to allocate other layers");
            linkLists_ = linkLists_new;

            max_elements_ = new_max_elements;
        }

        void saveIndex(const std::string &location) {
            std::ofstream output(location, std::ios::binary);
            std::streampos position;

            writeBinaryPOD(output, offsetLevel0_);
            writeBinaryPOD(output, max_elements_);
            writeBinaryPOD(output, cur_element_count);
            writeBinaryPOD(output, size_data_per_element_);
            writeBinaryPOD(output, label_offset_);
            writeBinaryPOD(output, offsetData_);
            writeBinaryPOD(output, maxlevel_);
            writeBinaryPOD(output, enterpoint_node_);
            writeBinaryPOD(output, maxM_);

            writeBinaryPOD(output, maxM0_);
            writeBinaryPOD(output, M_);
            writeBinaryPOD(output, mult_);
            writeBinaryPOD(output, ef_construction_);
            writeBinaryPOD(output, vecdim);

            output.write(data_level0_memory_, cur_element_count * size_data_per_element_);

            for (size_t i = 0; i < cur_element_count; i++) {
                unsigned int linkListSize = element_levels_[i] > 0 ? size_links_per_element_ * element_levels_[i] : 0;
                writeBinaryPOD(output, linkListSize);
                if (linkListSize)
                    output.write(linkLists_[i], linkListSize);
            }
            output.close();
        }

        void loadIndex(const std::string &location, size_t max_elements_i=0) {


            std::ifstream input(location, std::ios::binary);

            if (!input.is_open())
                throw std::runtime_error("Cannot open file");

            // get file size:
            input.seekg(0,input.end);
            std::streampos total_filesize=input.tellg();
            input.seekg(0,input.beg);

            readBinaryPOD(input, offsetLevel0_);
            readBinaryPOD(input, max_elements_);
            readBinaryPOD(input, cur_element_count);

            size_t max_elements=max_elements_i;
            if(max_elements < cur_element_count)
                max_elements = max_elements_;
            max_elements_ = max_elements;
            readBinaryPOD(input, size_data_per_element_);
            readBinaryPOD(input, label_offset_);
            readBinaryPOD(input, offsetData_);
            readBinaryPOD(input, maxlevel_);
            readBinaryPOD(input, enterpoint_node_);

            readBinaryPOD(input, maxM_);
            readBinaryPOD(input, maxM0_);
            readBinaryPOD(input, M_);
            readBinaryPOD(input, mult_);
            readBinaryPOD(input, ef_construction_);
            readBinaryPOD(input, vecdim);

            auto pos=input.tellg();

            /// Optional - check if index is ok:

            input.seekg(cur_element_count * size_data_per_element_,input.cur);
            for (size_t i = 0; i < cur_element_count; i++) {
                if(input.tellg() < 0 || input.tellg()>=total_filesize){
                    throw std::runtime_error("Index seems to be corrupted or unsupported");
                }

                unsigned int linkListSize;
                readBinaryPOD(input, linkListSize);
                if (linkListSize != 0) {
                    input.seekg(linkListSize,input.cur);
                }
            }

            // throw exception if it either corrupted or old index
            if(input.tellg()!=total_filesize)
                throw std::runtime_error("Index seems to be corrupted or unsupported");

            input.clear();

            /// Optional check end

            input.seekg(pos,input.beg);

            data_level0_memory_ = (char *) malloc(max_elements * size_data_per_element_);
            if (data_level0_memory_ == nullptr)
                throw std::runtime_error("Not enough memory: loadIndex failed to allocate level0");
            input.read(data_level0_memory_, cur_element_count * size_data_per_element_);

            size_links_per_element_ = maxM_ * sizeof(tableint) + sizeof(linklistsizeint);

            size_links_level0_ = maxM0_ * sizeof(tableint) + sizeof(linklistsizeint);
            std::vector<std::mutex>(max_elements).swap(link_list_locks_);
            std::vector<std::mutex>(max_update_element_locks).swap(link_list_update_locks_);

            visited_list_pool_ = new VisitedListPool(1, max_elements);

            linkLists_ = (char **) malloc(sizeof(void *) * max_elements);
            if (linkLists_ == nullptr)
                throw std::runtime_error("Not enough memory: loadIndex failed to allocate linklists");
            element_levels_ = std::vector<int>(max_elements);
            revSize_ = 1.0 / mult_;
            ef_ = 10;
            for (size_t i = 0; i < cur_element_count; i++) {
                label_lookup_[getExternalLabel(i)]=i;
                unsigned int linkListSize;
                readBinaryPOD(input, linkListSize);
                if (linkListSize == 0) {
                    element_levels_[i] = 0;

                    linkLists_[i] = nullptr;
                } else {
                    element_levels_[i] = linkListSize / size_links_per_element_;
                    linkLists_[i] = (char *) malloc(linkListSize);
                    if (linkLists_[i] == nullptr)
                        throw std::runtime_error("Not enough memory: loadIndex failed to allocate linklist");
                    input.read(linkLists_[i], linkListSize);
                }
            }

            has_deletions_=false;

            for (size_t i = 0; i < cur_element_count; i++) {
                if(isMarkedDeleted(i))
                    has_deletions_=true;
            }

            input.close();

            return;
        }

        unsigned short int getListCount(linklistsizeint * ptr) const {
            return *((unsigned short int *)ptr);
        }

        void setListCount(linklistsizeint * ptr, unsigned short int size) const {
            *((unsigned short int*)(ptr))=*((unsigned short int *)&size);
        }

        void addPoint(const void *data_point, std::shared_ptr<VectorIndex> p_index, labeltype label) {
            addPoint(data_point, p_index, label,-1);
        }

        tableint addPoint(const void *data_point, std::shared_ptr<VectorIndex> p_index, labeltype label, int level) {

            tableint cur_c = 0;
            {
                // Checking if the element with the same label already exists
                // if so, updating it *instead* of creating a new element.
                std::unique_lock <std::mutex> templock_curr(cur_element_count_guard_);
                auto search = label_lookup_.find(label);
                if (search != label_lookup_.end()) {
                    tableint existingInternalId = search->second;
                    templock_curr.unlock();

                    std::unique_lock <std::mutex> lock_el_update(link_list_update_locks_[(existingInternalId & (max_update_element_locks - 1))]);

                    if (isMarkedDeleted(existingInternalId)) {
                        unmarkDeletedInternal(existingInternalId);
                    }
                    // updatePoint(data_point, existingInternalId, 1.0);
                    printf("Unsupport updatePoint\n");
                    exit(1);

                    return existingInternalId;
                }

                if (cur_element_count >= max_elements_) {
                    throw std::runtime_error("The number of elements exceeds the specified limit");
                };

                cur_c = cur_element_count;
                cur_element_count++;
                label_lookup_[label] = cur_c;
            }

            // Take update lock to prevent race conditions on an element with insertion/update at the same time.
            std::unique_lock <std::mutex> lock_el_update(link_list_update_locks_[(cur_c & (max_update_element_locks - 1))]);
            std::unique_lock <std::mutex> lock_el(link_list_locks_[cur_c]);
            int curlevel;
#if PLATG
            curlevel = 0;
#else
            curlevel = getRandomLevel(mult_);
            if (level > 0)
                curlevel = level;
#endif
            element_levels_[cur_c] = curlevel;


            std::unique_lock <std::mutex> templock(global);
            int maxlevelcopy = maxlevel_;
            if (curlevel <= maxlevelcopy)
                templock.unlock();
            tableint currObj = enterpoint_node_;
            tableint enterpoint_copy = enterpoint_node_;


            memset(data_level0_memory_ + cur_c * size_data_per_element_ + offsetLevel0_, 0, size_data_per_element_);

            // Initialisation of the data and label
            memcpy(getExternalLabeLp(cur_c), &label, sizeof(labeltype));
            memcpy(getDataByInternalId(cur_c), data_point, data_size_);


            if (curlevel) {
                linkLists_[cur_c] = (char *) malloc(size_links_per_element_ * curlevel + 1);
                if (linkLists_[cur_c] == nullptr)
                    throw std::runtime_error("Not enough memory: addPoint failed to allocate linklist");
                memset(linkLists_[cur_c], 0, size_links_per_element_ * curlevel + 1);
            }

            if ((signed)currObj != -1) {

                if (curlevel < maxlevelcopy) {

                    float curdist = p_index->ComputeDistance(data_point, getDataByInternalId(currObj));
                    for (int level = maxlevelcopy; level > curlevel; level--) {


                        bool changed = true;
                        while (changed) {
                            changed = false;
                            unsigned int *data;
                            std::unique_lock <std::mutex> lock(link_list_locks_[currObj]);
                            data = get_linklist(currObj,level);
                            int size = getListCount(data);

                            tableint *datal = (tableint *) (data + 1);
                            for (int i = 0; i < size; i++) {
                                tableint cand = datal[i];
                                if (cand < 0 || cand > max_elements_)
                                    throw std::runtime_error("cand error");
                                float d = p_index->ComputeDistance(data_point, getDataByInternalId(cand));
                                if (d < curdist) {
                                    curdist = d;
                                    currObj = cand;
                                    changed = true;
                                }
                            }
                        }
                    }
                }

                bool epDeleted = isMarkedDeleted(enterpoint_copy);
                for (int level = std::min(curlevel, maxlevelcopy); level >= 0; level--) {
                    if (level > maxlevelcopy || level < 0)  // possible?
                        throw std::runtime_error("Level error");

                    std::priority_queue<std::pair<float, tableint>, std::vector<std::pair<float, tableint>>, CompareByFirst> top_candidates = searchBaseLayer(
                            currObj, p_index, data_point, level);
                    if (epDeleted) {
                        top_candidates.emplace(p_index->ComputeDistance(data_point, getDataByInternalId(enterpoint_copy)), enterpoint_copy);
                        if (top_candidates.size() > ef_construction_)
                            top_candidates.pop();
                    }
                    currObj = mutuallyConnectNewElement(data_point, cur_c, p_index, top_candidates, level, false);
                }


            } else {
                // Do nothing for the first element
                enterpoint_node_ = 0;
                maxlevel_ = curlevel;

            }

            //Releasing lock for the maximum level
            if (curlevel > maxlevelcopy) {
                enterpoint_node_ = cur_c;
                maxlevel_ = curlevel;
            }
            return cur_c;
        };


        /*
            using one queue to search
        */
        struct Neighbor {
            tableint id;
            float distance;
            bool flag;

            Neighbor() = default;
            Neighbor(tableint id, float distance, bool f) : id{id}, distance{distance}, flag(f) {}

            inline bool operator<(const Neighbor &other) const {
                return distance < other.distance;
            }
        };

        static inline int InsertIntoPool(Neighbor *addr, int L, Neighbor nn) {
            // find the location to insert
            int left = 0, right = L - 1;
            if (addr[left].distance > nn.distance){
                memmove((char *)&addr[left + 1], &addr[left], L * sizeof(Neighbor));
                addr[left] = nn;
                return left;
            }
            if (addr[right].distance < nn.distance){
                addr[L] = nn;
                return L;
            }
            while (left < right - 1){
                int mid = (left + right) / 2;
                if (addr[mid].distance > nn.distance)
                    right = mid;
                else
                    left = mid;
            }
            // check equal ID

            while (left > 0){
                if (addr[left].distance < nn.distance)
                    break;
                if (addr[left].id == nn.id)
                    return L + 1;
                left--;
            }
            if (addr[left].id == nn.id || addr[right].id == nn.id)
                return L + 1;
            memmove((char *)&addr[right + 1], &addr[right], (L - right) * sizeof(Neighbor));
            addr[right] = nn;
            return right;
        }


        /*
            支持rank-level的mapping
        */
#if RANKMAP
        // 不同rank内的起始搜索点
        std::vector<tableint> ept_rank;
        std::vector<int> interId_to_rankLabel;
        // 仅仅只是为了转换中心点
        std::vector<std::vector<tableint>> rankId_to_interId;

        void initRankMap(){
            int num_ranks = NUM_RANKS;
            ept_rank.resize(num_ranks);
            interId_to_rankLabel.resize(cur_element_count);
            rankId_to_interId.resize(num_ranks);

            // 图上的点到rank，暂时采用简单的mapping方式
            size_t num_max_rank = ceil(1.0 * cur_element_count / num_ranks);
            size_t num_pad_rank = num_max_rank * num_ranks - cur_element_count;
            std::vector<size_t> offest_rank_start(num_ranks);
            offest_rank_start[0] = 0;
            for (size_t i = 1; i < num_ranks; i++){
                if (i < (num_ranks - num_pad_rank + 1))
                    offest_rank_start[i] = offest_rank_start[i-1] + num_max_rank;
                else
                    offest_rank_start[i] = offest_rank_start[i-1] + (num_max_rank - 1);
            }

            for (tableint i = 0; i < cur_element_count; i++){
                for (int j = (num_ranks - 1); j >= 0; j--){
                    if (i >= offest_rank_start[j]){
                        interId_to_rankLabel[i] = j;
                        rankId_to_interId[j].push_back(i);
                        break;
                    }
                }
            }

            // 每个rank的起始点，暂时采用中心点
            // size_t vecdim = *(size_t *)(dist_func_param_);
            float* mass_comput = new float[num_max_rank * vecdim]();
            for (int i = 0; i < num_ranks; i++){
                unsigned rank_size = rankId_to_interId[i].size();
                for (int j = 0; j < rank_size; j++){
                    tableint cur_inter = rankId_to_interId[i][j];
                    memcpy(mass_comput + j * vecdim, getDataByInternalId(cur_inter), vecdim * sizeof(float));
                }
                unsigned center = compArrayCenter<float>(mass_comput, rank_size, vecdim);
                ept_rank[i] = rankId_to_interId[i][center];
            }
            delete[] mass_comput;
            std::vector<std::vector<tableint>>().swap(rankId_to_interId);

            stats = new QueryStats;
        }

        /*
            detail infomation
        */
        struct QueryStats {
            double hlc_us = 0;
            double rank_us = 0;
            double n_max_NDC = 0;
            double n_hops = 0;
        };

        QueryStats* stats = nullptr;


        /*
            input: query, k
            output: result
        */
        void searchParaRank(QueryResult& p_results, std::shared_ptr<VectorIndex> p_index, int K) const {
            // 相关参数
            // size_t l_search = std::max(ef_, (size_t)K);
            size_t l_search = (size_t)K * 2;
            size_t num_ranks = NUM_RANKS;

            // HLC level
            // visited list
            VisitedList *vl = visited_list_pool_->getFreeVisitedList();
            vl_type *visited_array = vl->mass;
            vl_type visited_array_tag = vl->curV;

            // priority queue
            std::vector<Neighbor> retset(l_search + 1);

            // rank buffer
            tableint* mem_rank_alloc = new tableint[num_ranks * maxM0_]();
            std::vector<std::pair<int, tableint*>> buffer_rank_alloc(num_ranks);
            for (int i = 0; i < num_ranks; i++){
                buffer_rank_alloc[i].first = 0;
                buffer_rank_alloc[i].second = mem_rank_alloc + i * maxM0_;
            }
            std::vector<std::stack<std::pair<float, tableint>>> buffer_rank_gather(num_ranks);

            // 本轮搜索的邻居list信息
            std::pair<size_t, int*> candidate_neighbors;

            // launch stage
            Timer clk_query = Timer();
            for (int i = 0; i < num_ranks; i++){
                tableint currObj = ept_rank[i];
                float curdist = p_index->ComputeDistance(p_results.GetTarget(), getDataByInternalId(currObj));

                buffer_rank_gather[i].push(std::make_pair(curdist, currObj));
            }
            if (stats != nullptr){
                stats->rank_us += clk_query.getElapsedTimeus();
                stats->n_max_NDC++;
                metric_distance_computations += num_ranks;
                clk_query.reset();
            }

            for (int i = 0; i < num_ranks; i++){
                tableint currObj = buffer_rank_gather[i].top().second;
                float curdist = buffer_rank_gather[i].top().first;
                buffer_rank_gather[i].pop();
                visited_array[currObj] = visited_array_tag;

                if (i < l_search){
                    retset[i].id = currObj;
                    retset[i].distance = curdist;
                    retset[i].flag = true;
                } else {
                    sort(retset.begin(), retset.begin() + l_search);
                    if (curdist >= retset[l_search - 1].distance)
                        continue;
                    retset[l_search - 1].id = currObj;
                    retset[l_search - 1].distance = curdist;
                    retset[l_search - 1].flag = true;
                }
            }
            sort(retset.begin(), retset.begin() + std::min(num_ranks, l_search));

            // running stage
            int k = 0;
            int cur_list_size = num_ranks;
            while (k < cur_list_size){
                int nk = cur_list_size;

                if (retset[k].flag){
                    retset[k].flag = false;

                    // 读取候选者的邻居信息
                    tableint current_node_id = retset[k].id;
                    candidate_neighbors.second = (int *) get_linklist0(current_node_id);
                    candidate_neighbors.first = getListCount((linklistsizeint*)candidate_neighbors.second);

                    // 分配邻居to rank
                    for (size_t j = 1; j <= candidate_neighbors.first; j++) {
                        int candidate_id = *(candidate_neighbors.second + j);

                        if (!(visited_array[candidate_id] == visited_array_tag)) {
                            visited_array[candidate_id] = visited_array_tag;
                            int rank_label = interId_to_rankLabel[candidate_id];
                            int len_buffer = buffer_rank_alloc[rank_label].first;
                            buffer_rank_alloc[rank_label].second[len_buffer] = candidate_id;
                            buffer_rank_alloc[rank_label].first++;
                        }
                    }

                    if (stats != nullptr){
                        stats->hlc_us += clk_query.getElapsedTimeus();
                        int n_max = 0;
                        for (std::pair<int, tableint*>& bra: buffer_rank_alloc){
                            n_max = std::max(n_max, bra.first);
                            metric_distance_computations += bra.first;
                        }
                        stats->n_max_NDC += n_max;
                        stats->n_hops++;

                        clk_query.reset();
                    }

                    for (int i = 0; i < num_ranks; i++){
                        for (int j = 0; j < buffer_rank_alloc[i].first; j++) {
                            tableint currObj = buffer_rank_alloc[i].second[j];

                            float curdist = p_index->ComputeDistance(p_results.GetTarget(), getDataByInternalId(currObj));
                            buffer_rank_gather[i].push(std::make_pair(curdist, currObj));
                        }
                        buffer_rank_alloc[i].first = 0;
                    }

                    if (stats != nullptr){
                        stats->rank_us += clk_query.getElapsedTimeus();
                        clk_query.reset();
                    }

                    for (int i = 0; i < num_ranks; i++){
                        while (!buffer_rank_gather[i].empty()){
                            float dist = buffer_rank_gather[i].top().first;
                            tableint candidate_id = buffer_rank_gather[i].top().second;
                            buffer_rank_gather[i].pop();
                            if (dist >= retset[cur_list_size - 1].distance && (cur_list_size == l_search))
                                continue;

                            Neighbor nn(candidate_id, dist, true);
                            int r = InsertIntoPool(retset.data(), cur_list_size, nn);
                            if (cur_list_size < l_search)
                                ++cur_list_size;
                            if (r < nk)
                                nk = r;
                        }
                    }
                }
                if (nk <= k)
                    k = nk;
                else
                    ++k;
            }

            delete[] mem_rank_alloc;
            visited_list_pool_->releaseVisitedList(vl);
            for (int i = 0; i < K; i++) {
                p_results.SetResult(i, (int32_t)getExternalLabel(retset[i].id), retset[i].distance);
            }
        };

#endif

        void checkIntegrity(){
            int connections_checked=0;
            std::vector <int > inbound_connections_num(cur_element_count,0);
            for(int i = 0;i < cur_element_count; i++){
                for(int l = 0;l <= element_levels_[i]; l++){
                    linklistsizeint *ll_cur = get_linklist_at_level(i,l);
                    int size = getListCount(ll_cur);
                    tableint *data = (tableint *) (ll_cur + 1);
                    std::unordered_set<tableint> s;
                    for (int j=0; j<size; j++){
                        assert(data[j] > 0);
                        assert(data[j] < cur_element_count);
                        assert (data[j] != i);
                        inbound_connections_num[data[j]]++;
                        s.insert(data[j]);
                        connections_checked++;

                    }
                    assert(s.size() == size);
                }
            }
            if(cur_element_count > 1){
                int min1=inbound_connections_num[0], max1=inbound_connections_num[0];
                for(int i=0; i < cur_element_count; i++){
                    assert(inbound_connections_num[i] > 0);
                    min1=std::min(inbound_connections_num[i],min1);
                    max1=std::max(inbound_connections_num[i],max1);
                }
                std::cout << "Min inbound: " << min1 << ", Max inbound:" << max1 << "\n";
            }
            std::cout << "integrity ok, checked " << connections_checked << " connections\n";

        }

    };

}
