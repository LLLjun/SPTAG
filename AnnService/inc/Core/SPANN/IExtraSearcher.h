// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

#ifndef _SPTAG_SPANN_IEXTRASEARCHER_H_
#define _SPTAG_SPANN_IEXTRASEARCHER_H_

#include "Options.h"

#include "inc/Core/VectorIndex.h"
#include "inc/Helper/AsyncFileReader.h"

#include <memory>
#include <vector>
#include <chrono>
#include <atomic>

namespace SPTAG {
    namespace SPANN {

        struct SearchStats
        {
            SearchStats()
                : m_check(0),
                m_exCheck(0),
                m_totalListElementsCount(0),
                m_diskIOCount(0),
                m_diskAccessCount(0),
                m_totalSearchLatency(0),
                m_totalLatency(0),
                m_exLatency(0),
                m_asyncLatency0(0),
                m_asyncLatency1(0),
                m_asyncLatency2(0),
                m_queueLatency(0),
                m_sleepLatency(0)
            {
            }

            int m_check;

            int m_exCheck;

            int m_totalListElementsCount;

            int m_diskIOCount;

            int m_diskAccessCount;

            double m_totalSearchLatency;

            double m_totalLatency;

            double m_exLatency;

            double m_asyncLatency0;

            double m_asyncLatency1;

            double m_asyncLatency2;

            double m_queueLatency;

            double m_sleepLatency;

            std::chrono::steady_clock::time_point m_searchRequestTime;

            int m_threadID;
        };

        template<typename T>
        class PageBuffer
        {
        public:
            PageBuffer()
                : m_pageBufferSize(0)
            {
            }

            void ReservePageBuffer(std::size_t p_size)
            {
                if (m_pageBufferSize < p_size)
                {
                    m_pageBufferSize = p_size;
                    m_pageBuffer.reset(static_cast<T*>(PAGE_ALLOC(sizeof(T) * m_pageBufferSize)), [=](T* ptr) { PAGE_FREE(ptr); });
                }
            }

            T* GetBuffer()
            {
                return m_pageBuffer.get();
            }

            std::size_t GetPageSize()
            {
                return m_pageBufferSize;
            }

        private:
            std::shared_ptr<T> m_pageBuffer;

            std::size_t m_pageBufferSize;
        };

        struct ExtraWorkSpace
        {
            ExtraWorkSpace() {}

            ~ExtraWorkSpace() {}

            ExtraWorkSpace(ExtraWorkSpace& other) {
                Initialize(other.m_deduper.MaxCheck(), other.m_deduper.HashTableExponent(), (int)other.m_pageBuffers.size(), (int)(other.m_pageBuffers[0].GetPageSize()), other.m_enableDataCompression);
                m_spaceID = g_spaceCount++;
            }

            void Initialize(int p_maxCheck, int p_hashExp, int p_internalResultNum, int p_maxPages, bool enableDataCompression) {
                m_postingIDs.reserve(p_internalResultNum);
                m_deduper.Init(p_maxCheck, p_hashExp);
                m_processIocp.reset(p_internalResultNum);
                m_pageBuffers.resize(p_internalResultNum);
                for (int pi = 0; pi < p_internalResultNum; pi++) {
                    m_pageBuffers[pi].ReservePageBuffer(p_maxPages);
                }
                m_diskRequests.resize(p_internalResultNum);
                m_enableDataCompression = enableDataCompression;
                if (enableDataCompression) {
                    m_decompressBuffer.ReservePageBuffer(p_maxPages);
                }
            }

            void Initialize(va_list& arg) {
                int maxCheck = va_arg(arg, int);
                int hashExp = va_arg(arg, int);
                int internalResultNum = va_arg(arg, int);
                int maxPages = va_arg(arg, int);
                bool enableDataCompression = bool(va_arg(arg, int));
                Initialize(maxCheck, hashExp, internalResultNum, maxPages, enableDataCompression);
            }

            static void Reset() { g_spaceCount = 0; }

            std::vector<int> m_postingIDs;

            COMMON::OptHashPosVector m_deduper;

            Helper::RequestQueue m_processIocp;

            std::vector<PageBuffer<std::uint8_t>> m_pageBuffers;

            bool m_enableDataCompression;
            PageBuffer<std::uint8_t> m_decompressBuffer;

            std::vector<Helper::AsyncReadRequest> m_diskRequests;

            int m_spaceID;

            static std::atomic_int g_spaceCount;
        };

#ifdef NMP_TRACE
        struct TraceLine {
            uint64_t Request_Arrival_Time;
            int Device_Number;
            uint64_t Starting_Logical_Sector_Address;
            uint64_t Request_Size_In_Sectors;
            int Type_of_Requests;

            TraceLine(uint64_t arrival_time, uint64_t logical_address, uint64_t size): 
                        Request_Arrival_Time(arrival_time), Device_Number(0), Starting_Logical_Sector_Address(logical_address),
                        Request_Size_In_Sectors(size), Type_of_Requests(1) {}
        };
        class NmpTrace {
        public:
            NmpTrace(int SetSize = 30000) {
                TracePool.reserve(SetSize);
                TraceSet.reserve(1000);
                page_per_channel = (uint64_t)ceil(1.0 * max_page / num_channel);
                if (isTypeSector)
                    sector_per_page = 8;
            }
            void addTrace(uint64_t arrival_time, uint64_t logical_address, uint64_t size) {
                uint64_t round_time = TraceSet.size();
                TraceSet.push_back(TraceLine(round_time, logical_address, size));
            }
            void clear() {
                std::vector<TraceLine>().swap(TraceSet);
            }
            void outputTrace(std::string dir_trace, bool is_end=false) {
                // add trace set to trace pool
                TracePool.push_back(TraceSet);
                clear();

                if (is_end) {
                    std::vector<std::vector<TraceLine>> TracePoolRefine;
                    std::vector<size_t> TraceSize;
                    for (std::vector<TraceLine>& trace_set: TracePool) {
                        std::vector<TraceLine> trace_set_refine;
                        allocLogicAddress(trace_set, trace_set_refine);
                        TracePoolRefine.push_back(trace_set_refine);
                        TraceSize.push_back(trace_set_refine.size());
                    }

                    // output info
                    std::string file_info = dir_trace + "/0_info.csv";
                    float t_avg = 1.0 * std::accumulate(TraceSize.begin(), TraceSize.end(), 0) / TraceSize.size();
                    int t_avg_i = findNearest(TraceSize, t_avg);
                    int t_max_i = std::max_element(TraceSize.begin(), TraceSize.end()) - TraceSize.begin();
                    int t_min_i = std::min_element(TraceSize.begin(), TraceSize.end()) - TraceSize.begin();
                    std::ofstream info_writer(file_info.c_str());
                    info_writer << "#channel:" << num_channel << "\n";
                    info_writer << "Type,Value,Range" << "\n";
                    info_writer << "Avg," << TraceSize[t_avg_i] << "," << std::to_string(t_avg_i) << "," << t_avg << "\n";
                    info_writer << "Max," << TraceSize[t_max_i] << "," << std::to_string(t_max_i) << "\n";
                    info_writer << "Min," << TraceSize[t_min_i] << "," << std::to_string(t_min_i) << "\n";
                    info_writer.close();
                    // output trace
                    std::vector<std::pair<std::string, int>> type_pair(3);
                    type_pair[0] = std::make_pair("avg", t_avg_i);
                    type_pair[1] = std::make_pair("max", t_max_i);
                    type_pair[2] = std::make_pair("min", t_min_i);
                    for (std::pair<std::string, int>& type: type_pair) {
                        std::string file_trace = dir_trace + "/" + type.first + ".txt";
                        int cur_batch = type.second;
                        std::ofstream file_writer(file_trace.c_str());
                        for (TraceLine& line: TracePoolRefine[cur_batch]) {
                            file_writer << line.Request_Arrival_Time << " "
                                        << line.Device_Number << " "
                                        << line.Starting_Logical_Sector_Address * sector_per_page << " "
                                        << line.Request_Size_In_Sectors * sector_per_page << " "
                                        << line.Type_of_Requests << "\n";
                        }
                        file_writer.close();
                    }
                    printf("[Trace] Output trace lines to %s success\n", dir_trace.c_str());
                }
            }
        private:
            std::vector<std::vector<TraceLine>> TracePool;
            std::vector<TraceLine> TraceSet;
            // sift1b: 213740788, spacev1b: 147442561
            uint64_t max_page = 213740788;
            uint64_t num_channel = 32;
            uint64_t page_per_channel;
            bool isTypeSector = true;
            uint64_t sector_per_page = 1;

            uint64_t total_page = 0;

            void allocLogicAddress(std::vector<TraceLine>& set, std::vector<TraceLine>& set_refine) {
                std::vector<int> allocChannel(set.size(), 0);
                std::vector<uint64_t> channelSize(num_channel, 0);
                for (size_t i = 0; i < set.size(); i++) {
                    uint64_t address = set[i].Starting_Logical_Sector_Address;
                    if (address > max_page) {printf("Error, logical address is out\n"); exit(1);}

                    int i_channel = (int)(address / page_per_channel);
                    channelSize[i_channel]++;
                    allocChannel[i] = i_channel;
                }
                // select max size channel
                int max_size_channel = std::max_element(channelSize.begin(), channelSize.end()) - channelSize.begin();
                for (size_t i = 0; i < allocChannel.size(); i++) {
                    if (allocChannel[i] == max_size_channel) {
                        set_refine.push_back(set[i]);
                        set_refine.back().Request_Arrival_Time = set_refine.size();
                        set_refine.back().Starting_Logical_Sector_Address = 
                                    set_refine.back().Starting_Logical_Sector_Address % page_per_channel;
                    }
                }
                total_page += set_refine.size();
            }
            int findNearest(std::vector<size_t>& list, float value) {
                int fid = 0;
                float diff = std::numeric_limits<float>::max();
                for (size_t i = 0; i < list.size(); i++) {
                    float tmp_diff = std::abs(value - list[i]);
                    if (tmp_diff < diff) {
                        diff = tmp_diff;
                        fid = i;
                    }
                }
                return fid;
            }
        };
#endif

        class IExtraSearcher
        {
        public:
            IExtraSearcher()
            {
            }

            virtual ~IExtraSearcher()
            {
            }

            virtual bool LoadIndex(Options& p_options) = 0;

            virtual void SearchIndex(ExtraWorkSpace* p_exWorkSpace,
                QueryResult& p_queryResults,
                std::shared_ptr<VectorIndex> p_index,
                SearchStats* p_stats,
                std::set<int>* truth = nullptr,
                std::map<int, std::set<int>>* found = nullptr) = 0;

            virtual bool BuildIndex(std::shared_ptr<Helper::VectorSetReader>& p_reader, 
                std::shared_ptr<VectorIndex> p_index, 
                Options& p_opt) = 0;

            virtual bool CheckValidPosting(SizeType postingID) = 0;
#ifdef NMP_TRACE
            std::shared_ptr<NmpTrace> m_trace = std::make_shared<NmpTrace>();
#endif
        };
    } // SPANN
} // SPTAG

#endif // _SPTAG_SPANN_IEXTRASEARCHER_H_