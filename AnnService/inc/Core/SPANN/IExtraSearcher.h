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
            NmpTrace(int SetSize = 100000) {
                TraceSet.reserve(SetSize);
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
            void outputTrace(std::string file) {
                std::vector<TraceLine> TraceSetRefine;
                TraceSetRefine.reserve(TraceSet.size());
                allocLogicAddress(TraceSetRefine);

                std::ofstream file_writer(file.c_str());
                for (TraceLine& line: TraceSetRefine) {
                    file_writer << line.Request_Arrival_Time << " "
                                << line.Device_Number << " "
                                << line.Starting_Logical_Sector_Address * sector_per_page << " "
                                << line.Request_Size_In_Sectors * sector_per_page << " "
                                << line.Type_of_Requests << "\n";
                }
                file_writer.close();
                printf("[Trace] Output %lu trace lines to %s success\n", TraceSetRefine.size(), file.c_str());
            }
        private:
            std::vector<TraceLine> TraceSet;
            uint64_t max_page = 2120109;
            uint64_t num_channel = 8;
            uint64_t page_per_channel;
            bool isTypeSector = true;
            uint64_t sector_per_page = 1;

            void allocLogicAddress(std::vector<TraceLine>& TraceSetRefine) {
                std::vector<int> allocChannel(TraceSet.size(), 0);
                std::vector<uint64_t> channelSize(num_channel, 0);
                for (size_t i = 0; i < TraceSet.size(); i++) {
                    uint64_t address = TraceSet[i].Starting_Logical_Sector_Address;
                    if (address > max_page) {printf("Error, logical address is out\n"); exit(1);}

                    int i_channel = (int)(address / page_per_channel);
                    channelSize[i_channel]++;
                    allocChannel[i] = i_channel;
                }
                // select max size channel
                int max_size_channel = std::max_element(channelSize.begin(), channelSize.end()) - channelSize.begin();
                for (size_t i = 0; i < allocChannel.size(); i++) {
                    if (allocChannel[i] == max_size_channel) {
                        TraceSetRefine.push_back(TraceSet[i]);
                        TraceSetRefine.back().Request_Arrival_Time = TraceSetRefine.size();
                    }
                }
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