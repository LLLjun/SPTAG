// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

#include "inc/Helper/VectorSetReader.h"
#include "inc/Helper/SimpleIniReader.h"
#include "inc/Helper/CommonHelper.h"
#include "inc/Helper/StringConvert.h"
#include "inc/Core/Common/CommonUtils.h"
#include "inc/Core/Common/TruthSet.h"
#include "inc/Core/Common/QueryResultSet.h"
#include "inc/Core/VectorIndex.h"
#include <algorithm>
#include <iomanip>
#include <set>
#include <atomic>
#include <ctime>
#include <thread>
#include <chrono>
#include "inc/Core/SPANN/IExtraSearcher.h"

using namespace SPTAG;

class SearcherOptions : public Helper::ReaderOptions
{
public:
    SearcherOptions() : Helper::ReaderOptions(VectorValueType::Float, 0, VectorFileType::TXT, "|", 32)
    {
        AddRequiredOption(m_queryFile, "-i", "--input", "Input raw data.");
        AddRequiredOption(m_indexFolder, "-x", "--index", "Index folder.");
        AddOptionalOption(m_truthFile, "-r", "--truth", "Truth file.");
        AddOptionalOption(m_resultFile, "-o", "--result", "Output result file.");
        AddOptionalOption(m_maxCheck, "-m", "--maxcheck", "MaxCheck for index.");
        AddOptionalOption(m_withMeta, "-a", "--withmeta", "Output metadata instead of vector id.");
        AddOptionalOption(m_K, "-k", "--KNN", "K nearest neighbors for search.");
        AddOptionalOption(m_truthK, "-tk", "--truthKNN", "truth set number.");
        AddOptionalOption(m_dataFile, "-df", "--data", "original data file.");
        AddOptionalOption(m_dataFileType, "-dft", "--dataFileType", "original data file type. (TXT, or DEFAULT)");
        AddOptionalOption(m_batch, "-b", "--batchsize", "Batch query size.");
        AddOptionalOption(m_genTruth, "-g", "--gentruth", "Generate truth file.");
        AddOptionalOption(m_debugQuery, "-q", "--debugquery", "Debug query number.");
        AddOptionalOption(m_enableADC, "-adc", "--adc", "Enable ADC Distance computation");
    }

    ~SearcherOptions() {}

    std::string m_queryFile;

    std::string m_indexFolder;

    std::string m_dataFile = "";

    std::string m_truthFile = "";

    std::string m_resultFile = "";

    std::string m_maxCheck = "8192";

    VectorFileType m_dataFileType = VectorFileType::DEFAULT;

    int m_withMeta = 0;

    int m_K = 32;

    int m_truthK = -1;

    int m_batch = 10000;

    int m_genTruth = 0;

    int m_debugQuery = -1;

    bool m_enableADC = false;
};

template <typename T>
int Process(std::shared_ptr<SearcherOptions> options, VectorIndex& index)
{
    std::ofstream log("Recall-result.out", std::ios::app);
    if (!log.is_open())
    {
        LOG(Helper::LogLevel::LL_Error, "ERROR: Cannot open logging file!\n");
        exit(-1);
    }

    auto vectorReader = Helper::VectorSetReader::CreateInstance(options);
    if (ErrorCode::Success != vectorReader->LoadFile(options->m_queryFile))
    {
        LOG(Helper::LogLevel::LL_Error, "Failed to read query file.\n");
        exit(1);
    }
    auto queryVectors = vectorReader->GetVectorSet(0, options->m_debugQuery);
    auto queryMetas = vectorReader->GetMetadataSet();

    std::shared_ptr<Helper::ReaderOptions> dataOptions(new Helper::ReaderOptions(queryVectors->GetValueType(), queryVectors->Dimension(), options->m_dataFileType));
    auto dataReader = Helper::VectorSetReader::CreateInstance(dataOptions);
    std::shared_ptr<VectorSet> dataVectors;
    if (options->m_dataFile != "")
    {
        if (ErrorCode::Success != dataReader->LoadFile(options->m_dataFile))
        {
            LOG(Helper::LogLevel::LL_Error, "Failed to read data file.\n");
            exit(1);
        }
        dataVectors = dataReader->GetVectorSet();
    }

    std::shared_ptr<Helper::DiskPriorityIO> ftruth;
    int truthDim = 0;
    if (options->m_truthFile != "")
    {
        if (options->m_genTruth) {
            if (dataVectors == nullptr) {
                LOG(Helper::LogLevel::LL_Error, "Cannot load data vectors to generate groundtruth! Please speicify data vector file by setting -df option.\n");
                exit(1);
            }
            COMMON::TruthSet::GenerateTruth<T>(queryVectors, dataVectors, options->m_truthFile, index.GetDistCalcMethod(), options->m_truthK,
                (options->m_truthFile.find("bin") != std::string::npos) ? TruthFileType::DEFAULT : TruthFileType::TXT);
        }

        ftruth = SPTAG::f_createIO();
        if (ftruth == nullptr || !ftruth->Initialize(options->m_truthFile.c_str(), std::ios::in | std::ios::binary)) {
            LOG(Helper::LogLevel::LL_Error, "ERROR: Cannot open %s for read!\n", options->m_truthFile.c_str());
            exit(1);
        }
        if (options->m_truthFile.find("bin") != std::string::npos) {
            LOG(Helper::LogLevel::LL_Info, "Load binary truth...\n");
        }
        else {
            LOG(Helper::LogLevel::LL_Info, "Load txt truth...\n");
        }
    }

    std::ofstream fp;
    if (options->m_resultFile != "")
    {
        fp.open(options->m_resultFile);
        if (!fp.is_open())
        {
            LOG(Helper::LogLevel::LL_Error, "ERROR: Cannot open %s for write!\n", options->m_resultFile.c_str());
        }
    }

    std::vector<std::string> maxCheck = Helper::StrUtils::SplitString(options->m_maxCheck, "#");
    if (options->m_truthK < 0) options->m_truthK = options->m_K;

    std::vector<std::set<SizeType>> truth(options->m_batch);
    int internalResultNum = options->m_K;
    if (index.GetIndexAlgoType() == IndexAlgoType::SPANN) {
        int SPANNInternalResultNum;
        if (SPTAG::Helper::Convert::ConvertStringTo<int>(index.GetParameter("SearchInternalResultNum", "BuildSSDIndex").c_str(), SPANNInternalResultNum))
            internalResultNum = max(internalResultNum, SPANNInternalResultNum);
    }
    std::vector<QueryResult> results(options->m_batch, QueryResult(NULL, internalResultNum, options->m_withMeta != 0));
    std::vector<float> latencies(options->m_batch, 0);
    int baseSquare = SPTAG::COMMON::Utils::GetBase<T>() * SPTAG::COMMON::Utils::GetBase<T>();

    // For Profiling
    std::vector<SPTAG::SPANN::SearchStats> tatalstats(maxCheck.size());
    std::vector<SPTAG::SPANN::SearchStats> stats(options->m_batch);

    // LOG(Helper::LogLevel::LL_Info, "[query]\t\t[maxcheck]\t[avg] \t[99%%] \t[95%] \t[recall] \t[qps] \t[mem]\n");
    // LOG(Helper::LogLevel::LL_Info, "[query]\t[maxcheck]\t[recall]\t[IO]\t[Access]\t[qps]\t[avg]\t[99%%]\t[extra(us)]\t[search(us)]\n");
    printf("[maxcheck]\t[recall]\t[IO]\t[Access]\t[qps]\t[avg]\t[99%%]\t[extra(us)]\t[search(us)]\n");
    std::vector<float> totalAvg(maxCheck.size(), 0.0), total99(maxCheck.size(), 0.0), total95(maxCheck.size(), 0.0), totalRecall(maxCheck.size(), 0.0), totalLatency(maxCheck.size(), 0.0);
    for (int startQuery = 0; startQuery < queryVectors->Count(); startQuery += options->m_batch)
    {
        int numQuerys = min(options->m_batch, queryVectors->Count() - startQuery);
        for (SizeType i = 0; i < numQuerys; i++) results[i].SetTarget(queryVectors->GetVector(startQuery + i));
        if (ftruth != nullptr) COMMON::TruthSet::LoadTruth(ftruth, truth, numQuerys, truthDim, options->m_truthK, (options->m_truthFile.find("bin") != std::string::npos)? TruthFileType::DEFAULT : TruthFileType::TXT);


        for (int mc = 0; mc < maxCheck.size(); mc++)
        {
            if (index.GetIndexAlgoType() == IndexAlgoType::SPANN) {
                // 其实这里是 索引返回的数量
                index.SetMaxCheck(std::atoi(maxCheck[mc].c_str()));
            } else {
                index.SetParameter("MaxCheck", maxCheck[mc].c_str());
            }

            // for (SizeType i = 0; i < numQuerys; i++) results[i].Reset();
            for (SizeType i = 0; i < numQuerys; i++) results[i].Resize(std::atoi(maxCheck[mc].c_str()));

            for (SizeType i = 0; i < numQuerys; i++) stats[i] = SPTAG::SPANN::SearchStats(std::atoi(maxCheck[mc].c_str()));

            std::atomic_size_t queriesSent(0);
            std::vector<std::thread> threads;
            auto func = [&]()
            {
                size_t qid = 0;
                while (true)
                {
                    qid = queriesSent.fetch_add(1);
                    if (qid < numQuerys)
                    {
                        auto t1 = std::chrono::high_resolution_clock::now();
                        index.SearchIndexStats(results[qid], (void *)(&stats[qid]));
                        auto t2 = std::chrono::high_resolution_clock::now();
                        latencies[qid] = (float)(std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count());
                    }
                    else
                    {
                        return;
                    }
                }
            };

            auto batchstart = std::chrono::high_resolution_clock::now();

            for (std::uint32_t i = 0; i < options->m_threadNum; i++) { threads.emplace_back(func); }
            for (auto& thread : threads) { thread.join(); }

            auto batchend = std::chrono::high_resolution_clock::now();
            float batchLatency = (float)(std::chrono::duration_cast<std::chrono::microseconds>(batchend - batchstart).count());

            float timeMean = 0, timeMin = MaxDist, timeMax = 0, timeStd = 0;
            for (int qid = 0; qid < numQuerys; qid++)
            {
                timeMean += latencies[qid];
                if (latencies[qid] > timeMax) timeMax = latencies[qid];
                if (latencies[qid] < timeMin) timeMin = latencies[qid];
            }
            timeMean /= numQuerys;
            for (int qid = 0; qid < numQuerys; qid++) timeStd += ((float)latencies[qid] - timeMean) * ((float)latencies[qid] - timeMean);
            timeStd = std::sqrt(timeStd / numQuerys);
            log << timeMean << " " << timeStd << " " << timeMin << " " << timeMax << " ";

            std::sort(latencies.begin(), latencies.begin() + numQuerys);
            float l99 = latencies[SizeType(numQuerys * 0.99)];
            float l95 = latencies[SizeType(numQuerys * 0.95)];

            float recall = 0;
            if (ftruth != nullptr)
            {
                recall = COMMON::TruthSet::CalculateRecall<T>(&index, results, truth, options->m_K, options->m_truthK, queryVectors, dataVectors, numQuerys, &log, options->m_debugQuery > 0);
            }

#ifndef _MSC_VER
            struct rusage rusage;
            getrusage(RUSAGE_SELF, &rusage);
            unsigned long long peakWSS = rusage.ru_maxrss * 1024 / 1000000000;
#else
            PROCESS_MEMORY_COUNTERS pmc;
            GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));
            unsigned long long peakWSS = pmc.PeakWorkingSetSize / 1000000000;
#endif
            // LOG(Helper::LogLevel::LL_Info, "%d-%d\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t\t%.4f\t\t%lluGB\n", startQuery, (startQuery + numQuerys), maxCheck[mc].c_str(), timeMean, l99, l95, recall, (numQuerys / batchLatency), peakWSS);
            totalAvg[mc] += timeMean * numQuerys;
            total95[mc] += l95 * numQuerys;
            total99[mc] += l99 * numQuerys;
            totalRecall[mc] += recall * numQuerys;
            totalLatency[mc] += batchLatency;

            // posting list detail
            // for (int qid = 0; qid < numQuerys; qid++){
            //     printf("query id: %d\n", qid);
            //     for (int i = 0; i < stats[qid].dist_num; i++)
            //         printf("%d\t%.3f\t%.3f\n", i, stats[qid].dist_heads[i], stats[qid].dist_post_max[i]);
            //     printf("\n");
            //     if (qid >= 21)
            //         exit(1);
            // }

            // posting list total
            int to = 0;
            for (int qid = 0; qid < numQuerys; qid++){
                float d_bound = stats[qid].dist_heads[9];
                int nu = 0;
                for (int i = 0; i < stats[qid].dist_num; i++){
                    if (stats[qid].dist_heads[i] == 0 && stats[qid].dist_post_max[i] == 0)
                        break;
                    if ((stats[qid].dist_heads[i] - stats[qid].dist_post_max[i]) >= d_bound)
                        nu++;
                }
                to += nu;
                // printf("qid: %d, over: %d\n", qid, nu);
            }
            // printf("total: %d\n", to);
            printf("%d\n", to);

            // exit(0);

            for (int qid = 0; qid < numQuerys; qid++){
                tatalstats[mc].m_exCheck += stats[qid].m_exCheck;
                tatalstats[mc].m_totalListElementsCount += stats[qid].m_totalListElementsCount;
                tatalstats[mc].m_diskIOCount += stats[qid].m_diskIOCount;
                tatalstats[mc].m_diskAccessCount += stats[qid].m_diskAccessCount;

                tatalstats[mc].m_totalSearchLatency += stats[qid].m_totalSearchLatency;
                tatalstats[mc].m_exLatency += stats[qid].m_exLatency;
                // tatalstats[mc].m_totalLatency += stats[qid].m_totalLatency;
            }
        }

        if (fp.is_open())
        {
            fp << std::setprecision(3) << std::fixed;
            for (SizeType i = 0; i < numQuerys; i++)
            {
                if (queryMetas != nullptr) {
                    ByteArray qmeta = queryMetas->GetMetadata(startQuery + i);
                    fp.write((const char*)qmeta.Data(), qmeta.Length());
                }
                else {
                    fp << i;
                }
                fp << ":";
                for (int j = 0; j < options->m_K; j++)
                {
                    if (results[i].GetResult(j)->VID < 0) {
                        fp << results[i].GetResult(j)->Dist << "@NULL" << std::endl;
                        continue;
                    }

                    if (!options->m_withMeta) {
                        fp << (results[i].GetResult(j)->Dist / baseSquare) << "@" << results[i].GetResult(j)->VID << std::endl;
                    }
                    else {
                        ByteArray vm = index.GetMetadata(results[i].GetResult(j)->VID);
                        fp << (results[i].GetResult(j)->Dist / baseSquare) << "@";
                        fp.write((const char*)vm.Data(), vm.Length());
                    }
                    fp << "|";
                }
                fp << std::endl;
            }
        }
    }       // end of batch
    for (int mc = 0; mc < maxCheck.size(); mc++)
        // LOG(Helper::LogLevel::LL_Info, "%d-%d\t%s\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", 0, queryVectors->Count(), 
        printf("%s\t%.4f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n",
            maxCheck[mc].c_str(), (totalRecall[mc] / queryVectors->Count()), 
            (1.0 * tatalstats[mc].m_diskIOCount / queryVectors->Count()),
            (1.0 * tatalstats[mc].m_diskAccessCount / queryVectors->Count()),
            (queryVectors->Count() / totalLatency[mc] * 1000000),
            (totalAvg[mc] / queryVectors->Count()), (total99[mc] / queryVectors->Count()),
            (1.0 * tatalstats[mc].m_exLatency / queryVectors->Count()),
            (1.0 * tatalstats[mc].m_totalSearchLatency / queryVectors->Count())
            );

    LOG(Helper::LogLevel::LL_Info, "Output results finish!\n");

    fp.close();
    log.close();
    return 0;
}

int main(int argc, char** argv)
{
    std::shared_ptr<SearcherOptions> options(new SearcherOptions);
    if (!options->Parse(argc - 1, argv + 1))
    {
        exit(1);
    }

    std::shared_ptr<SPTAG::VectorIndex> vecIndex;
    auto ret = SPTAG::VectorIndex::LoadIndex(options->m_indexFolder, vecIndex);
    if (SPTAG::ErrorCode::Success != ret || nullptr == vecIndex)
    {
        LOG(Helper::LogLevel::LL_Error, "Cannot open index configure file!");
        return -1;
    }
    if (SPTAG::COMMON::DistanceUtils::Quantizer)
    {
        COMMON::DistanceUtils::Quantizer->SetEnableADC(options->m_enableADC);
    }

    Helper::IniReader iniReader;
    for (int i = 1; i < argc; i++)
    {
        std::string param(argv[i]);
        size_t idx = param.find("=");
        if (idx == std::string::npos) continue;

        std::string paramName = param.substr(0, idx);
        std::string paramVal = param.substr(idx + 1);
        std::string sectionName;
        idx = paramName.find(".");
        if (idx != std::string::npos) {
            sectionName = paramName.substr(0, idx);
            paramName = paramName.substr(idx + 1);
        }
        iniReader.SetParameter(sectionName, paramName, paramVal);
        LOG(Helper::LogLevel::LL_Info, "Set [%s]%s = %s\n", sectionName.c_str(), paramName.c_str(), paramVal.c_str());
    }

    std::string sections[] = { "Base", "SelectHead", "BuildHead", "BuildSSDIndex", "Index" };
    for (int i = 0; i < 5; i++) {
        if (!iniReader.DoesParameterExist(sections[i], "NumberOfThreads")) {
            iniReader.SetParameter(sections[i], "NumberOfThreads", std::to_string(options->m_threadNum));
        }
        for (const auto& iter : iniReader.GetParameters(sections[i]))
        {
            vecIndex->SetParameter(iter.first.c_str(), iter.second.c_str(), sections[i]);
        }
    }

    vecIndex->UpdateIndex();

    switch (options->m_inputValueType)
    {
#define DefineVectorValueType(Name, Type) \
    case VectorValueType::Name: \
        Process<Type>(options, *(vecIndex.get())); \
        break; \

#include "inc/Core/DefinitionList.h"
#undef DefineVectorValueType

    default: break;
    }
    return 0;
}
