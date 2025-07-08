#include <cassert>
#include <chrono>
#include "stool/include/io/io.hpp"
#include "stool/include/strings/sa_bwt_lcp.hpp"

#include "stool/include/debug/print.hpp"
#include "stool/include/third_party/cmdline.h"
#include "stool/include/debug/debug.hpp"
#include "libdivsufsort/sa.hpp"
// #include "../module/rlbwt_iterator/src/include/rlbwt_iterator.hpp"
// #include "module/rlbwt_iterator/src/include/bwt.hpp"
#include "../include/include.hpp"
#include <sdsl/bit_vectors.hpp>

int main(int argc, char *argv[])
{
#if DEBUG
    std::cout << "\033[31m";
    std::cout << "DEBUG MODE" << std::endl;
    std::cout << "\033[39m" << std::endl;
#endif

    cmdline::parser p;
    p.add<std::string>("input_file", 'i', "input file name", true);
    p.add<std::string>("output_file", 'o', "output file path (default: input_file_path.mus)", false, "");
    p.add<int>("thread_num", 'p', "thread number", false, -1);

    p.parse_check(argc, argv);
    std::string inputFile = p.get<std::string>("input_file");
    // string mode = p.get<string>("mode");

    std::string outputFile = p.get<std::string>("output_file");
    int thread_num = p.get<int>("thread_num");
    if (thread_num < 0)
    {
        thread_num = std::thread::hardware_concurrency();
    }

    std::cout << "Thread number = " << thread_num << std::endl;

    std::ifstream ifs(inputFile);
    bool inputFileExist = ifs.is_open();
    if (!inputFileExist)
    {
        std::cout << inputFile << " cannot open." << std::endl;
        return -1;
    }

    if (outputFile.size() == 0)
    {
        outputFile = inputFile + ".mus";
    }

    stool::rlbwt2::BWTAnalysisResult analysisResult;
    stool::stnode_on_rlbwt::RLE<uint8_t> rlbwt;
    rlbwt.load(inputFile, analysisResult);

    uint64_t ms_count = 0;
    stool::stnode_on_rlbwt::STTreeAnalysisResult st_result;
    uint64_t ds_memory_usage = 0;

    using RDS = stool::stnode_on_rlbwt::RLEWaveletTree<uint64_t>;
    RDS ds = RDS(&rlbwt);
    ds_memory_usage = ds.get_using_memory();

    std::cout << "Enumerate LCP intervals..." << std::endl;
    stool::stnode_on_rlbwt::SuffixTreeNodes<uint64_t, RDS> stnodeTraverser;
    stnodeTraverser.initialize(thread_num, ds, true);

    auto depth_iterator = stnodeTraverser.begin();
    while (depth_iterator != stnodeTraverser.end())
    {
        //st_result.analyze(stnodeTraverser);
        std::vector<stool::LCPInterval<uint64_t>> r2;

        for (auto node_it = depth_iterator.begin(); node_it != depth_iterator.end(); node_it++)
        {
            bool b = stnodeTraverser.check_maximal_repeat(node_it);
            if (b)
            {
                uint64_t left = stnodeTraverser.get_left(node_it);
                uint64_t right = stnodeTraverser.get_right(node_it);
                stool::LCPInterval<uint64_t> intv(left, right, stnodeTraverser.get_current_lcp());
                std::cout << intv.to_string() << std::endl;
                //buffer.push_back(intv);
                //count++;
                /*
                if (buffer.size() >= 8000)
                {
                    out.write((char *)(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                    buffer.clear();
                }
                */
            }
        }
        depth_iterator++;
    }

    // debug(inputFile);
    // computeMUSs(inputFile, outputFile, thread_num);
}