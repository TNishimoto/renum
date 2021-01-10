#include <cassert>
#include <chrono>
#include "../module/stool/src/io.hpp"
#include "../module/stool/src/sa_bwt_lcp.hpp"

#include "../module/stool/src/print.hpp"
#include "../module/stool/src/cmdline.h"
#include "../module/stool/src/debug.hpp"
#include "../module/libdivsufsort/sa.hpp"
//#include "../module/rlbwt_iterator/src/include/rlbwt_iterator.hpp"
//#include "module/rlbwt_iterator/src/include/bwt.hpp"

#include "stool/src/io.hpp"
#include "stool/src/cmdline.h"
#include "stool/src/debug.hpp"

#include <sdsl/bit_vectors.hpp>

#include "../debug/beller_debug.hpp"
#include "../debug/naive_algorithms.hpp"
#include "../debug/stnode_checker.hpp"
#include "../debug/hyper_debug.hpp"

#include "../stnode_enumerator/application.hpp"
#include "../stnode_enumerator/rlcp_interval_enumerator.hpp"
#include "../rlbwt/io.hpp"
#include "../rlbwt/fpos_data_structure.hpp"
#include "../rlbwt/bwt_decompress.hpp"
#include "../main/common.hpp"

using namespace std;
using namespace stool;

using CHAR = char;
using INDEX = uint64_t;
using LCPINTV = stool::LCPInterval<uint64_t>;


void computeMaximalSubstrings(std::string inputFile, std::string outputFile, string mode, int thread_num)
{

    auto start = std::chrono::system_clock::now();
    string bit_size_mode = "UINT64_t";
    std::chrono::system_clock::time_point mid;

    std::ofstream out(outputFile, std::ios::out | std::ios::binary);
    if (!out)
    {
        throw std::runtime_error("Cannot open the output file!");
    }
    stool::rlbwt2::BWTAnalysisResult analysisResult;
    stool::lcp_on_rlbwt::RLE<uint8_t> rlbwt;
    rlbwt.load(inputFile, analysisResult);

    char MODE = mode[0];

    uint64_t ms_count = 0;
    stool::lcp_on_rlbwt::STTreeAnalysisResult st_result;


    if (analysisResult.str_size < UINT32_MAX - 10)
    {
        using RDS = stool::lcp_on_rlbwt::RLEWaveletTree<uint32_t>;
        RDS ds = RDS(&rlbwt, inputFile);
        mid = std::chrono::system_clock::now();

        std::cout << "Enumerate Maximal Substrings..." << std::endl;
        stool::lcp_on_rlbwt::SuffixTreeNodes<uint32_t, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds, false);
        ms_count = stool::lcp_on_rlbwt::Application::outputMaximalSubstrings(out, stnodeTraverser, st_result);
        bit_size_mode = "UINT32_t";
    }
    else
    {
        using RDS = stool::lcp_on_rlbwt::RLEWaveletTree<uint64_t>;
        RDS ds = RDS(&rlbwt, inputFile);
        mid = std::chrono::system_clock::now();

        std::cout << "Enumerate Maximal Substrings..." << std::endl;
        stool::lcp_on_rlbwt::SuffixTreeNodes<uint64_t, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds, false);
        ms_count = stool::lcp_on_rlbwt::Application::outputMaximalSubstrings(out, stnodeTraverser, st_result);
    }
    if (MODE == '0')
    {
        /*
        std::vector<uint64_t> lpos_vec;
        run_bits.to_vector(lpos_vec);

        assert(diff_char_vec.size() + 1 == lpos_vec.size());
        using LPOSDS = stool::lcp_on_rlbwt::RankSupportVectorWrapper<std::vector<uint64_t>>;
        LPOSDS lpos_vec_wrapper(lpos_vec);

        using FPOSDS = std::vector<uint64_t>;
        using RDS = stool::lcp_on_rlbwt::RLEWaveletTree<INDEX, LPOSDS, FPOSDS>;

        FPOSDS fposds = stool::lcp_on_rlbwt::FPosDataStructure::construct(diff_char_vec, lpos_vec_wrapper);

        RDS ds = RDS(diff_char_vec, wt, lpos_vec_wrapper, fposds);

        std::cout << "Enumerate Maximal Substrings..." << std::endl;
        stool::lcp_on_rlbwt::SuffixTreeNodes<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds, false);
        ms_count = stool::lcp_on_rlbwt::Application::outputMaximalSubstrings(out, stnodeTraverser, st_result);
        */
    }
    else
    {
        /*
        stool::EliasFanoVector lpos_vec;
        lpos_vec.build_from_builder(run_bits);
        std::cout << "LPOS Vec using memory = " << lpos_vec.get_using_memory() / 1000 << "[KB]" << std::endl;
        data_structure_bytes += lpos_vec.get_using_memory();

        //lpos_vec.build_from_bit_vector(run_bits);
        using LPOSDS = stool::EliasFanoVector;
        using FPOSDS = stool::lcp_on_rlbwt::LightFPosDataStructure;
        FPOSDS fposds = stool::lcp_on_rlbwt::LightFPosDataStructure(diff_char_vec, lpos_vec, wt);
        std::cout << "FPOS Vec using memory = " << fposds.get_using_memory() / 1000 << "[KB]" << std::endl;
        data_structure_bytes += fposds.get_using_memory();

        mid = std::chrono::system_clock::now();
        construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start).count();
        std::cout << "Construction time: " << construction_time << "[ms]" << std::endl;
        std::cout << "Data structure Size \t\t\t : " << (data_structure_bytes / 1000) << "[KB]" << std::endl;
        */
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double enumeration_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - mid).count();
    double construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start).count();

    double bps = ((double)analysisResult.str_size / ((double)elapsed / 1000)) / 1000;

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "RLBWT File \t\t\t\t : " << inputFile << std::endl;
    std::cout << "Output \t\t\t\t\t : " << outputFile << std::endl;
    std::cout << "LPOS and FPos Vector type \t\t\t\t : " << (MODE == '0' ? "std::vector<uint64_t>" : "EliasFano") << std::endl;
    std::cout << "Peak children count \t\t\t : " << st_result.max_nodes_at_level << std::endl;
    //std::cout << "Data structure Size \t\t\t : " << (data_structure_bytes / 1000) << "[KB]" << std::endl;

    std::cout << "Thread number \t\t\t\t : " << thread_num << std::endl;
    std::cout << "Integer Type \t\t\t\t : " << bit_size_mode << std::endl;

    std::cout << "The length of the input text \t\t : " << analysisResult.str_size << std::endl;
    std::cout << "The number of runs on BWT \t\t : " << analysisResult.run_count << std::endl;
    std::cout << "The number of maximum substrings \t : " << ms_count << std::endl;
    std::cout << "Excecution time \t\t\t : " << elapsed << "[ms]" << std::endl;
    std::cout << "Character per second \t\t\t : " << bps << "[KB/s]" << std::endl;
    std::cout << "\t Preprocessing time \t\t : " << construction_time << "[ms]" << std::endl;
    std::cout << "\t Enumeration time \t\t : " << enumeration_time << "[ms]" << std::endl;

    std::cout << "_______________________________________________________" << std::endl;

    std::cout << "\033[39m" << std::endl;
}

int main(int argc, char *argv[])
{
#if DEBUG
    std::cout << "\033[31m";
    std::cout << "DEBUG MODE" << std::endl;
    std::cout << "\033[39m" << std::endl;
#endif

    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file name", true);
    //p.add<bool>("mode", 'm', "mode", false, false);
    p.add<string>("output_file", 'o', "output file name", false, "");
    p.add<string>("mode", 'm', "mode", false, "1");
    p.add<int>("thread_num", 'p', "thread number", false, -1);

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    string mode = p.get<string>("mode");

    string outputFile = p.get<string>("output_file");
    string format = "binary";
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
        if (format == "csv")
        {
            outputFile = inputFile + ".max.csv";
        }
        else
        {
            outputFile = inputFile + ".max";
        }
    }

    computeMaximalSubstrings(inputFile, outputFile, mode, thread_num);
}
