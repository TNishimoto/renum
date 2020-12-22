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

//#include "hpp/bwt.hpp"
//#include "hpp/fmindex.hpp"
#include "../beller/beller_debug.hpp"
#include "../hyper/application.hpp"
#include "../main/common.hpp"
#include "../test/naive_algorithms.hpp"
#include "../hyper/rlcp_interval_enumerator.hpp"
#include "../rlbwt/io.hpp"
#include <sdsl/bit_vectors.hpp>
#include "../rlbwt/fpos_data_structure.hpp"

#include "../rlbwt/bwt_decompress.hpp"
#include "../rlbwt/light_rlbwt.hpp"
#include "../test/stnode_checker.hpp"
#include "../debug/hyper_debug.hpp"

using namespace std;
using namespace stool;

using CHAR = char;
using INDEX = uint64_t;
using LCPINTV = stool::LCPInterval<uint64_t>;

/*
stool::rlbwt2::BWTAnalysisResult LoadBWT(std::string inputFile, sdsl::int_vector<> &diff_char_vec, stool::EliasFanoVectorBuilder &run_bits)
{
    auto bwtAnalysis = stool::rlbwt2::load_RLBWT_from_file(inputFile, diff_char_vec, run_bits);
    std::cout << "BWT using memory = " << sdsl::size_in_bytes(diff_char_vec) / 1000 << "[KB]" << std::endl;
    std::cout << "Run bits using memory = " << run_bits.get_using_memory() / 1000 << "[KB]" << std::endl;
    return bwtAnalysis;
}
template<typename T>
void constructDataStructure(sdsl::int_vector<> &diff_char_vec, stool::EliasFanoVectorBuilder &run_bits){

}
*/

void computeMaximalSubstrings(std::string inputFile, std::string outputFile, string mode, int thread_num)
{
    auto start = std::chrono::system_clock::now();

    std::ofstream out(outputFile, std::ios::out | std::ios::binary);
    if (!out)
    {
        throw std::runtime_error("Cannot open the output file!");
    }

    sdsl::int_vector<> diff_char_vec;
    stool::EliasFanoVectorBuilder run_bits;
    //std::vector<bool> run_bits;
    auto bwtAnalysis = stool::rlbwt2::load_RLBWT_from_file(inputFile, diff_char_vec, run_bits);
    std::cout << "BWT using memory = " << sdsl::size_in_bytes(diff_char_vec) / 1000 << "[KB]" << std::endl;
    std::cout << "Run bits using memory = " << run_bits.get_using_memory() / 1000 << "[KB]" << std::endl;

    sdsl::store_to_file(diff_char_vec, inputFile + ".tmp");

    stool::WT wt;
    construct(wt, inputFile + ".tmp");

    string bit_size_mode = "UINT64_t";

    std::cout << "WT using memory = " << sdsl::size_in_bytes(wt) / 1000 << "[KB]" << std::endl;

    char MODE = mode[0];
    //char FPOSMODE = mode[1];

    //std::cout << "LPOS Data Structure: " << MODE << std::endl;
    //std::cout << "FPOS Data Structure: " << FPOSMODE << std::endl;

    uint64_t ms_count = 0;
    stool::lcp_on_rlbwt::STTreeAnalysisResult st_result;

    double construction_time = 0;
    std::chrono::system_clock::time_point mid;
    if (MODE == '0')
    {
        std::vector<uint64_t> lpos_vec;
        run_bits.to_vector(lpos_vec);

        assert(diff_char_vec.size() + 1 == lpos_vec.size());
        using LPOSDS = stool::lcp_on_rlbwt::RankSupportVectorWrapper<std::vector<uint64_t>>;
        LPOSDS lpos_vec_wrapper(lpos_vec);

        using FPOSDS = std::vector<uint64_t>;
        using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<INDEX, LPOSDS, FPOSDS>;

        FPOSDS fposds = stool::lcp_on_rlbwt::FPosDataStructure::construct(diff_char_vec, lpos_vec_wrapper);

        RDS ds = RDS(diff_char_vec, wt, lpos_vec_wrapper, fposds);

        std::cout << "Enumerate Maximal Substrings..." << std::endl;
        stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds);
        ms_count = stool::lcp_on_rlbwt::Application<RDS>::outputMaximalSubstrings(out, stnodeTraverser, st_result);
    }
    else
    {
        stool::EliasFanoVector lpos_vec;
        lpos_vec.build_from_builder(run_bits);
        std::cout << "LPOS Vec using memory = " << lpos_vec.get_using_memory() / 1000 << "[KB]" << std::endl;

        //lpos_vec.build_from_bit_vector(run_bits);
        using LPOSDS = stool::EliasFanoVector;
        using FPOSDS = stool::lcp_on_rlbwt::LightFPosDataStructure;
        FPOSDS fposds = stool::lcp_on_rlbwt::LightFPosDataStructure(diff_char_vec, lpos_vec, wt);
        std::cout << "FPOS Vec using memory = " << fposds.get_using_memory() / 1000 << "[KB]" << std::endl;

        mid = std::chrono::system_clock::now();
        construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start).count();
        std::cout << "Construction time: " << construction_time << "[ms]" << std::endl;

        if (bwtAnalysis.str_size < UINT32_MAX - 10)
        {
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<uint32_t, LPOSDS, FPOSDS>;
            RDS ds = RDS(diff_char_vec, wt, lpos_vec, fposds);

            std::cout << "Enumerate Maximal Substrings..." << std::endl;
            stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<uint32_t, RDS> stnodeTraverser;
            stnodeTraverser.initialize(thread_num, ds);
            ms_count = stool::lcp_on_rlbwt::Application<RDS>::outputMaximalSubstrings(out, stnodeTraverser, st_result);
            bit_size_mode = "UINT32_t";
        }
        else
        {
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<uint64_t, LPOSDS, FPOSDS>;
            RDS ds = RDS(diff_char_vec, wt, lpos_vec, fposds);

            std::cout << "Enumerate Maximal Substrings..." << std::endl;
            stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<uint64_t, RDS> stnodeTraverser;
            stnodeTraverser.initialize(thread_num, ds);
            ms_count = stool::lcp_on_rlbwt::Application<RDS>::outputMaximalSubstrings(out, stnodeTraverser, st_result);
        }
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double enumeration_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - mid).count();

    double bps = ((double)bwtAnalysis.str_size / ((double)elapsed / 1000)) / 1000;

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "RLBWT File \t\t\t\t : " << inputFile << std::endl;
    std::cout << "Output \t\t\t\t\t : " << outputFile << std::endl;
    std::cout << "LPOS and FPos Vector type \t\t\t\t : " << (MODE == '0' ? "std::vector<uint64_t>" : "EliasFano") << std::endl;
    std::cout << "Peak children count \t\t\t : " << st_result.max_nodes_at_level << std::endl;

    std::cout << "Thread number \t\t\t\t : " << thread_num << std::endl;
    std::cout << "Integer Type \t\t\t\t : " << bit_size_mode << std::endl;

    std::cout << "The length of the input text \t\t : " << bwtAnalysis.str_size << std::endl;
    std::cout << "The number of runs on BWT \t\t : " << bwtAnalysis.run_count << std::endl;
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

    p.add<bool>("debug", 'd', "debug", false, false);

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    bool debug = p.get<bool>("debug");
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

    if (!debug)
    {
        computeMaximalSubstrings(inputFile, outputFile, mode, thread_num);
    }
    else
    {
        testMaximalSubstrings<INDEX>(inputFile, mode, thread_num);

        //test2(inputFile, mode);
    }
}
