#include <cassert>
#include <chrono>
#include "stool/include/io.hpp"
#include "stool/include/sa_bwt_lcp.hpp"

#include "stool/include/print.hpp"
#include "stool/include/cmdline.h"
#include "stool/include/debug.hpp"
#include "../module/libdivsufsort/sa.hpp"
//#include "../module/rlbwt_iterator/src/include/rlbwt_iterator.hpp"
//#include "module/rlbwt_iterator/src/include/bwt.hpp"


#include <sdsl/bit_vectors.hpp>

#include "../debug/beller_debug.hpp"
#include "../debug/naive_algorithms.hpp"
#include "../debug/stnode_checker.hpp"

#include "../stnode_enumerator/application.hpp"
#include "../stnode_enumerator/rlcp_interval_enumerator.hpp"
#include "../rlbwt/io.hpp"
#include "../rlbwt/fpos_data_structure.hpp"
#include "../rlbwt/bwt_decompress.hpp"


using namespace std;
using namespace stool;

using CHAR = char;
using INDEX = uint64_t;
using LCPINTV = stool::LCPInterval<uint64_t>;

void computeMaximalSubstrings(std::string inputFile, std::string dataFile, string mode, int thread_num, bool is_count)
{
    auto start = std::chrono::system_clock::now();

    sdsl::int_vector<> diff_char_vec;
    stool::EliasFanoVectorBuilder run_bits;
    //std::vector<bool> run_bits;
    auto bwtAnalysis = stool::rlbwt2::load_RLBWT_from_file(inputFile, diff_char_vec, run_bits);
    std::cout << "BWT using memory = " << sdsl::size_in_bytes(diff_char_vec) / 1000 << "[KB]" << std::endl;
    std::cout << "Run bits using memory = " << run_bits.get_using_memory() / 1000 << "[KB]" << std::endl;

    uint64_t data_structure_bytes = 0;
    sdsl::store_to_file(diff_char_vec, inputFile + ".tmp");

    stool::WT wt;
    construct(wt, inputFile + ".tmp");

    string bit_size_mode = "UINT64_t";
    data_structure_bytes += sdsl::size_in_bytes(wt);
    std::cout << "WT using memory = " << sdsl::size_in_bytes(wt) / 1000 << "[KB]" << std::endl;

    //char MODE = mode[0];
    //char FPOSMODE = mode[1];

    //std::cout << "LPOS Data Structure: " << MODE << std::endl;
    //std::cout << "FPOS Data Structure: " << FPOSMODE << std::endl;

    //uint64_t ms_count = 0;
    //stool::renum::STTreeAnalysisResult st_result;

    double construction_time = 0;
    std::chrono::system_clock::time_point mid;

    stool::EliasFanoVector lpos_vec;
    lpos_vec.build_from_builder(run_bits);
    std::cout << "LPOS Vec using memory = " << lpos_vec.get_using_memory() / 1000 << "[KB]" << std::endl;
    data_structure_bytes += lpos_vec.get_using_memory();

    //lpos_vec.build_from_bit_vector(run_bits);
    using LPOSDS = stool::EliasFanoVector;
    using FPOSDS = stool::renum::LightFPosDataStructure;
    FPOSDS fposds = stool::renum::LightFPosDataStructure(diff_char_vec, lpos_vec, wt);
    std::cout << "FPOS Vec using memory = " << fposds.get_using_memory() / 1000 << "[KB]" << std::endl;
    data_structure_bytes += fposds.get_using_memory();

    mid = std::chrono::system_clock::now();
    construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start).count();
    std::cout << "Construction time: " << construction_time << "[ms]" << std::endl;
    std::cout << "Data structure Size \t\t\t : " << (data_structure_bytes / 1000) << "[KB]" << std::endl;

    using RDS = stool::renum::RLEWaveletTree<uint32_t, LPOSDS, FPOSDS>;
    RDS ds = RDS(diff_char_vec, wt, lpos_vec, fposds);

    std::cout << "Enumerate Maximal Substrings..." << std::endl;
    stool::renum::STNodeTraverser<uint32_t, RDS> stnodeTraverser;
    stnodeTraverser.initialize(thread_num, ds);

    if (dataFile.size() != 0)
    {
        std::ifstream inp;
        inp.open(dataFile, std::ios::binary);
        stnodeTraverser.load(inp);
        inp.close();
    }

    uint64_t input_children_count = stnodeTraverser.child_count();
    uint64_t output_children_count = 0;
    if (is_count)
    {
        output_children_count = stnodeTraverser.parallel_count();
    }
    else
    {
        stnodeTraverser.process();
        output_children_count = stnodeTraverser.child_count();
    }

    uint64_t current_lcp = stnodeTraverser.get_current_lcp();
    std::string outputFile = inputFile + "." + std::to_string(current_lcp) + ".data";
    if (!is_count)
    {
        std::ofstream out(outputFile, std::ios::out | std::ios::binary);
        if (!out)
        {
            throw std::runtime_error("Cannot open the output file!");
        }

        stnodeTraverser.write(out);
    }
    bit_size_mode = "UINT32_t";

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double enumeration_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - mid).count();

    double bps = ((double)bwtAnalysis.str_size / ((double)elapsed / 1000)) / 1000;

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "RLBWT File \t\t\t\t : " << inputFile << std::endl;
    if(!is_count){
        std::cout << "Output \t\t\t\t\t : " << outputFile << std::endl;
    }
    std::cout << "Input chilren count \t\t\t : " << input_children_count << std::endl;
    std::cout << "output chilren count \t\t\t : " << output_children_count << std::endl;

    std::cout << "Data structure Size \t\t\t : " << (data_structure_bytes / 1000) << "[KB]" << std::endl;

    std::cout << "Thread number \t\t\t\t : " << thread_num << std::endl;
    std::cout << "Integer Type \t\t\t\t : " << bit_size_mode << std::endl;

    std::cout << "The length of the input text \t\t : " << bwtAnalysis.str_size << std::endl;
    std::cout << "The number of runs on BWT \t\t : " << bwtAnalysis.run_count << std::endl;
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
    p.add<string>("data_file", 'x', "data file name", false, "");

    //p.add<bool>("mode", 'm', "mode", false, false);
    //p.add<string>("output_file", 'o', "output file name", false, "");
    p.add<string>("mode", 'm', "mode", false, "1");
    p.add<int>("thread_num", 'p', "thread number", false, -1);

    p.add<bool>("debug", 'd', "debug", false, false);
    p.add<bool>("count", 'c', "count", false, false);

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    string dataFile = p.get<string>("data_file");

    //bool debug = p.get<bool>("debug");
    bool is_count = p.get<bool>("count");

    string mode = p.get<string>("mode");

    //string outputFile = p.get<string>("output_file");
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

    computeMaximalSubstrings(inputFile, dataFile, mode, thread_num, is_count);
}
