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

using namespace std;
using namespace stool;


using CHAR = char;
using INDEX = uint64_t;
using LCPINTV = stool::LCPInterval<uint64_t>;

/*
void testLCPIntervals(std::string inputFile, bool lightWeight, bool correctCheck)
{
    stool::rlbwt::RLBWT<std::vector<CHAR>, std::vector<INDEX>> rlestr = stool::rlbwt::Constructor::load_RLBWT_from_file<CHAR, INDEX>(inputFile);
    using FPOSDS = std::vector<uint64_t>;
    using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<RLBWT<>, uint64_t, FPOSDS>;
    FPOSDS fposds = stool::lcp_on_rlbwt::FPosDataStructure::construct(rlestr);
    RDS ds = RDS(rlestr, fposds);

    auto testVec = stool::lcp_on_rlbwt::Application<RLBWT<>, RDS>::constructLCPIntervals(rlestr, &ds);

    if (correctCheck)
    {

        std::vector<char> text2 = stool::load_char_vec_from_file(inputFile, true);
        vector<INDEX> sa = stool::construct_suffix_array(text2);
        auto correctLCP = stool::constructLCP(text2, sa);
        std::cout << "Correct" << std::endl;
        std::vector<LCPINTV> correct_intervals = stool::beller::naive_compute_complete_lcp_intervals<uint64_t>(sa, correctLCP);
        stool::beller::equal_check_lcp_intervals(testVec, correct_intervals);
        std::cout << "OK!" << std::endl;
    }
    std::cout << "rlbwt = " << rlestr.rle_size() << std::endl;
}
*/

/*
void test2(std::string inputFile, bool lightWeight)
{

    using RLBWT_STR = stool::rlbwt::RLBWT<std::vector<CHAR>, std::vector<INDEX>>;
    RLBWT_STR rlestr = stool::rlbwt::Constructor::load_RLBWT_from_file<CHAR, INDEX>(inputFile);
    lcp_on_rlbwt::RLBWTDataStructures<RLBWT_STR, uint64_t> rlbwtds(rlestr, lightWeight);
    lcp_on_rlbwt::RLCPIntervalEnumerator<RLBWT_STR, uint64_t, RDS> enumerator(&rlbwtds);

    for(auto w : enumerator){
        w.first.print();
    }
}
*/

void testMaximalSubstrings(std::string inputFile, string mode, int thread_num)
{

    sdsl::int_vector<> diff_char_vec;
    stool::EliasFanoVectorBuilder run_bits;
    auto bwtAnalysis =  stool::rlbwt2::load_RLBWT_from_file(inputFile, diff_char_vec, run_bits);
    stool::WT wt;
    construct_im(wt, diff_char_vec);

    /*
    sdsl::wt_huff_int<> wt2;
    construct_im(wt2, diff_char_vec);
    */


//DEBUG
    if(diff_char_vec.size() < 100){
        std::cout << "Run heads: ";
        for(uint64_t i=0;i<diff_char_vec.size();i++){
            std::cout << (char)bwtAnalysis.id_to_character_vec[diff_char_vec[i]];
        }
        std::cout << std::endl;
        /*
        for(uint64_t i=0;i<wt2.size();i++){
            std::cout << (char)bwtAnalysis.id_to_character_vec[wt2[i]];
        }
        std::cout << std::endl;
        */
    }

    stool::lcp_on_rlbwt::STNodeChecker stnc;
    stnc.initialize(inputFile);
    /*
    std::vector<char> text = stool::bwt::decompress_bwt(inputFile);
    vector<uint64_t> sa = stool::construct_suffix_array(text);    
    vector<stool::LCPInterval<uint64_t>> correct_intervals = stool::esaxx::naive_compute_maximal_substrings<char, uint64_t>(text, sa);
    vector<stool::LCPInterval<uint64_t>> correct_lcp_intervals = stool::esaxx::naive_compute_complete_lcp_intervals<char, uint64_t>(text, sa);
    */


    //uint64_t input_text_size = ds.str_size();
    std::vector<stool::LCPInterval<uint64_t>> test_Intervals;

    char LPOSMODE = mode[0];
    char FPOSMODE = mode[1];

    std::cout << "LPOS Data Structure: " << LPOSMODE << std::endl;
    std::cout << "FPOS Data Structure: " << FPOSMODE << std::endl;

    if (LPOSMODE == '0')
    {
        std::vector<uint64_t> lpos_vec;
        std::cout << "Building lpos_vec ..." << std::flush;
        run_bits.to_vector(lpos_vec);
        std::cout << "[Finished]" << std::endl;

        assert(diff_char_vec.size() + 1 == lpos_vec.size());
        using LPOSDS = stool::lcp_on_rlbwt::RankSupportVectorWrapper<std::vector<uint64_t>>;
        LPOSDS lpos_vec_wrapper(lpos_vec);

        if (FPOSMODE == '0')
        {
            using FPOSDS = std::vector<uint64_t>;
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<uint64_t, LPOSDS, FPOSDS>;
            FPOSDS fposds = stool::lcp_on_rlbwt::FPosDataStructure::construct(diff_char_vec, lpos_vec_wrapper);

            RDS ds = RDS(diff_char_vec, wt, lpos_vec_wrapper, fposds);
            ds.stnc = &stnc;
            stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<INDEX, RDS> stnodeTraverser;

            stnodeTraverser.initialize(thread_num, ds);
            std::vector<stool::LCPInterval<uint64_t>> tmp = stool::lcp_on_rlbwt::Application<RDS>::testLCPIntervals(stnodeTraverser);
            test_Intervals.swap(tmp);
        }
        else
        {
            using FPOSDS = stool::lcp_on_rlbwt::LightFPosDataStructure;
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<uint64_t, LPOSDS, FPOSDS>;
            FPOSDS fposds = stool::lcp_on_rlbwt::LightFPosDataStructure(diff_char_vec, lpos_vec_wrapper, wt);
            RDS ds = RDS(diff_char_vec, wt, lpos_vec_wrapper, fposds);
            ds.stnc = &stnc;

            stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<INDEX, RDS> stnodeTraverser;
            stnodeTraverser.initialize(thread_num, ds);
            std::vector<stool::LCPInterval<uint64_t>> tmp = stool::lcp_on_rlbwt::Application<RDS>::testLCPIntervals(stnodeTraverser);
            test_Intervals.swap(tmp);
        }
    }
    else
    {

        stool::EliasFanoVector lpos_vec;
        lpos_vec.build_from_builder(run_bits);
        //lpos_vec.build_from_bit_vector(run_bits);
        using LPOSDS = stool::EliasFanoVector;



        if (FPOSMODE == '0')
        {
            using FPOSDS = std::vector<uint64_t>;
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<uint64_t, LPOSDS, FPOSDS>;
            FPOSDS fposds = stool::lcp_on_rlbwt::FPosDataStructure::construct(diff_char_vec, lpos_vec);
            RDS ds = RDS(diff_char_vec, wt, lpos_vec, fposds);
                        ds.stnc = &stnc;

            stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<INDEX, RDS> stnodeTraverser;
            stnodeTraverser.initialize(thread_num, ds);
            std::vector<stool::LCPInterval<uint64_t>> tmp = stool::lcp_on_rlbwt::Application<RDS>::testLCPIntervals(stnodeTraverser);
            test_Intervals.swap(tmp);
        }
        else
        {
            using FPOSDS = stool::lcp_on_rlbwt::LightFPosDataStructure;
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<uint64_t, LPOSDS, FPOSDS>;
            FPOSDS fposds = stool::lcp_on_rlbwt::LightFPosDataStructure(diff_char_vec, lpos_vec, wt);
            RDS ds = RDS(diff_char_vec, wt, lpos_vec, fposds);
                        ds.stnc = &stnc;

            stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<INDEX, RDS> stnodeTraverser;
            stnodeTraverser.initialize(thread_num, ds);
            std::vector<stool::LCPInterval<uint64_t>> tmp = stool::lcp_on_rlbwt::Application<RDS>::testLCPIntervals(stnodeTraverser);
            test_Intervals.swap(tmp);
        }
    }


   
    
    stool::beller::equal_check_lcp_intervals(test_Intervals, stnc.lcp_intervals);
    std::cout << "Maximal repeats check OK!" << std::endl;
    

    
}

void computeMaximalSubstrings(std::string inputFile,std::string outputFile, string mode, int thread_num)
{
    auto start = std::chrono::system_clock::now();


    sdsl::int_vector<> diff_char_vec;
    stool::EliasFanoVectorBuilder run_bits;
    //std::vector<bool> run_bits;
    auto bwtAnalysis = stool::rlbwt2::load_RLBWT_from_file(inputFile, diff_char_vec, run_bits);


    stool::WT wt;
    construct_im(wt, diff_char_vec);
    /*
    throw -1;

    sdsl::wt_huff_int<> wt2;
    construct_im(wt2, diff_char_vec);
    */




    //DEBUG
    if(diff_char_vec.size() < 100){
        std::cout << "Run heads: ";
        for(uint64_t i=0;i<diff_char_vec.size();i++){
            std::cout << (char)diff_char_vec[i];
        }
        std::cout << std::endl;
    }

    //uint64_t input_text_size = ds.str_size();
    std::vector<stool::LCPInterval<uint64_t>> test_Intervals;

    char LPOSMODE = mode[0];
    char FPOSMODE = mode[1];

    std::cout << "LPOS Data Structure: " << LPOSMODE << std::endl;
    std::cout << "FPOS Data Structure: " << FPOSMODE << std::endl;

    std::ofstream out(outputFile, std::ios::out | std::ios::binary);
    if (!out)
    {
        throw std::runtime_error("Cannot open the output file!");
    }
    uint64_t ms_count = 0;
    stool::lcp_on_rlbwt::STTreeAnalysisResult st_result;

    if (LPOSMODE == '0')
    {
        std::vector<uint64_t> lpos_vec;
        run_bits.to_vector(lpos_vec);

        assert(diff_char_vec.size() + 1 == lpos_vec.size());
        using LPOSDS = stool::lcp_on_rlbwt::RankSupportVectorWrapper<std::vector<uint64_t>>;
        LPOSDS lpos_vec_wrapper(lpos_vec);

        //using LPOSDS = std::vector<uint64_t>;

        if (FPOSMODE == '0')
        {
            using FPOSDS = std::vector<uint64_t>;
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<uint64_t, LPOSDS, FPOSDS>;

            FPOSDS fposds = stool::lcp_on_rlbwt::FPosDataStructure::construct(diff_char_vec, lpos_vec_wrapper);

            RDS ds = RDS(diff_char_vec, wt, lpos_vec_wrapper, fposds);

            std::cout << "Enumerate Maximal Substrings..." << std::endl;
            stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<INDEX, RDS> stnodeTraverser;
            stnodeTraverser.initialize(thread_num, ds);
            ms_count = stool::lcp_on_rlbwt::Application<RDS>::outputMaximalSubstrings(out, stnodeTraverser, st_result);
        }
        else
        {
            using FPOSDS = stool::lcp_on_rlbwt::LightFPosDataStructure;
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<uint64_t, LPOSDS, FPOSDS>;
            FPOSDS fposds = stool::lcp_on_rlbwt::LightFPosDataStructure(diff_char_vec, lpos_vec_wrapper, wt);
            RDS ds = RDS(diff_char_vec, wt, lpos_vec_wrapper, fposds);

            std::cout << "Enumerate Maximal Substrings..." << std::endl;
            stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<INDEX, RDS> stnodeTraverser;
            stnodeTraverser.initialize(thread_num, ds);
            ms_count = stool::lcp_on_rlbwt::Application<RDS>::outputMaximalSubstrings(out, stnodeTraverser, st_result);
        }
    }
    else
    {
        stool::EliasFanoVector lpos_vec;
        lpos_vec.build_from_builder(run_bits);
        
        //lpos_vec.build_from_bit_vector(run_bits);
        using LPOSDS = stool::EliasFanoVector;
        if (FPOSMODE == '0')
        {
            using FPOSDS = std::vector<uint64_t>;
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<uint64_t, LPOSDS, FPOSDS>;
            FPOSDS fposds = stool::lcp_on_rlbwt::FPosDataStructure::construct(diff_char_vec, lpos_vec);
            RDS ds = RDS(diff_char_vec, wt, lpos_vec, fposds);

            std::cout << "Enumerate Maximal Substrings..." << std::endl;
            stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<INDEX, RDS> stnodeTraverser;
            stnodeTraverser.initialize(thread_num, ds);
            ms_count = stool::lcp_on_rlbwt::Application<RDS>::outputMaximalSubstrings(out, stnodeTraverser, st_result);
        }
        else
        {
            using FPOSDS = stool::lcp_on_rlbwt::LightFPosDataStructure;
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<uint64_t, LPOSDS, FPOSDS>;

            FPOSDS fposds = stool::lcp_on_rlbwt::LightFPosDataStructure(diff_char_vec, lpos_vec, wt);
            RDS ds = RDS(diff_char_vec, wt, lpos_vec, fposds);

            std::cout << "Enumerate Maximal Substrings..." << std::endl;
            stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<INDEX, RDS> stnodeTraverser;
            
            stnodeTraverser.initialize(thread_num, ds);

            
            ms_count = stool::lcp_on_rlbwt::Application<RDS>::outputMaximalSubstrings(out, stnodeTraverser, st_result);
            
        }
    }

    /*
    std::cout << "Enumerate Maximal Substrings..." << std::endl;
    uint64_t ms_count = 0;
    if (input_text_size - 10 < UINT32_MAX)
    {
        ms_count = stool::lcp_on_rlbwt::Application<RDS>::outputMaximalSubstrings(out, &ds);
    }
    else
    {
        ms_count = stool::lcp_on_rlbwt::Application<RDS>::outputMaximalSubstrings(out, &ds);
    }
    */
    auto end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    double bps = ((double)bwtAnalysis.str_size / ((double)elapsed / 1000)) / 1000;

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "RLBWT File \t\t\t\t : " << inputFile << std::endl;
    std::cout << "Output \t\t\t\t\t : " << outputFile << std::endl;
    std::cout << "LPOS Vector \t\t\t\t : " << (LPOSMODE == '0' ? "std::vector<uint64_t>" : "EliasFano") << std::endl;
    std::cout << "FPOS Vector \t\t\t\t : " << (FPOSMODE == '0' ? "std::vector<uint64_t>" : "EliasFano") << std::endl;
    std::cout << "Peak children count \t\t\t : " << st_result.max_nodes_at_level << std::endl;

    std::cout << "Thread number \t\t\t\t : " << thread_num << std::endl;


    std::cout << "The length of the input text \t\t : " << bwtAnalysis.str_size << std::endl;
    std::cout << "The number of runs on BWT \t\t : " << bwtAnalysis.run_count << std::endl;
    std::cout << "The number of maximum substrings \t : " << ms_count << std::endl;
    std::cout << "Excecution time \t\t\t : " << elapsed << "[ms]" << std::endl;
    std::cout << "Character per second \t\t\t : " << bps << "[KB/s]" << std::endl;

    std::cout << "_______________________________________________________" << std::endl;
    std::cout << "\033[39m" << std::endl;
}

int main(int argc, char *argv[])
{
    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file name", true);
    //p.add<bool>("mode", 'm', "mode", false, false);
    p.add<string>("output_file", 'o', "output file name", false, "");
    p.add<string>("mode", 'm', "mode", false, "00");
    p.add<int>("thread_num", 'p', "thread number", false, -1);

    p.add<bool>("debug", 'd', "debug", false, false);

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    bool debug = p.get<bool>("debug");
    string mode = p.get<string>("mode");

    string outputFile = p.get<string>("output_file");
    string format = "binary";
    int thread_num = p.get<int>("thread_num");
    if(thread_num < 0){
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
        testMaximalSubstrings(inputFile, mode, thread_num);

        //test2(inputFile, mode);
    }
}
