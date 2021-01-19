#include <cassert>
#include <chrono>
#include <stdio.h>
#include "../module/stool/src/io.hpp"
#include "../module/stool/src/sa_bwt_lcp.hpp"

#include "../module/stool/src/print.hpp"
#include "../module/stool/src/cmdline.h"
#include "../module/stool/src/debug.hpp"
#include "../module/libdivsufsort/sa.hpp"

//#include "hpp/bwt.hpp"
#include "../basic/interval_search_data_structure.hpp"
//#include "../beller/beller_interval.hpp"
#include "../debug/beller_debug.hpp"


#include "../debug/naive_algorithms.hpp"
#include "../stnode_enumerator/single/single_stnode_traverser.hpp"
#include "../stnode_enumerator/application.hpp"

#include <sdsl/wt_algorithm.hpp>
#include "../basic/sdsl_functions.hpp"

//#include "../postorder_maximal_substring_intervals.hpp"
//#include "../forward_bwt.hpp"

using namespace std;
//using namespace stool;
//using namespace stool::rlbwt;

using CHAR = char;
using INDEX = uint64_t;
using LCPINTV = stool::LCPInterval<uint64_t>;
class LCPIntervalTest
{
public:
    template <typename STNODES>
    static std::vector<stool::LCPInterval<uint64_t>> testLCPIntervals(STNODES &stnodeSequencer)
    {

        std::vector<stool::LCPInterval<uint64_t>> r;

        auto it = stnodeSequencer.begin();
        while (it != stnodeSequencer.end())
        {
            std::vector<stool::LCPInterval<uint64_t>> r2;

            for (auto interval : it)
            {
                stool::LCPInterval<uint64_t> copy;
                copy.i = interval.first;
                copy.j = interval.second;
                copy.lcp = stnodeSequencer.get_current_lcp();
                r2.push_back(copy);
            }
            for (auto &it : r2)
            {
                r.push_back(it);
            }
            it++;
        }

        std::cout << "STOP" << std::endl;
        return r;
    }
};

void computeLCPIntervals(std::string inputFile, bool correctCheck)
{

    //string text = "";
    std::cout << "Loading : " << inputFile << std::endl;
    std::vector<uint8_t> text = stool::load_text_from_file(inputFile, true);
    vector<INDEX> sa = stool::construct_suffix_array(text);
    sdsl::int_vector<> bwt;
    stool::FMIndex::constructBWT(text, sa, bwt);

    /*
    sdsl::bit_vector bv;
    bv.resize(bwt.size());
    bool b = true;
    //std::cout << bwt.size() << std::endl;
    for (uint64_t i = 0; i < bwt.size(); i++)
    {
        if (i > 0 && bwt[i] != bwt[i - 1])
        {
            b = !b;
        }
        bv[i] = b;
    }
    sdsl::bit_vector::rank_1_type bwt_bit_rank1(&bv);
    */

    std::vector<uint64_t> C;
    stool::FMIndex::constructC(bwt, C);

    wt_huff<> wt;
    construct_im(wt, bwt);

    uint64_t lastChar = bwt[bwt.size() - 1];

    stool::IntervalSearchDataStructure<uint8_t> range;
    range.initialize(&wt, &C, lastChar);

    stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<uint32_t> wsearch;
    wsearch.initialize(&range, bwt.size());
    stool::stnode_on_rlbwt::SingleSTNodeTraverser<uint32_t, stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<uint32_t>> traverser;
    traverser.initialize(&wsearch, false);
    auto test_Intervals = LCPIntervalTest::testLCPIntervals(traverser);

    //auto test_Intervals = stool::beller::computeLCPIntervals<uint64_t>(range);
    //test_Intervals.push_back(LCPINTV(0, text.size() - 1, 0));

    if (correctCheck)
    {
        auto correctLCP = stool::constructLCP(text, sa);
        std::cout << "Correct" << std::endl;
        std::vector<LCPINTV> correct_intervals = stool::esaxx::NaiveAlgorithms::naive_compute_lcp_intervals(text, sa);
        stool::beller::equal_check_lcp_intervals(test_Intervals, correct_intervals);
        std::cout << "OK!" << std::endl;
    }
}

void computeMaximalSubstrings(std::string inputFile, std::string outputFile)
{

    //string text = "";
    auto start = std::chrono::system_clock::now();
    //sdsl::bit_vector bv;

    std::cout << "Loading : " << inputFile << std::endl;
    std::cout << "Constructing dbit array for maximal repeats : " << std::endl;
    //stool::stnode_on_rlbwt::SDSLFunction::constructDBitArray(inputFile, bv);
    //sdsl::bit_vector::rank_1_type bwt_bit_rank1(&bv);


    std::cout << "Constructing Wavelet Tree..." << std::endl;
    wt_huff<> wt;
    std::vector<uint64_t> C;
    uint8_t lastChar = stool::stnode_on_rlbwt::SDSLFunction::load_wavelet_tree(inputFile, wt, C);



    uint64_t ms_count = 0;

    stool::IntervalSearchDataStructure<uint8_t> range;
    range.initialize(&wt, &C, lastChar);

    std::ofstream out(outputFile, std::ios::out | std::ios::binary);
    if (!out)
    {
        throw std::runtime_error("Cannot open the output file!");
    }
    uint64_t input_text_size = wt.size();

    auto mid = std::chrono::system_clock::now();

    std::cout << "Enumerating..." << std::endl;
            stool::stnode_on_rlbwt::STTreeAnalysisResult st_result;

    if (input_text_size - 10 < UINT32_MAX)
    {
        using INDEX_TYPE = uint32_t;

        stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<INDEX_TYPE> wsearch;
        wsearch.initialize(&range, input_text_size );
        stool::stnode_on_rlbwt::SingleSTNodeTraverser<INDEX_TYPE, stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<INDEX_TYPE>> traverser;
        traverser.initialize(&wsearch, false);
        ms_count = stool::stnode_on_rlbwt::Application::outputMaximalSubstrings(out, traverser, st_result);

    }
    else
    {
        using INDEX_TYPE = uint64_t;

        stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<INDEX_TYPE> wsearch;
        wsearch.initialize(&range, input_text_size );
        stool::stnode_on_rlbwt::SingleSTNodeTraverser<INDEX_TYPE, stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<INDEX_TYPE>> traverser;
        traverser.initialize(&wsearch, false);
        ms_count = stool::stnode_on_rlbwt::Application::outputMaximalSubstrings(out, traverser, st_result);

    }
    auto end = std::chrono::system_clock::now();
    double enumeration_time = std::chrono::duration_cast<std::chrono::seconds>(end - mid).count();
    double construction_time = std::chrono::duration_cast<std::chrono::seconds>(mid - start).count();
    //std::cout << "Construction time: " << construction_time << "[s]" << std::endl;

    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    double bps = ((double)input_text_size / ((double)elapsed)) / 1000;

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "BWT File \t\t\t\t : " << inputFile << std::endl;
    std::cout << "Output File \t\t\t\t : " << outputFile << std::endl;
    std::cout << "The length of the input text \t\t : " << input_text_size << std::endl;
    std::cout << "The number of maximum substrings \t : " << ms_count << std::endl;
    std::cout << "Peak count \t\t\t\t : " << st_result.max_nodes_at_level << std::endl;
    std::cout << "The memory usage of Wavelet tree \t : " << sdsl::size_in_bytes(wt) / 1000 << "[KB]" << std::endl;
    std::cout << "The memory usage of Queue \t\t : " << st_result.peak_memory_of_queue / 1000 << "[KB]" << std::endl;

    std::cout << "Excecution time \t\t\t : " << elapsed << "[s]" << std::endl;
    std::cout << "Character per second \t\t\t : " << bps << "[KB/s]" << std::endl;

    std::cout << "\t Preprocessing time \t\t : " << construction_time << "[s]" << std::endl;
    std::cout << "\t Enumeration time \t\t : " << enumeration_time << "[s]" << std::endl;

    std::cout << "_______________________________________________________" << std::endl;
    std::cout << "\033[39m" << std::endl;

}

int main(int argc, char *argv[])
{
    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file name", true);
    p.add<string>("output_file", 'o', "output file name", false, "");

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    string outputFile = p.get<string>("output_file");

    std::ifstream ifs(inputFile);
    bool inputFileExist = ifs.is_open();
    if (!inputFileExist)
    {
        std::cout << inputFile << " cannot open." << std::endl;
        return -1;
    }

    if (outputFile.size() == 0)
    {
            outputFile = inputFile + ".max";
    }
    computeMaximalSubstrings(inputFile, outputFile);


}
