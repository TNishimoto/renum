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

#include "stool/src/io.hpp"
#include "stool/src/cmdline.h"
#include "stool/src/debug.hpp"

#include <sdsl/bit_vectors.hpp>

#include "../debug/beller_debug.hpp"
#include "../debug/naive_algorithms.hpp"
#include "../debug/stnode_checker.hpp"

#include "../stnode_enumerator/application.hpp"
#include "../rlbwt/bwt_analysis_result.hpp"
#include "../rlbwt/rle.hpp"

#include "../rlbwt/fpos_data_structure.hpp"
#include "../rlbwt/bwt_decompress.hpp"

#include "../stnode_enumerator/explicit_weiner_link_computer.hpp"

#include "../stnode_enumerator/application.hpp"

#include "../stnode_enumerator/depth_first_search/dfs_traverser.hpp"

using namespace std;
using namespace stool;

using CHAR = char;
using INDEX = uint64_t;
using LCPINTV = stool::LCPInterval<uint64_t>;

class LCPIntervalTest
{
public:
    template <typename STNODES>
    static std::vector<stool::LCPInterval<uint64_t>> testLCPIntervals(STNODES &stnodeSequencer, stool::stnode_on_rlbwt::STNodeChecker *checker)
    {

        std::vector<stool::LCPInterval<uint64_t>> r;

        auto it = stnodeSequencer.begin();
        while (it != stnodeSequencer.end())
        {
            std::vector<stool::LCPInterval<uint64_t>> r2;
            std::cout << "Test, LCP = " << it.get_depth() << std::endl;
            for (auto node_it = it.begin(); node_it != it.end(); node_it++)
            {
                if (checker != nullptr && checker->bwt.size() < 500)
                {
                    node_it.print();
                }
                if (stnodeSequencer.has_edge_characters())
                {

                    uint64_t child_count = node_it.get_children_count();
                    uint64_t depth = it.get_depth();
                    for (uint64_t i = 0; i < child_count; i++)
                    {
                        uint64_t left = node_it.get_child_left_boundary(i);
                        uint64_t right = node_it.get_child_right_boundary(i);
                        char c = node_it.get_edge_character(i);
                        if (checker != nullptr)
                        {
                            checker->check_edge_character(left, right, c, depth);
                        }
                    }
                }

                /*
                */

                stool::LCPInterval<uint64_t> copy;
                uint64_t left = stnodeSequencer.get_left(node_it);
                uint64_t right = stnodeSequencer.get_right(node_it);

                copy.i = left;
                copy.j = right;
                copy.lcp = stnodeSequencer.get_current_lcp();
                r2.push_back(copy);

                if (checker != nullptr)
                {
                    bool b = stnodeSequencer.check_maximal_repeat(node_it);
                    checker->check_maximal_repeat(left, right, b);
                }
            }
            //std::cout << std::endl;
            if (checker != nullptr)
            {
                checker->increment(r2);
            }
            for (auto &it : r2)
            {
                r.push_back(it);
            }
            it++;
        }
        std::cout
            << "STOP" << std::endl;
        return r;
    }

    template <typename INDEX_SIZE, typename INTERVAL_SEARCH, typename CHAR = uint8_t>
    static std::vector<stool::LCPInterval<uint64_t>> testLCPIntervals(stool::stnode_on_rlbwt::DFSTraverser<INDEX_SIZE, INTERVAL_SEARCH, CHAR> &stnodeSequencer, stool::stnode_on_rlbwt::STNodeChecker *checker)
    {

        std::vector<stool::LCPInterval<uint64_t>> r;

        auto it = stnodeSequencer.begin();
        while (it != stnodeSequencer.end())
        {
            uint64_t left = it.get_left();
            uint64_t right = it.get_right();
            uint64_t lcp = it.get_lcp();
            auto intv = stool::LCPInterval<uint64_t>(left, right, lcp);
            if(checker != nullptr){
                checker->check_lcp_interval(left, right, lcp);
            }
            if (stnodeSequencer.get_input_text_length() < 100)
            {
                std::cout << intv.to_string() << std::endl;
            }
            if(intv.lcp > stnodeSequencer.get_input_text_length()){
                assert(false);
                throw -1;
            }
            r.push_back(intv);

            it++;
        }
        std::cout
            << "STOP" << std::endl;
        return r;
    }
};
template <typename INDEX>
void testLCPIntervals(std::string inputFile, string mode, int thread_num)
{
    stool::rlbwt2::BWTAnalysisResult analysisResult;
    stool::stnode_on_rlbwt::RLE<uint8_t> rlbwt;
    rlbwt.load(inputFile, analysisResult);

    stool::stnode_on_rlbwt::STNodeChecker stnc;
    stnc.initialize(inputFile);

    std::vector<stool::LCPInterval<uint64_t>> test_Intervals;
    using RDS = stool::stnode_on_rlbwt::RLEWaveletTree<INDEX>;
    RDS ds = RDS(&rlbwt);

    if (mode == "A")
    {

        std::cout << "General Test" << std::endl;
        stool::stnode_on_rlbwt::SuffixTreeNodes<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds, true);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest::testLCPIntervals(stnodeTraverser, &stnc);
        test_Intervals.swap(tmp);
    }
    else if (mode == "B")
    {

        std::cout << "Standard Test" << std::endl;
        stool::stnode_on_rlbwt::STNodeTraverser<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds, true);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest::testLCPIntervals(stnodeTraverser, &stnc);
        test_Intervals.swap(tmp);
    }
    else if (mode == "C")
    {

        std::cout << "Fast Test" << std::endl;
        stool::stnode_on_rlbwt::FastSTNodeTraverser<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds, true);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest::testLCPIntervals(stnodeTraverser, &stnc);
        test_Intervals.swap(tmp);
    }
    else
    {

        std::cout << "Single Test" << std::endl;
        using INTERVAL_SEARCH = stool::stnode_on_rlbwt::ExplicitWeinerLinkComputerOnRLBWT<RDS>;
        stool::stnode_on_rlbwt::SingleSTNodeTraverser<INDEX, INTERVAL_SEARCH> stnodeTraverser;
        INTERVAL_SEARCH em;
        em.initialize(&ds);
        stnodeTraverser.initialize(&em, true);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest::testLCPIntervals(stnodeTraverser, &stnc);
        test_Intervals.swap(tmp);
    }

    //std::vector<stool::LCPInterval<uint64_t>> tmp = stool::stnode_on_rlbwt::Application<RDS>::testLCPIntervals(stnodeTraverser);

    stool::beller::equal_check_lcp_intervals(test_Intervals, stnc.lcp_intervals);
    std::cout << "LCP interval check OK!" << std::endl;
}
void computeLCPIntervals_beller(std::string inputFile)
{

    //string text = "";
    std::cout << "Loading : " << inputFile << std::endl;

    stool::stnode_on_rlbwt::STNodeChecker stnc;
    stnc.initialize(inputFile);

    std::vector<uint8_t> bwt2;
    stool::bwt::load(inputFile, bwt2);
    sdsl::int_vector<> bwt;
    bwt.width(8);
    bwt.resize(bwt2.size());
    for (uint64_t i = 0; i < bwt2.size(); i++)
    {
        bwt[i] = bwt2[i];
    }
    /*
    sdsl::bit_vector bv;
    bv.resize(bwt.size());
    bool b = true;
    for (uint64_t i = 0; i < bwt.size(); i++)
    {
        if (i > 0 && bwt[i] != bwt[i - 1])
        {
            b = !b;
        }
        bv[i] = b;
        //std::cout << (uint64_t)bv[i];
    }
    std::cout << std::endl;
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
    traverser.initialize(&wsearch, true);

    auto test_Intervals = LCPIntervalTest::testLCPIntervals(traverser, &stnc);
}

void computeLCPIntervals_DFS(std::string inputFile)
{
    std::cout << "DFS Test" << std::endl;
    //string text = "";
    std::cout << "Loading : " << inputFile << std::endl;

    stool::stnode_on_rlbwt::STNodeChecker stnc;
    stnc.initialize(inputFile);

    std::vector<uint8_t> bwt2;
    stool::bwt::load(inputFile, bwt2);
    sdsl::int_vector<> bwt;
    bwt.width(8);
    bwt.resize(bwt2.size());
    for (uint64_t i = 0; i < bwt2.size(); i++)
    {
        bwt[i] = bwt2[i];
    }
    /*
    sdsl::bit_vector bv;
    bv.resize(bwt.size());
    bool b = true;
    for (uint64_t i = 0; i < bwt.size(); i++)
    {
        if (i > 0 && bwt[i] != bwt[i - 1])
        {
            b = !b;
        }
        bv[i] = b;
        //std::cout << (uint64_t)bv[i];
    }
    std::cout << std::endl;
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
    stool::stnode_on_rlbwt::DFSTraverser<uint32_t, stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<uint32_t>> traverser;
    traverser.initialize(&wsearch, true);

    auto test_Intervals = LCPIntervalTest::testLCPIntervals<uint32_t, stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<uint32_t>, uint8_t>(traverser, &stnc);
    stool::beller::equal_check_lcp_intervals(test_Intervals, stnc.lcp_intervals);
    std::cout << "OK!" << std::endl;
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
    p.add<string>("mode", 'm', "mode", false, "A");
    p.add<int>("thread_num", 'p', "thread number", false, -1);

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    string mode = p.get<string>("mode");

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
    if (mode == "beller")
    {
        computeLCPIntervals_beller(inputFile);
    }
    else if (mode == "dfs")
    {
        computeLCPIntervals_DFS(inputFile);
    }
    else
    {
        testLCPIntervals<INDEX>(inputFile, mode, thread_num);
    }
}
