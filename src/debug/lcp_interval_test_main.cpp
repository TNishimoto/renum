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
#include "../stnode_enumerator/weiner_link_search.hpp"

#include "../stnode_enumerator/application.hpp"

using namespace std;
using namespace stool;

using CHAR = char;
using INDEX = uint64_t;
using LCPINTV = stool::LCPInterval<uint64_t>;

template <typename RLBWTDS>
class LCPIntervalTest
{
public:
    using CHAR = typename RLBWTDS::CHAR;
    //using CHARVEC = typename RLBWT_STR::char_vec_type;
    using INDEX_SIZE = typename RLBWTDS::INDEX;
    using UCHAR = typename std::make_unsigned<CHAR>::type;

    template <typename STNODES, typename BWT>
    static std::vector<stool::LCPInterval<uint64_t>> testLCPIntervals(STNODES &stnodeSequencer, RLBWTDS *ds, BWT &bwt)
    {

        std::vector<stool::LCPInterval<uint64_t>> r;

        auto it = stnodeSequencer.begin();
        while (it != stnodeSequencer.end())
        {
            std::vector<stool::LCPInterval<uint64_t>> r2;

            for (auto node_it = it.begin(); node_it != it.end(); node_it++)
            {

                stool::LCPInterval<uint64_t> copy;
                INDEX_SIZE left = stnodeSequencer.get_left(node_it);
                INDEX_SIZE right = stnodeSequencer.get_right(node_it);

                copy.i = left;
                copy.j = right;
                copy.lcp = stnodeSequencer.get_current_lcp();
                r2.push_back(copy);
                bool b = stnodeSequencer.check_maximal_repeat(node_it);
                if (b)
                {
                    bool mb = false;
                    for (uint64_t i = left + 1; i <= right; i++)
                    {
                        uint8_t l = bwt[i - 1];
                        uint8_t r = bwt[i];

                        if (l != r)
                        {
                            mb = true;
                            break;
                        }
                    }
                    if (mb != b)
                    {
                        std::cout << mb << "/" << b << std::endl;
                        std::cout << "[" << left << ", " << right << "]" << std::endl;

                        std::cout << "BWT on Range: ";
                        for (uint64_t i = left; i <= right; i++)
                        {
                            std::cout << (char)bwt[i];
                        }
                        std::cout << std::endl;
                    }
                    assert(mb == b);
                }
            }
            if (ds != nullptr && ds->stnc != nullptr)
            {
                ds->stnc->increment(r2);
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
};
template <typename INDEX>
void testLCPIntervals(std::string inputFile, string mode, int thread_num)
{

    sdsl::int_vector<> diff_char_vec;
    stool::EliasFanoVectorBuilder run_bits;
    auto bwtAnalysis = stool::rlbwt2::load_RLBWT_from_file(inputFile, diff_char_vec, run_bits);
    stool::WT wt;
    construct_im(wt, diff_char_vec);

    std::cout << "BWT using memory = " << sdsl::size_in_bytes(diff_char_vec) / 1000 << "[KB]" << std::endl;
    std::cout << "Run bits using memory = " << run_bits.get_using_memory() / 1000 << "[KB]" << std::endl;

    //DEBUG
    if (diff_char_vec.size() < 100)
    {
        std::cout << "Run heads: ";
        for (uint64_t i = 0; i < diff_char_vec.size(); i++)
        {
            std::cout << (char)bwtAnalysis.id_to_character_vec[diff_char_vec[i]];
        }
        std::cout << std::endl;
    }

    stool::lcp_on_rlbwt::STNodeChecker stnc;
    stnc.initialize(inputFile);

    std::vector<stool::LCPInterval<uint64_t>> test_Intervals;

    stool::EliasFanoVector lpos_vec;
    lpos_vec.build_from_builder(run_bits);
    //lpos_vec.build_from_bit_vector(run_bits);
    using LPOSDS = stool::EliasFanoVector;

    using FPOSDS = stool::lcp_on_rlbwt::LightFPosDataStructure;
    using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<INDEX, LPOSDS, FPOSDS>;
    FPOSDS fposds = stool::lcp_on_rlbwt::LightFPosDataStructure(diff_char_vec, lpos_vec, wt);
    RDS ds = RDS(diff_char_vec, wt, lpos_vec, fposds);
    ds.stnc = &stnc;

    std::vector<uint8_t> plain_bwt;
    stool::bwt::load(inputFile, plain_bwt);

    if (mode == "A")
    {

        std::cout << "General Test" << std::endl;
        stool::lcp_on_rlbwt::SuffixTreeNodes<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest<RDS>::testLCPIntervals(stnodeTraverser, &ds, plain_bwt);
        test_Intervals.swap(tmp);
    }
    else if (mode == "B")
    {

        std::cout << "Standard Test" << std::endl;
        stool::lcp_on_rlbwt::STNodeTraverser<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest<RDS>::testLCPIntervals(stnodeTraverser, &ds, plain_bwt);
        test_Intervals.swap(tmp);
    }
    else if (mode == "C")
    {

        std::cout << "Fast Test" << std::endl;
        stool::lcp_on_rlbwt::FastSTNodeTraverser<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest<RDS>::testLCPIntervals(stnodeTraverser, &ds, plain_bwt);
        test_Intervals.swap(tmp);
    }
    else
    {

        std::cout << "Single Test" << std::endl;
        using INTERVAL_SEARCH = stool::lcp_on_rlbwt::ExplicitWeinerLinkEmulator<INDEX, RDS>;
        stool::lcp_on_rlbwt::SingleSTNodeTraverser<INDEX, INTERVAL_SEARCH> stnodeTraverser;
        INTERVAL_SEARCH em;
        em.initialize(&ds);
        stnodeTraverser.initialize(&em);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest<RDS>::testLCPIntervals(stnodeTraverser, &ds, plain_bwt);
        test_Intervals.swap(tmp);
    }

    //std::vector<stool::LCPInterval<uint64_t>> tmp = stool::lcp_on_rlbwt::Application<RDS>::testLCPIntervals(stnodeTraverser);

    stool::beller::equal_check_lcp_intervals(test_Intervals, stnc.lcp_intervals);
    std::cout << "LCP interval check OK!" << std::endl;
}
void computeLCPIntervals_beller(std::string inputFile)
{

    using LPOSDS = stool::EliasFanoVector;

    using FPOSDS = stool::lcp_on_rlbwt::LightFPosDataStructure;
    using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<INDEX, LPOSDS, FPOSDS>;
    //string text = "";
    std::cout << "Loading : " << inputFile << std::endl;
    std::vector<uint8_t> text = stool::load_text_from_file(inputFile, true);
    vector<INDEX> sa = stool::construct_suffix_array(text);
    sdsl::int_vector<> bwt;
    stool::FMIndex::constructBWT(text, sa, bwt);

    sdsl::bit_vector bv;
    bv.resize(bwt.size());
    bool b = true;
    //std::cout << bwt.size() << std::endl;
    //std::cout << "BV: ";
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

    std::vector<uint64_t> C;
    stool::FMIndex::constructC(bwt, C);

    wt_huff<> wt;
    construct_im(wt, bwt);

    uint64_t lastChar = bwt[bwt.size() - 1];

    stool::IntervalSearchDataStructure range;
    range.initialize(&wt, &C, lastChar);

    stool::lcp_on_rlbwt::ExplicitWeinerLinkSearch<uint32_t> wsearch;
    wsearch.initialize(&range, &bwt_bit_rank1, bwt.size());
    stool::lcp_on_rlbwt::SingleSTNodeTraverser<uint32_t, stool::lcp_on_rlbwt::ExplicitWeinerLinkSearch<uint32_t>> traverser;
    traverser.initialize(&wsearch);
    auto test_Intervals = LCPIntervalTest<RDS>::testLCPIntervals(traverser, nullptr, bwt);

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
    if(mode == "E"){
        computeLCPIntervals_beller(inputFile);
    }else{
    testLCPIntervals<INDEX>(inputFile, mode, thread_num);
    }
}
