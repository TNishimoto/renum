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

template <typename RLBWTDS>
class LCPIntervalTest
{
public:
    using CHAR = typename RLBWTDS::CHAR;
    //using CHARVEC = typename RLBWT_STR::char_vec_type;
    using INDEX_SIZE = typename RLBWTDS::INDEX;
    using UCHAR = typename std::make_unsigned<CHAR>::type;

    template <typename STNODES>
    static std::vector<stool::LCPInterval<uint64_t>> testLCPIntervals(STNODES &stnodeSequencer, RLBWTDS *ds)
    {

        std::vector<stool::LCPInterval<uint64_t>> r;

        while (!stnodeSequencer.isStop())
        {
            std::vector<stool::LCPInterval<uint64_t>> r2;

            stnodeSequencer.process();

            stool::LCPInterval<uint64_t> it;
            it.i = 0;
            it.j = 0;
            it.lcp = 0;

            stnodeSequencer.get_lcp_intervals(r2);
            /*
            if(ds->str_size() < 1000){
                for(auto &it : r2){
                    std::cout << it.to_string();
                }
                std::cout << std::endl;
            }
            */

            if (ds->stnc != nullptr)
            {
                ds->stnc->increment(r2);
            }
            for (auto &it : r2)
            {
                r.push_back(it);
            }
        }

        std::cout << "STOP" << std::endl;
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

    if (mode == "A")
    {
        std::cout << "General Test" << std::endl;
        stool::lcp_on_rlbwt::SuffixTreeNodes<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest<RDS>::testLCPIntervals(stnodeTraverser, &ds);
        test_Intervals.swap(tmp);
    }
    else if (mode == "B")
    {
        std::cout << "Standard Test" << std::endl;
        stool::lcp_on_rlbwt::STNodeTraverser<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest<RDS>::testLCPIntervals(stnodeTraverser, &ds);
        test_Intervals.swap(tmp);
    }
    else if (mode == "C")
    {
        std::cout << "Fast Test" << std::endl;
        stool::lcp_on_rlbwt::FastSTNodeTraverser<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest<RDS>::testLCPIntervals(stnodeTraverser, &ds);
        test_Intervals.swap(tmp);
    }

    else
    {
        std::cout << "Single Test" << std::endl;
        stool::lcp_on_rlbwt::SingleSTNodeTraverser<INDEX, RDS> stnodeTraverser;
        stnodeTraverser.initialize(&ds);
        std::vector<stool::LCPInterval<uint64_t>> tmp = LCPIntervalTest<RDS>::testLCPIntervals(stnodeTraverser, &ds);
        test_Intervals.swap(tmp);

    }

    //std::vector<stool::LCPInterval<uint64_t>> tmp = stool::lcp_on_rlbwt::Application<RDS>::testLCPIntervals(stnodeTraverser);

    stool::beller::equal_check_lcp_intervals(test_Intervals, stnc.lcp_intervals);
    std::cout << "LCP interval check OK!" << std::endl;
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

    testLCPIntervals<INDEX>(inputFile, mode, thread_num);
}
