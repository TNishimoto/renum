#pragma once

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

#include "../hyper/application.hpp"
#include "../hyper/rlcp_interval_enumerator.hpp"
#include "../rlbwt/io.hpp"
#include "../rlbwt/fpos_data_structure.hpp"
#include "../rlbwt/bwt_decompress.hpp"
#include "../rlbwt/light_rlbwt.hpp"
#include "../main/common.hpp"



template <typename INDEX>
void testMaximalSubstrings(std::string inputFile, string mode, int thread_num)
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
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<INDEX, LPOSDS, FPOSDS>;
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
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<INDEX, LPOSDS, FPOSDS>;
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
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<INDEX, LPOSDS, FPOSDS>;
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
            using RDS = stool::lcp_on_rlbwt::RLBWTDataStructures<INDEX, LPOSDS, FPOSDS>;
            FPOSDS fposds = stool::lcp_on_rlbwt::LightFPosDataStructure(diff_char_vec, lpos_vec, wt);
            RDS ds = RDS(diff_char_vec, wt, lpos_vec, fposds);
            ds.stnc = &stnc;

            stool::lcp_on_rlbwt::ParallelSTNodeWTraverser<INDEX, RDS> stnodeTraverser;
            stnodeTraverser.initialize(thread_num, ds);
            std::cout << "TEST" << std::endl;

            std::vector<stool::LCPInterval<uint64_t>> tmp = stool::lcp_on_rlbwt::Application<RDS>::testLCPIntervals(stnodeTraverser);
            test_Intervals.swap(tmp);
        }
    }

    stool::beller::equal_check_lcp_intervals(test_Intervals, stnc.lcp_intervals);
    std::cout << "LCP interval check OK!" << std::endl;
}