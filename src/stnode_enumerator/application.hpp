#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>
//#include "next_rinterval_storage_constructor.hpp"
#include "../rlbwt/rlbwt_data_structures.hpp"

#include "suffix_tree_nodes.hpp"
#include "weiner_link_emulator.hpp"
#include <thread>
#include "../debug/stnode_checker.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        class STTreeAnalysisResult
        {
        public:
            uint64_t max_nodes_at_level;
        };

        template <typename RLBWTDS>
        class Application
        {
        public:
            using CHAR = typename RLBWTDS::CHAR;
            //using CHARVEC = typename RLBWT_STR::char_vec_type;
            using INDEX_SIZE = typename RLBWTDS::INDEX;
            using UCHAR = typename std::make_unsigned<CHAR>::type;
            using RINTERVAL = RInterval<INDEX_SIZE>;

            static uint64_t outputMaximalSubstrings(std::ofstream &out, SuffixTreeNodes<INDEX_SIZE, RLBWTDS> &stnodeSequencer, STTreeAnalysisResult &analysis)
            {

                uint64_t count = 0;

                while (!stnodeSequencer.is_finished())
                {

                    stnodeSequencer.succ();
                    count += stnodeSequencer.write_maximal_repeats(out);
                }
                std::cout << "Enumerated" << std::endl;
                analysis.max_nodes_at_level = stnodeSequencer.peak_child_count;

                return count;
            }
            static std::vector<stool::LCPInterval<uint64_t>> testLCPIntervals(SuffixTreeNodes<INDEX_SIZE, RLBWTDS> &stnodeSequencer)
            {

                std::vector<stool::LCPInterval<uint64_t>> r;

                while (!stnodeSequencer.is_finished())
                {
                    std::vector<stool::LCPInterval<uint64_t>> r2;

                    stnodeSequencer.succ();

                    stool::LCPInterval<uint64_t> it;
                    it.i = 0;
                    it.j = 0;
                    it.lcp = 0;

                    stnodeSequencer.get_lcp_intervals(r2);

                    if (stnodeSequencer._RLBWTDS->stnc != nullptr)
                    {
                        stnodeSequencer._RLBWTDS->stnc->increment(r2);
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

    } // namespace lcp_on_rlbwt
} // namespace stool