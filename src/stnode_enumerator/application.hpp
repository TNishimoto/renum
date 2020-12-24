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

#include "stnode_enumerator.hpp"
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

            static uint64_t outputMaximalSubstrings(std::ofstream &out, STNodeEnumerator<INDEX_SIZE, RLBWTDS> &stnodeSequencer, STTreeAnalysisResult &analysis)
            {

                uint64_t count = 0;

                while (!stnodeSequencer.isStop())
                {

                    stnodeSequencer.process();
                    count += stnodeSequencer.write_maximal_repeats(out);

                }
                std::cout << "Enumerated" << std::endl;
                analysis.max_nodes_at_level = stnodeSequencer.peak_child_count;

                return count;
            }
            static std::vector<stool::LCPInterval<uint64_t>> testLCPIntervals(STNodeEnumerator<INDEX_SIZE, RLBWTDS> &stnodeSequencer)
            {

                std::vector<stool::LCPInterval<uint64_t>> r;

                while (!stnodeSequencer.isStop())
                {
                    stnodeSequencer.process();

                    stool::LCPInterval<uint64_t> it;
                    it.i = 0;
                    it.j = 0;
                    it.lcp = 0;

                    uint64_t L = 0;

                    for (uint64_t i = 0; i < stnodeSequencer.node_count(); i++)
                    {

                        L = stnodeSequencer.get_stnode(L, it);
                        r.push_back(it);
                    }
                    if (stnodeSequencer._RLBWTDS->stnc != nullptr)
                    {
                        stnodeSequencer._RLBWTDS->stnc->increment(stnodeSequencer.node_count());
                    }
                }
                std::cout << "STOP" << std::endl;
                return r;
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool