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

        class Application
        {
        public:

            template <typename STNODES>
            static uint64_t outputMaximalSubstrings(std::ofstream &out, STNODES &stnodeSequencer, STTreeAnalysisResult &analysis)
            {
                using INDEX_SIZE = typename STNODES::index_type;

                uint64_t count = 0;

                std::vector<stool::LCPInterval<INDEX_SIZE>> buffer;

                auto it = stnodeSequencer.begin();
                while (it != stnodeSequencer.end())
                {
                    if(analysis.max_nodes_at_level < stnodeSequencer.child_count()){
                        analysis.max_nodes_at_level = stnodeSequencer.child_count();
                    }
                    std::vector<stool::LCPInterval<uint64_t>> r2;

                    for (auto node_it = it.begin(); node_it != it.end(); node_it++)
                    {
                        bool b = stnodeSequencer.check_maximal_repeat(node_it);
                        if (b)
                        {
                            INDEX_SIZE left = stnodeSequencer.get_left(node_it);
                            INDEX_SIZE right = stnodeSequencer.get_right(node_it);
                            stool::LCPInterval<INDEX_SIZE> intv(left, right, stnodeSequencer.get_current_lcp());
                            buffer.push_back(intv);
                            count++;
                            if (buffer.size() >= 8000)
                            {
                                out.write(reinterpret_cast<const char *>(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                                buffer.clear();
                            }
                        }
                    }
                    it++;
                }
                
                if (buffer.size() >= 1)
                {
                    out.write(reinterpret_cast<const char *>(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                    buffer.clear();
                }
                std::cout << "Enumerated" << std::endl;


                return count;
            }
            
        };

    } // namespace lcp_on_rlbwt
} // namespace stool