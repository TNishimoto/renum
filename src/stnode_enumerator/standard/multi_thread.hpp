#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <deque>
#include <vector>
#include <type_traits>
#include "../weiner_link_emulator.hpp"
#include "stnode_sub_traverser.hpp"
#include "../single/single_stnode_traverser.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {

        std::mutex mtx;

        template <typename INDEX_SIZE, typename RLBWTDS>
        void distribute(std::vector<STNodeSubTraverser<INDEX_SIZE, RLBWTDS> *> &trees, STNodeVector<INDEX_SIZE> &tmp, std::queue<STNodeSubTraverser<INDEX_SIZE, RLBWTDS> *> &uqueue, uint64_t limit)
        {
            bool store_edge_chars = trees[0]->has_edge_characters();
            while (tmp.size() > 0)
            {
                uint64_t w = tmp.get_last_width();
                if (uqueue.size() > 0)
                {
                    auto top = uqueue.front();
                    if (top->childvec_size() + w <= top->capacity() )
                    {
                        top->move_push(tmp);
                    }
                    else
                    {
                        uqueue.pop();
                    }
                }
                else
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    auto st = new STNodeSubTraverser<INDEX_SIZE, RLBWTDS>(limit, store_edge_chars);
                    trees.push_back(st);
                    uqueue.push(st);
                }
            }
        }

        template <typename INDEX_SIZE, typename RLBWTDS>
        void parallel_succ_stnodes(std::vector<STNodeSubTraverser<INDEX_SIZE, RLBWTDS> *> &trees, uint64_t fst_position,
                                      std::stack<uint64_t> &position_stack, ExplicitWeinerLinkEmulator<RLBWTDS> &em, uint64_t limit)
        {
            STNodeVector<INDEX_SIZE> tmp;
            //tmp.store_edge_chars = trees[0]->has_edge_characters();
            std::queue<STNodeSubTraverser<INDEX_SIZE, RLBWTDS> *> uqueue;
            uint64_t pos = fst_position;
            while (pos != UINT64_MAX)
            {

                trees[pos]->computeNextSTNodes(em, tmp);
                trees[pos]->clear();
                uqueue.push(trees[pos]);
                distribute(trees, tmp, uqueue, limit);

                //assert(trees[pos]->children_count() <= limit);
                pos = UINT64_MAX;
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    if (position_stack.size() > 0)
                    {
                        pos = position_stack.top();
                        position_stack.pop();
                    }
                }
            }
        }

        template <typename INDEX_SIZE, typename RLBWTDS>
        void parallel_count_stnodes(std::vector<STNodeSubTraverser<INDEX_SIZE, RLBWTDS> *> &trees, uint64_t fst_position,
                                    std::stack<uint64_t> &position_stack, ExplicitWeinerLinkEmulator<RLBWTDS> &em, uint64_t &output_children_count)
        {
            output_children_count = 0;
            uint64_t pos = fst_position;
            while (pos != UINT64_MAX)
            {
                auto pair = trees[pos]->countNextLCPIntervalSet(em);
                output_children_count += pair.second;
                pos = UINT64_MAX;
                {
                    std::lock_guard<std::mutex> lock(mtx);
                    if (position_stack.size() > 0)
                    {
                        pos = position_stack.top();
                        position_stack.pop();
                    }
                }
            }
        }
    } // namespace lcp_on_rlbwt
} // namespace stool