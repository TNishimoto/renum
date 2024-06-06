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
#include "../explicit_weiner_link_computer_on_rlbwt.hpp"
#include "stool/include/elias_fano_vector.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace stnode_on_rlbwt
    {
        /*
            void build_succinct()
            {
                std::cout << "BUILD SUccinct" << std::endl;
                std::vector<std::pair<uint64_t, uint64_t>> children;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    for (uint64_t j = 0; j < this->sub_trees[i]->children_count(); j++)
                    {
                        children.push_back(std::pair<uint64_t, uint64_t>(i, j));
                    }
                }
                sort(children.begin(), children.end(), [&](const std::pair<uint64_t, uint64_t> &lhs, const std::pair<uint64_t, uint64_t> &rhs) {
                    auto &left = sub_trees[lhs.first]->childVec[lhs.second];
                    auto &right = sub_trees[rhs.first]->childVec[rhs.second];
                    uint64_t begin_pos1 = _RLBWTDS->_fposDS[left.beginIndex] + left.beginDiff;
                    uint64_t begin_pos2 = _RLBWTDS->_fposDS[right.beginIndex] + right.beginDiff;
                    return begin_pos1 < begin_pos2;
                });
                auto wBuilder = new SuccinctSortedSTChildrenBuilder<INDEX_SIZE, RLBWTDS>();
                std::vector<SuccinctSortedSTChildrenBuilder<INDEX_SIZE, RLBWTDS> *> wBuilders;

                wBuilders.push_back(wBuilder);
                //SuccinctSortedSTChildrenBuilder<INDEX_SIZE, RLBWTDS> wBuilder;
                wBuilder->initialize(0, this->node_count(), this->_RLBWTDS, nullptr);
                wBuilder->set(this->_RLBWTDS->str_size(), children.size());

                for (auto &it : children)
                {
                    auto &item = sub_trees[it.first]->childVec[it.second];
                    uint8_t c = this->_RLBWTDS->bwt[item.beginIndex];
                    uint64_t begin_pos = this->_RLBWTDS->_fposDS[item.beginIndex] + item.beginDiff;
                    uint64_t end_pos = this->_RLBWTDS->_fposDS[item.endIndex] + item.endDiff;
                    bool isLeft = sub_trees[it.first]->w_builder[it.second];

                    if (this->current_lcp == 1)
                    {
                        isLeft = begin_pos == 0;
                    }

                    LightweightInterval newIntv(begin_pos, end_pos, isLeft);
                    wBuilder->push(newIntv, c);
                }

                SuccinctSortedSTChildren<INDEX_SIZE, RLBWTDS> *succ = new SuccinctSortedSTChildren<INDEX_SIZE, RLBWTDS>();

                SuccinctSortedSTChildrenBuilder<INDEX_SIZE, RLBWTDS>::merge(*succ, wBuilders);
                //wBuilder.buildSuccinctSortedSTChildren(*succ);
                this->succinct_sub_trees.push_back(succ);
                std::cout << "Memory: " << (succ->get_using_memory() / 1000) << "[KB]" << std::endl;
                delete wBuilder;
            }
            */
           /*
           void lightEnumerate()
            {
#if DEBUG
                if (this->_RLBWTDS->str_size() < 100)
                {
                    this->succinct_sub_trees[0]->print();
                }
#endif
                std::vector<SuccinctSortedSTChildrenBuilder<INDEX_SIZE, RLBWTDS> *> wBuilders;

                uint64_t psize = (this->node_count() / 2) + 1;
                //uint64_t psize = this->node_count() + 1;

                uint64_t px = 0;
                while (px < this->node_count())
                {
                    uint64_t xsize = px + psize <= this->node_count() ? psize : this->node_count() - px;
                    auto wBuilder = new SuccinctSortedSTChildrenBuilder<INDEX_SIZE, RLBWTDS>();
                    wBuilder->initialize(px, px + xsize - 1, this->_RLBWTDS, this->succinct_sub_trees[0]);
                    wBuilders.push_back(wBuilder);
                    px += xsize;
                }

                uint64_t next_child_count = 0;
                for (auto &it : wBuilders)
                {

                    next_child_count += it->countNextSTNodes(this->ems[0]);

                    it->set();

                    it->computeNextSTNodes(this->ems[0]);
                }

                SuccinctSortedSTChildren<INDEX_SIZE, RLBWTDS> *succ = new SuccinctSortedSTChildren<INDEX_SIZE, RLBWTDS>();
                SuccinctSortedSTChildrenBuilder<INDEX_SIZE, RLBWTDS>::merge(*succ, wBuilders);
                for (uint64_t i = 0; i < wBuilders.size(); i++)
                {
                    delete wBuilders[i];
                }
                delete this->succinct_sub_trees[0];

                this->succinct_sub_trees[0] = succ;

                assert(false);
                //this->recompute_node_counter();

                assert(total_counter <= this->_RLBWTDS->str_size());

                assert(this->child_count() > 0);
                assert(this->node_count() > 0);
            }
            */
    }
}