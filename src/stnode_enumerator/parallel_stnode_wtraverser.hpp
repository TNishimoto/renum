#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>
#include <thread>
#include "stool/src/byte.hpp"
#include <cmath>
#include "../rlbwt/range_distinct/succinct_range_distinct.hpp"

#include "standard/stnode_wtraverser.hpp"
#include "succinct/succinct_sorted_stchildren_builder.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {

        template <typename INDEX_SIZE, typename RLBWTDS>
        class ParallelSTNodeWTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using STNODE_WTRAVERSER = STNodeSubWTraverser<INDEX_SIZE, RLBWTDS>;
            using SUCCINCT_STNODE_WTRAVERSER = SuccinctSortedSTChildren<INDEX_SIZE, RLBWTDS>;

            STNodeWTraverser<INDEX_SIZE, RLBWTDS> standard_st_traverser;

            std::vector<SUCCINCT_STNODE_WTRAVERSER *> succinct_sub_trees;

            //std::stack<STNODE_WTRAVERSER> empty_trees;

            //std::vector<STNODE_WTRAVERSER> sub_tmp_trees;
            std::vector<ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> ems;
            std::vector<LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE>> lightRDs;
            std::vector<SuccinctRangeDistinctDataStructure<INDEX_SIZE>> heavyRDs;

            uint64_t print_interval = 100;
            uint64_t print_interval_counter = 0;

        public:
            //uint64_t current_lcp = 0;
            uint64_t total_counter = 0;
            //uint64_t strSize = 0;
            //uint64_t node_count = 0;


            //uint64_t child_count = 0;
            uint64_t peak_child_count = 0;
            uint64_t alpha = 2;
            uint64_t _expected_peak_memory_bits = 0;
            uint64_t _switch_threshold = 0;
            //uint64_t thread_count = 1;
            uint64_t debug_peak_memory = 0;

            bool succinct_mode = false;

            RLBWTDS *_RLBWTDS;
            //bool is_parallel = true;
            /*
            uint64_t count_maximal_repeats()
            {
                return 0;
                assert(false);

                throw -1;
            }
            */
            /*
            STNODE_WTRAVERSER *get_sub_tree()
            {
                return &this->sub_tree;
            }
            */
            uint64_t node_count() const
            {
                if (this->succinct_mode)
                {
                    uint64_t k = 0;
                    for (auto &it : this->succinct_sub_trees)
                    {
                        k += it->node_count();
                    }
                    return k;
                }
                else
                {
                    return this->standard_st_traverser.node_count();
                }

                //return this->sub_tree.node_count();
            }
            uint64_t child_count() const {
                return this->standard_st_traverser.child_count();
            }

            uint64_t expected_peak_memory_bits()
            {
                return this->_expected_peak_memory_bits;
            }

            uint64_t switch_threshold()
            {
                return this->_switch_threshold;
            }

            void initialize(uint64_t size, RLBWTDS &_RLBWTDS)
            {
                this->_RLBWTDS = &_RLBWTDS;
                //this->strSize = _RLBWTDS.str_size();
                this->standard_st_traverser.initialize(size, _RLBWTDS);

                /*
                double ratio = (double)this->_RLBWTDS->str_size() / (double)this->_RLBWTDS->rle_size();
                double d = std::log2(ratio);
                this->_expected_peak_memory_bits = this->_RLBWTDS->rle_size() * d;
                this->_switch_threshold = this->alpha * (this->expected_peak_memory_bits() / (sizeof(uint64_t) * 8));
                */
            }
            uint64_t write_maximal_repeats(std::ofstream &out)
            {
                return this->standard_st_traverser.write_maximal_repeats(out);
            }
            uint64_t current_lcp(){
                return this->standard_st_traverser.get_current_lcp();
            }
            uint64_t get_stnode2(uint64_t L, stool::LCPInterval<uint64_t> &output)
            {
                if (!this->succinct_mode)
                {
                    return this->standard_st_traverser.get_stnode2(L, output);
                }
                else
                {
                    uint64_t p = 0;
                    for (uint64_t i = 0; i < this->succinct_sub_trees.size(); i++)
                    {
                        uint64_t csize = this->succinct_sub_trees[i]->children_count();

                        if (p <= L && L < p + csize)
                        {
                            return p + this->succinct_sub_trees[i]->get_stnode2(L - p, output, this->current_lcp() - 1);
                        }
                        else
                        {
                            p += csize;
                        }
                    }
                    return UINT64_MAX;
                }
            }

            uint64_t kk = 0;
#if DEBUG
            uint64_t prev_child_count = 0;
#endif

            void lightEnumerate()
            {
#if DEBUG
                std::cout << "LIGHT LCP = " << this->current_lcp() << std::endl;
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

            void process()
            {
                if (this->total_counter > 0)
                {
                    uint64_t ccc = this->_RLBWTDS->str_size() / this->print_interval;
                    uint64_t pp_num = ccc * this->print_interval_counter;

                    if (this->debug_peak_memory < this->get_using_memory())
                    {
                        this->debug_peak_memory = this->get_using_memory();
                    }

                    if (this->total_counter >= pp_num)
                    {
                        std::cout << "[" << (this->print_interval_counter) << "/" << this->print_interval << "] ";
                        std::cout << "LCP = " << this->current_lcp();
                        std::cout << ", Peak = " << this->peak_child_count;
                        std::cout << ", Current = " << this->child_count();
                        std::cout << ", current_total = " << (this->_RLBWTDS->str_size() - this->total_counter);
                        std::cout << ", Peak = " << this->debug_peak_memory / 1000 << "[KB]" << std::endl;
                        //std::cout << "Peak = " << debug_peak_counter << "[KB]" << std::endl;
                        /*
                        if(this->print_interval_counter == 2){
                            throw -1;
                        }
                        */

                        this->print_interval_counter++;
                    }
                }
#if DEBUG
                if (prev_child_count != this->child_count())
                {
                    std::cout << "LCP = " << current_lcp() << "/" << this->child_count() << "/" << std::endl;
                    prev_child_count = this->child_count();
                }
                if (this->_RLBWTDS->str_size() < 100)
                {
                    std::cout << "Start Enumerate" << std::endl;
                    this->print();
                }
#endif

                if (this->current_lcp() == 0)
                {
                    this->standard_st_traverser.process();
                }
                else
                {

                    bool isHeavy = true;

                    if (isHeavy)
                    {
                        this->standard_st_traverser.process();
                    }
                    else
                    {
                        this->lightEnumerate();
                    }
                }
                this->update_info();
#if DEBUG

                if (this->_RLBWTDS->str_size() < 100)
                {
                    std::cout << "Enumerate END" << std::endl;
                    this->print();
                }

#endif
            }
            void update_info()
            {
                this->total_counter += this->child_count() - this->node_count();
                if(this->peak_child_count < this->child_count()){
                    this->peak_child_count = this->child_count();
                }


            }

            bool isStop()
            {
                return this->standard_st_traverser.isStop();
                //return total_counter == strSize - 1;
            }

            void print()
            {
                this->standard_st_traverser.print();
            }

        private:
            uint64_t get_using_memory()
            {
                return this->standard_st_traverser.get_using_memory();
            }
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
        };

    } // namespace lcp_on_rlbwt
} // namespace stool