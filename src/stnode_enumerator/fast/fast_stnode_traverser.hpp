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
#include "fast_stnode_sub_traverser.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        std::mutex mtx3;
        template <typename INDEX_SIZE, typename RLBWTDS>
        void parallel_process_fast_stnodes(std::vector<FastSTNodeSubTraverser<INDEX_SIZE, RLBWTDS> *> &trees, uint64_t fst_position,
                                           std::stack<uint64_t> &position_stack, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em, uint64_t limit)
        {

            uint64_t pos = fst_position;
            while (pos != UINT64_MAX)
            {
                trees[pos]->computeNextLCPIntervalSet(em, limit);
                pos = UINT64_MAX;
                {
                    std::lock_guard<std::mutex> lock(mtx3);
                    if (position_stack.size() > 0)
                    {
                        pos = position_stack.top();
                        position_stack.pop();
                    }
                }
            }
        }
        template <typename INDEX_SIZE, typename RLBWTDS>
        class FastSTNodeTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using STNODE_SUB_TRAVERSER = FastSTNodeSubTraverser<INDEX_SIZE, RLBWTDS>;

            std::vector<STNODE_SUB_TRAVERSER *> sub_trees;

            std::deque<std::deque<INDEX_SIZE>> children_intervals;
            std::deque<std::deque<bool>> first_child_flag_vec;
            std::deque<std::deque<bool>> maximal_repeat_check_vec;

            std::vector<ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> ems;
            //uint64_t minimum_child_count = 1000;
            //uint64_t sub_tree_limit_size = 2000;
            uint64_t mmax = 0;
            uint64_t peak_mmax = 10000;

            uint64_t sub_tree_max_count = 100;

            uint64_t current_lcp = 0;
            uint64_t _child_count = 0;
            uint64_t _node_count = 0;

            uint64_t thread_count = 1;
            std::stack<uint64_t> position_stack;
            uint64_t kk = 0;
#if DEBUG
            uint64_t prev_child_count = 0;
#endif

        public:
            RLBWTDS *_RLBWTDS;

            void set_peak(uint64_t _peak){
                this->peak_mmax = _peak;
            }

            uint64_t get_current_lcp() const
            {
                return current_lcp;
            }
            uint64_t child_count() const
            {
                return this->_child_count;
            }

            uint64_t node_count() const
            {
                return this->_node_count;
            }

            void initialize(uint64_t size, RLBWTDS &__RLBWTDS)
            {
                this->_RLBWTDS = &__RLBWTDS;

                this->thread_count = size;

                auto st = new STNODE_SUB_TRAVERSER(this->_RLBWTDS);
                sub_trees.push_back(st);
                //sub_trees.resize(256);

                for (uint64_t i = 0; i < this->thread_count; i++)
                {
                    ems.push_back(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>());
                    ems[ems.size() - 1].initialize(&__RLBWTDS);
                }
            }
            uint64_t write_maximal_repeats(std::ofstream &out)
            {
                uint64_t count = 0;
                std::vector<stool::LCPInterval<INDEX_SIZE>> buffer;

                uint64_t size = this->node_count();
                uint64_t L = 0;
                RINTERVAL intv;
                for (uint64_t i = 0; i < size; i++)
                {
                    L = this->get_stnode(L, intv);
                    if (this->maximal_repeat_check_vec[0][i])
                    {
                        uint64_t _left = this->_RLBWTDS->get_lpos(intv.beginIndex) + intv.beginDiff;
                        uint64_t _right = this->_RLBWTDS->get_lpos(intv.endIndex) + intv.endDiff;

                        stool::LCPInterval<INDEX_SIZE> newLCPIntv(_left, _right, this->current_lcp - 1);
                        buffer.push_back(newLCPIntv);
                        count++;
                        if (buffer.size() >= 8000)
                        {
                            out.write(reinterpret_cast<const char *>(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                            buffer.clear();
                        }
                    }
                }

                if (buffer.size() >= 1)
                {
                    out.write(reinterpret_cast<const char *>(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                    buffer.clear();
                }
                return count;
            }
            uint64_t get_stnode(uint64_t L, RINTERVAL &output)
            {
                assert(this->first_child_flag_vec.size() > 0);
                uint64_t R = L + 1;
                while (R < this->first_child_flag_vec[0].size() && !this->first_child_flag_vec[0][R])
                {
                    R++;
                }

                R--;

                assert(this->children_intervals.size() > 0);
                assert(this->children_intervals[0].size() > (R * 4) + 3);

                output.beginIndex = this->children_intervals[0][L * 4];
                output.beginDiff = this->children_intervals[0][(L * 4) + 1];
                output.endIndex = this->children_intervals[0][(R * 4) + 2];
                output.endDiff = this->children_intervals[0][(R * 4) + 3];

                assert(output.beginIndex < this->_RLBWTDS->rle_size());

                return R + 1;
            }

            void process()
            {
                bool isSingleProcess = false;

                if (this->first_child_flag_vec.size() > 0)
                {
                    mmax -= this->first_child_flag_vec[0].size();
                    this->children_intervals.pop_front();
                    this->first_child_flag_vec.pop_front();
                    this->maximal_repeat_check_vec.pop_front();
                }

                if (this->first_child_flag_vec.size() == 0)
                {
                    for (auto &it : this->sub_trees)
                    {
                        it->set_parent_current_lcp(this->current_lcp);
                    }

                    if (current_lcp == 0)
                    {
                        this->sub_trees[0]->first_compute(ems[0]);
                    }

                    isSingleProcess = this->thread_count == 1;
                    if (isSingleProcess)
                    {

                        this->single_process();
                    }
                    else
                    {

                        this->parallel_process();
                    }
                    this->rep();
                    this->remove_empty_trees();
                    if ((double)this->sub_trees.size() * 2 < (double)this->sub_trees.capacity())
                    {
                        this->sub_trees.shrink_to_fit();
                    }
                    this->split();
                }
                this->recompute_node_counter();
            }
            void rep()
            {
                uint64_t k = UINT64_MAX;
                assert(this->children_intervals.size() == 0);

                if (this->sub_trees.size() == 0)
                    return;

                if (this->current_lcp == 0)
                {
                    this->children_intervals.push_back(std::deque<INDEX_SIZE>());
                    this->first_child_flag_vec.push_back(std::deque<bool>());
                    this->maximal_repeat_check_vec.push_back(std::deque<bool>());

                    for (auto &it : this->sub_trees)
                    {
                        it->pop_level(this->children_intervals[0], this->first_child_flag_vec[0], this->maximal_repeat_check_vec[0]);
                    }
                }
                else
                {
                    for (auto &it : this->sub_trees)
                    {
                        uint64_t x = (it->get_current_lcp() - this->current_lcp) + it->finished_level_count();
                        if (x < k)
                        {
                            k = x;
                        }
                    }
                    assert(k != UINT64_MAX);
                    for (uint64_t i = 0; i < k; i++)
                    {
                        this->children_intervals.push_back(std::deque<INDEX_SIZE>());
                        this->first_child_flag_vec.push_back(std::deque<bool>());
                        this->maximal_repeat_check_vec.push_back(std::deque<bool>());
                    }
                    for (auto &it : this->sub_trees)
                    {
                        for (uint64_t i = 0; i < k; i++)
                        {
                            uint64_t klcp = this->current_lcp + i;
                            if (klcp == it->get_current_lcp())
                            {
                                it->pop_level(this->children_intervals[i], this->first_child_flag_vec[i], this->maximal_repeat_check_vec[i]);
                            }
                        }
                    }
                    for (uint64_t i = 0; i < k; i++)
                    {
                        mmax += this->first_child_flag_vec[i].size();
                    }
                }
                assert(this->first_child_flag_vec.size() > 0);
            }
            bool isStop()
            {
                return this->current_lcp > 0 && this->child_count() == 0;
                //return total_counter == strSize - 1;
            }

            void print()
            {
                std::cout << "PRINT PTREE" << std::endl;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    if (this->sub_trees[i]->get_current_lcp() == (this->current_lcp - 1))
                    {
                        this->sub_trees[i]->print_info();
                    }
                }

                std::cout << "[END]" << std::endl;
            }

            uint64_t get_using_memory()
            {
                uint64_t k = 0;

                for (auto &it : this->sub_trees)
                {
                    k += it->get_using_memory();
                }

                return k;
            }

            void import(std::deque<INDEX_SIZE> &item1, std::deque<bool> &item2, std::deque<bool> &item3, uint64_t lcp)
            {
                uint64_t k = 1000;
                while (this->sub_trees.size() > 0)
                {
                    delete this->sub_trees[this->sub_trees.size() - 1];
                    this->sub_trees.pop_back();
                }
                while (item3.size() > 0)
                {
                    auto st = new STNODE_SUB_TRAVERSER(this->_RLBWTDS);
                    sub_trees.push_back(st);
                    uint64_t num = item3.size() < k ? item3.size() : k;
                    st->import(item1, item2, item3, lcp, num);
                    
                }
                assert(item1.size() == 0);
                assert(item2.size() == 0);
                assert(item3.size() == 0);

                this->current_lcp = lcp;
            }

        private:
            void recompute_node_counter()
            {
                if (this->first_child_flag_vec.size() == 0)
                {
                    _node_count = 0;
                    _child_count = 0;
                }
                else
                {
                    _node_count = this->maximal_repeat_check_vec[0].size();
                    _child_count = this->first_child_flag_vec[0].size();
                }

                current_lcp++;
            }
            void parallel_process()
            {
#if DEBUG
                std::cout << "PARALLEL PROCESS" << std::endl;
#endif
                std::vector<uint64_t> fst_pos_vec;
                fst_pos_vec.resize(this->thread_count, UINT64_MAX);
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    if (i < this->thread_count)
                    {
                        fst_pos_vec[i] = i;
                    }
                    else
                    {
                        position_stack.push(i);
                    }
                }

                //auto start = std::chrono::system_clock::now();
                uint64_t limit = (peak_mmax - this->mmax) / (this->sub_trees.size() + 1);

                std::vector<thread> threads;
                for (uint64_t i = 0; i < this->thread_count; i++)
                {
                    threads.push_back(thread(parallel_process_fast_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees), fst_pos_vec[i], ref(position_stack), ref(ems[i]), limit));
                }

                for (thread &t : threads)
                    t.join();
                    //auto end = std::chrono::system_clock::now();
#if DEBUG
                std::cout << "PARALLEL PROCESS[END]" << std::endl;
#endif

                assert(position_stack.size() == 0);
            }
            void single_process()
            {
                if (this->sub_trees.size() == 0)
                    return;
                uint64_t fst_pos = 0;
                for (uint64_t i = 1; i < this->sub_trees.size(); i++)
                {
                    position_stack.push(i);
                }


                uint64_t limit = peak_mmax - this->mmax;

                parallel_process_fast_stnodes<INDEX_SIZE, RLBWTDS>(ref(sub_trees), fst_pos, ref(position_stack), ref(ems[0]), limit);
                assert(position_stack.size() == 0);
            }
            void split()
            {
                //return;
                if (this->sub_trees.size() == 0)
                    return;
                sort(this->sub_trees.begin(), this->sub_trees.end(), [&](const STNODE_SUB_TRAVERSER *lhs, const STNODE_SUB_TRAVERSER *rhs) {
                    int x1 = lhs->last_node_count() == 0 ? 1 : 0;
                    int x2 = rhs->last_node_count() == 0 ? 1 : 0;
                    if (x1 < x2)
                    {
                        return true;
                    }
                    else if (x1 > x2)
                    {
                        return false;
                    }
                    else
                    {
                        return lhs->get_last_lcp() < rhs->get_last_lcp();
                    }
                });

                //assert(this->sub_trees[0]->last_node_count() > 0);

                while (this->sub_trees[0]->last_node_count() > 1 && this->sub_trees.size() < this->sub_tree_max_count)
                {
                    auto st = new STNODE_SUB_TRAVERSER(this->_RLBWTDS);
                    this->sub_trees[0]->split(*st);
                    this->sub_trees.push_back(st);
                }
                //assert(this->sub_trees[0]->last_node_count() > 0);
            }

            void remove_empty_trees()
            {

                uint64_t nonEmptyCount = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    if (!this->sub_trees[i]->is_empty())
                    {
                        if (i != nonEmptyCount)
                        {
                            auto tmp = this->sub_trees[nonEmptyCount];
                            this->sub_trees[nonEmptyCount] = this->sub_trees[i];
                            this->sub_trees[i] = tmp;
                        }
                        nonEmptyCount++;
                    }
                }
                for (uint64_t i = nonEmptyCount; i < this->sub_trees.size(); i++)
                {
                    this->sub_trees[i]->clear();
                    delete this->sub_trees[i];
                    this->sub_trees[i] = nullptr;
                }
                this->sub_trees.resize(nonEmptyCount);
            }
        }; // namespace lcp_on_rlbwt

    } // namespace lcp_on_rlbwt
} // namespace stool