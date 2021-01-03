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
        void parallel_succ_fast_stnodes(std::vector<FastSTNodeSubTraverser<INDEX_SIZE, RLBWTDS> *> &trees, uint64_t fst_position,
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

            using ITERATOR = STNodeIterator<FastSTNodeTraverser>;
            using DEPTH_ITERATOR = STDepthIterator<FastSTNodeTraverser>;
            std::vector<STNODE_SUB_TRAVERSER *> sub_trees;

            std::deque<INDEX_SIZE> children_intervals;
            std::deque<bool> first_child_flag_vec;
            std::deque<bool> maximal_repeat_check_vec;
            std::deque<uint64_t> child_width_vec;
            std::deque<uint64_t> st_width_vec;

            std::vector<ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> ems;
            //uint64_t minimum_child_count = 1000;
            //uint64_t sub_tree_limit_size = 2000;
            uint64_t mmax = 0;
            uint64_t peak_mmax = 10000;

            uint64_t sub_tree_max_count = 100;

            int64_t current_lcp = -1;
            //uint64_t _child_count = 0;
            //uint64_t _node_count = 0;

            uint64_t thread_count = 1;
            std::stack<uint64_t> position_stack;
            uint64_t kk = 0;
#if DEBUG
            uint64_t prev_child_count = 0;
#endif

        public:
            using index_type = INDEX_SIZE;
            RLBWTDS *_RLBWTDS;

            void set_peak(uint64_t _peak)
            {
                this->peak_mmax = _peak;
            }
            void clear()
            {

                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    this->sub_trees[i]->clear();
                    delete this->sub_trees[i];
                    this->sub_trees[i] = nullptr;
                }
                this->sub_trees.clear();
                this->first_child_flag_vec.clear();
                this->maximal_repeat_check_vec.clear();
                this->child_width_vec.clear();
                this->st_width_vec.clear();
                this->current_lcp = -1;
            }

            uint64_t get_current_lcp() const
            {
                return current_lcp;
            }
            uint64_t child_count() const
            {
                return this->child_width_vec.size() == 0 ? 0 : this->child_width_vec[0];
            }

            uint64_t node_count() const
            {
                return this->st_width_vec.size() == 0 ? 0 : this->st_width_vec[0];
            }

            void initialize(uint64_t size, RLBWTDS &__RLBWTDS)
            {
                this->_RLBWTDS = &__RLBWTDS;

                this->thread_count = size;

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
                    if (this->maximal_repeat_check_vec[i])
                    {
                        uint64_t _left = this->_RLBWTDS->get_lpos(intv.beginIndex) + intv.beginDiff;
                        uint64_t _right = this->_RLBWTDS->get_lpos(intv.endIndex) + intv.endDiff;

                        stool::LCPInterval<INDEX_SIZE> newLCPIntv(_left, _right, this->current_lcp);
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

            void get_lcp_intervals(std::vector<stool::LCPInterval<uint64_t>> &output)
            {

                std::vector<stool::LCPInterval<INDEX_SIZE>> buffer;

                uint64_t size = this->node_count();
                uint64_t L = 0;
                RINTERVAL intv;
                for (uint64_t i = 0; i < size; i++)
                {
                    L = this->get_stnode(L, intv);
                    uint64_t _left = this->_RLBWTDS->get_lpos(intv.beginIndex) + intv.beginDiff;
                    uint64_t _right = this->_RLBWTDS->get_lpos(intv.endIndex) + intv.endDiff;
                    stool::LCPInterval<INDEX_SIZE> newLCPIntv(_left, _right, this->current_lcp);
                    output.push_back(newLCPIntv);
                }
            }
            uint64_t get_stnode(uint64_t L, RINTERVAL &output) const
            {
                assert(this->first_child_flag_vec.size() > 0);
                uint64_t R = L + 1;
                while (R < this->first_child_flag_vec.size() && !this->first_child_flag_vec[R])
                {
                    R++;
                }

                R--;

                assert(this->children_intervals.size() > 0);
                assert(this->children_intervals.size() > (R * 4) + 3);

                output.beginIndex = this->children_intervals[L * 4];
                output.beginDiff = this->children_intervals[(L * 4) + 1];
                output.endIndex = this->children_intervals[(R * 4) + 2];
                output.endDiff = this->children_intervals[(R * 4) + 3];

                assert(output.beginIndex < this->_RLBWTDS->rle_size());

                return R + 1;
            }
            void pop()
            {
                uint64_t k1 = this->node_count();
                uint64_t k2 = this->child_count();

                for (uint64_t i = 0; i < k1; i++)
                {
                    this->maximal_repeat_check_vec.pop_back();
                }
                for (uint64_t i = 0; i < k2; i++)
                {
                    this->children_intervals.pop_front();
                    this->children_intervals.pop_front();
                    this->children_intervals.pop_front();
                    this->children_intervals.pop_front();
                    this->first_child_flag_vec.pop_front();
                }

                mmax -= k2;
                this->child_width_vec.pop_front();
                this->st_width_vec.pop_front();
                //this->children_intervals.pop_front();
                //this->first_child_flag_vec.pop_front();
                //this->maximal_repeat_check_vec.pop_front();
            }

            bool succ()
            {
                if (this->is_finished())
                    return false;
                bool isSingleProcess = false;

                if (this->first_child_flag_vec.size() > 0)
                {
                    this->pop();
                }

                if (this->first_child_flag_vec.size() == 0)
                {
                    for (auto &it : this->sub_trees)
                    {
                        it->set_parent_current_lcp(this->current_lcp + 1);
                    }

                    if (current_lcp == -1)
                    {
                        auto st = new STNODE_SUB_TRAVERSER(this->_RLBWTDS);
                        sub_trees.push_back(st);

                        this->sub_trees[0]->first_compute(ems[0]);
                    }

                    isSingleProcess = this->thread_count == 1;
                    if (isSingleProcess)
                    {

                        this->single_succ();
                    }
                    else
                    {

                        this->parallel_succ();
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
                return true;
            }
            void rep()
            {
                uint64_t finished_level_count = UINT64_MAX;
                assert(this->children_intervals.size() == 0);

                if (this->sub_trees.size() == 0)
                    return;

                if (this->current_lcp == -1)
                {

                    uint64_t _child_count = this->first_child_flag_vec.size();
                    uint64_t _st_count = this->maximal_repeat_check_vec.size();
                    for (auto &it : this->sub_trees)
                    {
                        it->pop_level(this->children_intervals, this->first_child_flag_vec, this->maximal_repeat_check_vec);
                    }
                    this->child_width_vec.push_back(this->first_child_flag_vec.size() - _child_count);
                    this->st_width_vec.push_back(this->maximal_repeat_check_vec.size() - _st_count);
                }
                else
                {
#if DEBUG
                    for (auto &it : this->sub_trees)
                    {
                        assert(it->finished_level_count() > 0);
                    }

#endif

                    for (auto &it : this->sub_trees)
                    {
                        uint64_t maxFinishedLCP = (it->get_current_lcp() - (this->current_lcp + 1)) + it->finished_level_count();
                        if (maxFinishedLCP < finished_level_count)
                        {
                            finished_level_count = maxFinishedLCP;
                        }
                    }
                    assert(finished_level_count != UINT64_MAX);

                    std::queue<uint64_t> tmp_queue;
                    for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                    {
                        tmp_queue.push(i);
                    }
                    for (uint64_t i = 0; i < finished_level_count; i++)
                    {
                        uint64_t k = tmp_queue.size();
                        uint64_t klcp = (this->current_lcp + 1 + i);
                        uint64_t _child_count = this->first_child_flag_vec.size();
                        uint64_t _st_count = this->maximal_repeat_check_vec.size();

                        for (uint64_t j = 0; j < k; j++)
                        {
                            uint64_t top = tmp_queue.front();
                            tmp_queue.pop();
                            auto &it = this->sub_trees[top];
                            if (klcp == it->get_current_lcp())
                            {
                                assert(it->finished_level_count() > 0);
                                it->pop_level(this->children_intervals, this->first_child_flag_vec, this->maximal_repeat_check_vec);

                                tmp_queue.push(top);
                            }
                            else if (klcp < it->get_current_lcp())
                            {
                                tmp_queue.push(top);
                            }
                        }
                        this->child_width_vec.push_back(this->first_child_flag_vec.size() - _child_count);
                        this->st_width_vec.push_back(this->maximal_repeat_check_vec.size() - _st_count);
                    }
                    /*
                    for (auto &it : this->sub_trees)
                    {
                    }
                    */
                    mmax += this->first_child_flag_vec.size();
                }
                assert(this->first_child_flag_vec.size() > 0);
            }
            bool is_finished() const
            {
                return this->current_lcp >= 0 && this->child_count() == 0;
                //return total_counter == strSize - 1;
            }

            void print()
            {
                std::cout << "PRINT PTREE" << std::endl;

                stool::Printer::print_bits("first_child_flag_vec", this->first_child_flag_vec);
                stool::Printer::print_bits("maximal_repeat_check_vec", this->maximal_repeat_check_vec);
                stool::Printer::print("child_width_vec", this->child_width_vec);
                stool::Printer::print("stnode_count_vec", this->st_width_vec);

                /*
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    if (this->sub_trees[i]->get_current_lcp() == (this->current_lcp + 1))
                    {
                        this->sub_trees[i]->print_info();
                    }
                }
                */

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

            void import(STNodeVector<INDEX_SIZE> &item, uint64_t lcp)
            {
                uint64_t k = 1000;
                while (this->sub_trees.size() > 0)
                {
                    delete this->sub_trees[this->sub_trees.size() - 1];
                    this->sub_trees.pop_back();
                }
                while (item.size() > 0)
                {
                    auto st = new STNODE_SUB_TRAVERSER(this->_RLBWTDS);
                    sub_trees.push_back(st);
                    uint64_t num = item.size() < k ? item.size() : k;
                    st->import(item, lcp, num);
                }
                assert(item.size() == 0);

                this->current_lcp = lcp - 1;
            }

        private:
            uint64_t _count_st_size()
            {
                assert(false);
                throw -1;
            }
            uint64_t _count_children_size()
            {
                assert(false);
                throw -1;
            }

            void recompute_node_counter()
            {

                current_lcp++;
            }
            void parallel_succ()
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
                    threads.push_back(thread(parallel_succ_fast_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees), fst_pos_vec[i], ref(position_stack), ref(ems[i]), limit));
                }

                for (thread &t : threads)
                    t.join();
                    //auto end = std::chrono::system_clock::now();
#if DEBUG
                std::cout << "PARALLEL PROCESS[END]" << std::endl;
#endif

                assert(position_stack.size() == 0);
            }
            void single_succ()
            {
                if (this->sub_trees.size() == 0)
                    return;
                uint64_t fst_pos = 0;
                for (uint64_t i = 1; i < this->sub_trees.size(); i++)
                {
                    position_stack.push(i);
                }

                uint64_t limit = peak_mmax - this->mmax;

                parallel_succ_fast_stnodes<INDEX_SIZE, RLBWTDS>(ref(sub_trees), fst_pos, ref(position_stack), ref(ems[0]), limit);
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

        public:
            void increment(ITERATOR &iter) const
            {
                this->increment(iter.child_index, iter.node_index, iter.array_index);
            }
            void increment(INDEX_SIZE &child_index, INDEX_SIZE &node_index, INDEX_SIZE &array_index) const
            {
                RINTERVAL output;
                child_index = this->get_stnode(child_index, output);

                if (child_index >= this->child_width_vec[0] )
                {
                    child_index = std::numeric_limits<INDEX_SIZE>::max();
                    node_index = std::numeric_limits<INDEX_SIZE>::max();
                    array_index = std::numeric_limits<INDEX_SIZE>::max();

                }
                else
                {
                    node_index++;
                }
            }
            INDEX_SIZE get_left(const ITERATOR &iter) const
            {
                return this->get_left(iter.child_index);
            }
            INDEX_SIZE get_left(INDEX_SIZE child_index) const
            {
                RINTERVAL output;
                this->get_stnode(child_index, output);
                uint64_t _left = this->_RLBWTDS->get_lpos(output.beginIndex) + output.beginDiff;
                return _left;
            }

            INDEX_SIZE get_right(const ITERATOR &iter) const
            {
                return this->get_right(iter.child_index);
            }
            INDEX_SIZE get_right(INDEX_SIZE child_index) const
            {
                RINTERVAL output;
                this->get_stnode(child_index, output);
                uint64_t _right = this->_RLBWTDS->get_lpos(output.endIndex) + output.endDiff;
                return _right;
            }

            DEPTH_ITERATOR begin()
            {
                this->clear();
                this->succ();
                return DEPTH_ITERATOR(this, true);
            }
            DEPTH_ITERATOR end()
            {
                return DEPTH_ITERATOR(this, false);
            }
            void set_current_first_iterator(ITERATOR &it) const
            {
                this->set_current_first_iterator(it.child_index, it.node_index, it.array_index);
            }
            void set_current_first_iterator(INDEX_SIZE &child_index, INDEX_SIZE &node_index, INDEX_SIZE &array_index) const
            {
                if (this->is_finished())
                {
                    node_index = std::numeric_limits<INDEX_SIZE>::max();
                    child_index = std::numeric_limits<INDEX_SIZE>::max();
                    array_index = std::numeric_limits<INDEX_SIZE>::max();
                }
                else
                {
                    child_index = 0;
                    node_index = 0;
                    array_index = std::numeric_limits<INDEX_SIZE>::max();
                }
            }

        }; // namespace lcp_on_rlbwt

    } // namespace lcp_on_rlbwt
} // namespace stool