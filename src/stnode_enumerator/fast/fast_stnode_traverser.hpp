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
                                           std::stack<uint64_t> &position_stack, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
        {

            uint64_t pos = fst_position;
            while (pos != UINT64_MAX)
            {

                trees[pos]->computeNextLCPIntervalSet(em, 1000);
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

            std::vector<ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> ems;
            uint64_t minimum_child_count = 1000;
            uint64_t sub_tree_limit_size = 2000;

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

                if (this->thread_count == 1)
                {
                    this->sub_tree_limit_size = UINT64_MAX;
                }

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
                for (auto &it : this->sub_trees)
                {
                    uint64_t size = it->current_node_count();
                    uint64_t L = 0;
                    RINTERVAL intv;
                    for (uint64_t i = 0; i < size; i++)
                    {
                        L = it->read_st_node(L, intv);
                        if (it->check_maximal_repeat(i))
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

                uint64_t p = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    uint64_t csize = this->sub_trees[i]->current_children_count();

                    if (p <= L && L < p + csize)
                    {
                        return p + this->sub_trees[i]->read_st_node(L - p, output);
                    }
                    else
                    {
                        p += csize;
                    }
                }
                return UINT64_MAX;
            }

            void process()
            {
                bool isSingleProcess = false;

                if (current_lcp > 0)
                {

                    isSingleProcess = this->child_count() < minimum_child_count || this->thread_count == 1;
                    //bool b = true;
                    if (isSingleProcess)
                    {
                        this->single_process();
                    }
                    else
                    {

                        this->parallel_process();
                    }
                }
                else
                {

                    this->sub_trees[0]->first_compute(ems[0]);
                }

                //this->merge();
                this->remove_empty_trees();

                if ((double)this->sub_trees.size() * 2 < (double)this->sub_trees.capacity())
                {
                    this->sub_trees.shrink_to_fit();
                }
                this->recompute_node_counter();
            }
            bool isStop()
            {
                return this->current_lcp > 0 && this->child_count() == 0;
                //return total_counter == strSize - 1;
            }

            void print()
            {
                std::cout << "PRINT PTREE" << std::endl;
                //this->sub_tree.print_info();

                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    this->sub_trees[i]->print_info();
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

        private:
            void recompute_node_counter()
            {

                uint64_t current_child_count = 0;
                uint64_t current_node_count = 0;

                for (auto &it : this->sub_trees)
                {
                    while(it->get_current_lcp() < this->current_lcp){
                        it->pop_level();
                    }
                    current_child_count += it->current_children_count();
                    current_node_count += it->current_node_count();
                }

                _node_count = current_node_count;
                _child_count = current_child_count;

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
                std::vector<thread> threads;
                for (uint64_t i = 0; i < this->thread_count; i++)
                {
                    threads.push_back(thread(parallel_process_fast_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees), fst_pos_vec[i], ref(position_stack), ref(ems[i])));
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
                if (this->child_count() == 0)
                    return;
                uint64_t fst_pos = 0;
                for (uint64_t i = 1; i < this->sub_trees.size(); i++)
                {
                    position_stack.push(i);
                }

                assert(this->child_count() > 0);

                parallel_process_fast_stnodes<INDEX_SIZE, RLBWTDS>(ref(sub_trees), fst_pos, ref(position_stack), ref(ems[0]));
                assert(position_stack.size() == 0);
            }
            /*
            void merge()
            {

                uint64_t mergeIndex = 0;
                uint64_t mergeCount = 0;

                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    uint64_t k = this->sub_trees[i]->children_count();
                    if (k < (this->sub_tree_limit_size / 10))
                    {
                        while (mergeIndex < i)
                        {
                            if (this->sub_trees[mergeIndex]->children_count() > 0 && this->sub_trees[mergeIndex]->children_count() + k < this->sub_tree_limit_size)
                            {
                                mergeCount++;
                                this->sub_trees[mergeIndex]->merge(*this->sub_trees[i]);
                                break;
                            }
                            else
                            {
                                mergeIndex++;
                            }
                        }
                    }
                }
            }
            */
            void remove_empty_trees()
            {
                uint64_t nonEmptyCount = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    if (this->sub_trees[i]->current_node_count() > 0)
                    {
                        if (i != nonEmptyCount)
                        {

                            this->sub_trees[nonEmptyCount]->swap(*this->sub_trees[i]);

                            //this->sub_trees[i] = STNODE_SUB_TRAVERSER();
                            //sub_trees[i]._RLBWTDS = this->_RLBWTDS;
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
        };

    } // namespace lcp_on_rlbwt
} // namespace stool