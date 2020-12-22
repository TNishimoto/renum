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
#include "stnode_sub_wtraverser.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        struct ParallelData
        {
            uint64_t start_index;
            uint64_t width;
            uint64_t rank;

            ParallelData()
            {
            }
            ParallelData(uint64_t _start_index, uint64_t _rank)
            {
                this->start_index = _start_index;
                this->rank = _rank;
                this->width = 0;
            }
        };

        template <typename INDEX_SIZE, typename RLBWTDS>
        bool checkMaximalRepeat(const RInterval<INDEX_SIZE> &lcpIntv, RLBWTDS &_RLBWTDS)
        {
            RInterval<INDEX_SIZE> it = _RLBWTDS.getIntervalOnL(lcpIntv);
            uint8_t fstChar = _RLBWTDS.get_char_by_run_index(it.beginIndex);
            uint8_t lstChar = _RLBWTDS.get_char_by_run_index(it.endIndex);
            if (fstChar == lstChar)
            {

                if (it.beginIndex != it.endIndex)
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }
            else
            {
                return true;
            }
        }
        std::mutex mtx;
        template <typename INDEX_SIZE, typename RLBWTDS>
        void parallel_process_stnodes(std::vector<STNodeSubWTraverser<INDEX_SIZE, RLBWTDS> *> &trees, uint64_t fst_position,
                                      std::stack<uint64_t> &position_stack, std::vector<STNodeSubWTraverser<INDEX_SIZE, RLBWTDS> *> &new_trees, uint64_t limit, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
        {
            uint64_t pos = fst_position;
            bool b = false;
            while (pos != UINT64_MAX)
            {

                bool b2 = trees[pos]->computeNextLCPIntervalSet(em, new_trees, limit);
                b = b || b2;
                assert(trees[pos]->children_count() <= limit);
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
        class STNodeWTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using STNODE_SUB_WTRAVERSER = STNodeSubWTraverser<INDEX_SIZE, RLBWTDS>;

            std::vector<STNODE_SUB_WTRAVERSER *> sub_trees;
            std::vector<std::vector<STNODE_SUB_WTRAVERSER *>> new_trees;

            std::vector<ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> ems;
            std::vector<LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE>> lightRDs;
            std::vector<SuccinctRangeDistinctDataStructure<INDEX_SIZE>> heavyRDs;
            uint64_t minimum_child_count = 1000;
            uint64_t sub_tree_limit_size = 2000;


            uint64_t current_lcp = 0;
            uint64_t _child_count = 0;
            uint64_t _node_count = 0;

            uint64_t thread_count = 1;
            std::stack<uint64_t> position_stack;


        public:
            RLBWTDS *_RLBWTDS;

            uint64_t get_current_lcp() const {
                return current_lcp;
            }
            uint64_t child_count() const {
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

                //assert(size == 1);

                //this->sub_trees.resize(size);

                auto st = new STNODE_SUB_WTRAVERSER(this->_RLBWTDS);
                sub_trees.push_back(st);
                //sub_trees.resize(256);

                this->new_trees.resize(size);

                uint8_t lastChar = __RLBWTDS.bwt[__RLBWTDS.bwt.size() - 1];

                for (uint64_t i = 0; i < this->thread_count; i++)
                {

                    lightRDs.push_back(LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE>());
                    lightRDs[lightRDs.size() - 1].preprocess(&__RLBWTDS.bwt);

                    heavyRDs.push_back(SuccinctRangeDistinctDataStructure<INDEX_SIZE>());
                    heavyRDs[heavyRDs.size() - 1].initialize(&__RLBWTDS.wt, lastChar);

                    ems.push_back(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>());
                    ems[ems.size() - 1].initialize(&__RLBWTDS);
                }
                for (uint64_t i = 0; i < this->thread_count; i++)
                {

                    ems[i].lightDS = &lightRDs[i];
                    ems[i].heavyDS = &heavyRDs[i];
                }
                /*
                double ratio = (double)this->_RLBWTDS->str_size() / (double)this->_RLBWTDS->rle_size();
                double d = std::log2(ratio);
                this->_expected_peak_memory_bits = this->_RLBWTDS->rle_size() * d;
                this->_switch_threshold = this->alpha * (this->expected_peak_memory_bits() / (sizeof(uint64_t) * 8));
                */
            }
            uint64_t write_maximal_repeats(std::ofstream &out)
            {
                uint64_t count = 0;
                std::vector<stool::LCPInterval<INDEX_SIZE>> buffer;
                for (auto &it : this->sub_trees)
                {
                    uint64_t size = it->node_count();
                    uint64_t L = 0;
                    uint64_t _left = 0, _right = 0;
                    for (uint64_t i = 0; i < size; i++)
                    {
                        L = it->increment(L, _left, _right);
                        if (it->check_maximal_repeat(i))
                        {
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
            uint64_t get_stnode2(uint64_t L, stool::LCPInterval<uint64_t> &output)
            {
                uint64_t p = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    uint64_t csize = this->sub_trees[i]->children_count();

                    if (p <= L && L < p + csize)
                    {
                        return p + this->sub_trees[i]->get_stnode2(L - p, output, this->current_lcp - 1);
                    }
                    else
                    {
                        p += csize;
                    }
                }
                return UINT64_MAX;
            }
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
            void remove_empty_trees()
            {
                uint64_t nonEmptyCount = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    if (this->sub_trees[i]->node_count() > 0)
                    {
                        if (i != nonEmptyCount)
                        {

                            this->sub_trees[nonEmptyCount]->swap(*this->sub_trees[i]);

                            //this->sub_trees[i] = STNODE_SUB_WTRAVERSER();
                            //sub_trees[i]._RLBWTDS = this->_RLBWTDS;
                        }
                        nonEmptyCount++;
                    }
                }
                for (uint64_t i = nonEmptyCount; i < this->sub_trees.size(); i++)
                {
                    this->sub_trees[i]->clear();
                    delete this->sub_trees[i];
                }
                this->sub_trees.resize(nonEmptyCount);
            }

            uint64_t kk = 0;
#if DEBUG
            uint64_t prev_child_count = 0;
#endif
            void heavyEnumerate()
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

                for (auto &it : this->new_trees)
                {
                    while (it.size() > 0)
                    {

                        this->sub_trees.push_back(it[it.size() - 1]);
                        it.pop_back();
                    }
                }

                this->merge();
                this->remove_empty_trees();

                if ((double)this->sub_trees.size() * 2 < (double)this->sub_trees.capacity())
                {
                    this->sub_trees.shrink_to_fit();
                }



            }

            void process()
            {
                this->heavyEnumerate();
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

            void single_process()
            {
                if(this->child_count() == 0) return;
                uint64_t fst_pos = 0;
                for (uint64_t i = 1; i < this->sub_trees.size(); i++)
                {
                    position_stack.push(i);
                }

                assert(this->child_count() > 0);

                parallel_process_stnodes<INDEX_SIZE, RLBWTDS>(ref(sub_trees), fst_pos, ref(position_stack), ref(new_trees[0]), this->sub_tree_limit_size, ref(ems[0]));
                assert(position_stack.size() == 0);

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
                    threads.push_back(thread(parallel_process_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees), fst_pos_vec[i], ref(position_stack), ref(new_trees[i]), this->sub_tree_limit_size, ref(ems[i])));
                }

                for (thread &t : threads)
                    t.join();
                //auto end = std::chrono::system_clock::now();

                assert(position_stack.size() == 0);
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
            void recompute_node_counter()
            {
                
                uint64_t current_child_count = 0;
                uint64_t current_node_count = 0;
                for (auto &it : this->sub_trees)
                {
                    current_child_count += it->children_count();
                    current_node_count += it->node_count();
                }
                _node_count = current_node_count;
                _child_count = current_child_count;

                //total_counter += (current_child_count - current_node_count);
                //this->child_count = current_child_count;
                //this->node_count = current_node_count;
                /*
                if (current_child_count > this->peak_child_count)
                {
                    this->peak_child_count = current_child_count;
                }
                */
                current_lcp++;
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool