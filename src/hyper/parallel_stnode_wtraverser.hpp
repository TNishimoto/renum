#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>
#include "stnode_wtraverser.hpp"
#include <thread>
#include "stool/src/byte.hpp"
#include <cmath>
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
        void parallel_process_stnodes(std::vector<STNodeWTraverser<INDEX_SIZE, RLBWTDS> *> &trees, uint64_t fst_position,
                                      std::stack<uint64_t> &position_stack, std::vector<STNodeWTraverser<INDEX_SIZE, RLBWTDS> *> &new_trees, uint64_t limit, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
        {
#if DEBUG
            //std::cout << "Run " << data.start_index << ", " << data.width << ", " <<  data.rank << ", " << tree.node_count() << std::endl;
            //tree.print2();
#endif
            uint64_t pos = fst_position;
            bool b = false;
            while (pos != UINT64_MAX)
            {

                bool b2 = trees[pos]->computeNextLCPIntervalSet(em, new_trees, limit);
                b = b || b2;
                assert(trees[pos]->children_count() <= limit);
                /*
                std::cout << "A! " << trees.size() << std::endl;
                if (trees.size() > 1)
                {
                    std::cout << "SP!" << std::endl;
                    throw -1;
                }
                */

                //if (data.width > 0)
                //{
                //for (uint64_t i = data.start_index; i < data.start_index + data.width; i++)
                //{
                //}

                /*
                tmp_tree.maximal_repeat_check_vec.resize(tmp_tree.node_count(), false);

                uint64_t size = tmp_tree.w_builder.size();
                RInterval<INDEX_SIZE> it;
                uint64_t L = 0;
                uint64_t x = 0;
                for (uint64_t i = 1; i < size; i++)
                {
                    if (tmp_tree.w_builder[i])
                    {
                        uint64_t R = i - 1;
                        tmp_tree.get_stnode(L, R, it);
                        bool b = checkMaximalRepeat(it, *(em._RLBWTDS));
                        tmp_tree.maximal_repeat_check_vec[x++] = b;
                        L = i;
                    }
                }
                */
                //}

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
        class ParallelSTNodeWTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using STNODE_WTRAVERSER = STNodeWTraverser<INDEX_SIZE, RLBWTDS>;

            std::vector<STNODE_WTRAVERSER *> sub_trees;
            std::vector<std::vector<STNODE_WTRAVERSER *>> new_trees;

            //std::stack<STNODE_WTRAVERSER> empty_trees;

            //std::vector<STNODE_WTRAVERSER> sub_tmp_trees;
            std::vector<ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> ems;
            std::vector<LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE>> lightRDs;
            std::vector<SuccinctRangeDistinctDataStructure<INDEX_SIZE>> heavyRDs;
            uint64_t minimum_child_count = 1000;
            uint64_t sub_tree_limit_size = 1000;

        public:
            uint64_t current_lcp = 0;
            uint64_t total_counter = 0;
            uint64_t strSize = 0;
            //uint64_t node_count = 0;
            uint64_t child_count = 0;
            uint64_t peak_child_count = 0;
            uint64_t alpha = 2;
            uint64_t _expected_peak_memory_bits = 0;
            uint64_t _switch_threshold = 0;
            uint64_t thread_count = 1;
            std::stack<uint64_t> position_stack;

            RLBWTDS *_RLBWTDS;
            //bool is_parallel = true;
            uint64_t count_maximal_repeats()
            {
                return 0;
                assert(false);

                throw -1;
            }
            /*
            STNODE_WTRAVERSER *get_sub_tree()
            {
                return &this->sub_tree;
            }
            */
            uint64_t node_count() const
            {
                uint64_t k = 0;
                for (auto &it : this->sub_trees)
                {
                    k += it->node_count();
                }
                return k;

                //return this->sub_tree.node_count();
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
                this->strSize = _RLBWTDS.str_size();
                this->thread_count = size;

                //assert(size == 1);

                //this->sub_trees.resize(size);

                auto st = new STNODE_WTRAVERSER(this->_RLBWTDS);
                sub_trees.push_back(st);
                //sub_trees.resize(256);

                this->new_trees.resize(size);

                for (uint64_t i = 0; i < this->thread_count; i++)
                {

                    lightRDs.push_back(LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE>());
                    lightRDs[lightRDs.size() - 1].preprocess(&_RLBWTDS.bwt);

                    heavyRDs.push_back(SuccinctRangeDistinctDataStructure<INDEX_SIZE>());
                    heavyRDs[heavyRDs.size() - 1].initialize(&_RLBWTDS.wt, &_RLBWTDS.bwt);

                    ems.push_back(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>());
                    ems[ems.size() - 1].initialize(&_RLBWTDS);
                }
                for (uint64_t i = 0; i < this->thread_count; i++)
                {

                    ems[i].lightDS = &lightRDs[i];
                    ems[i].heavyDS = &heavyRDs[i];
                }
                double ratio = (double)this->_RLBWTDS->str_size() / (double)this->_RLBWTDS->rle_size();
                double d = std::log2(ratio);
                this->_expected_peak_memory_bits = this->_RLBWTDS->rle_size() * d;
                this->_switch_threshold = this->alpha * (this->expected_peak_memory_bits() / (sizeof(uint64_t) * 8));
            }

            void get_stnode(uint64_t index, RINTERVAL &output)
            {
                std::cout << "get ST NODE" << std::endl;
                throw -1;
                /*
                uint64_t p = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    if (p <= index && index < p + this->sub_trees[i].node_count())
                    {
                        this->sub_trees[i].get_stnode(index - p, output);

                        break;
                    }
                    else
                    {
                        p += this->sub_trees[i].node_count();
                    }
                }
                */
            }
            uint64_t get_stnode2(uint64_t L, RINTERVAL &output)
            {

                uint64_t p = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    uint64_t csize = this->sub_trees[i]->children_count();

                    if (p <= L && L < p + csize)
                    {
                        return p + this->sub_trees[i]->get_stnode2(L - p, output);
                    }
                    else
                    {
                        p += csize;
                    }
                }
                return UINT64_MAX;
            }

            /*
            RINTERVAL &get_child(uint64_t index)
            {
                uint64_t x = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    uint64_t k = this->sub_trees[i].weinerCount;
                    if (x <= index && index < x + k)
                    {
                        return this->sub_trees[i].weinerVec[index - x];
                    }
                    else
                    {
                        x += k;
                    }
                }
                assert(false);
                throw -1;
            }
            */
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
                /*
                if (mergeCount > 3)
                {
                    std::cout << "M" << mergeCount << std::endl;
                }
                */
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

                            //this->sub_trees[i] = STNODE_WTRAVERSER();
                            //sub_trees[i]._RLBWTDS = this->_RLBWTDS;
                        }
                        nonEmptyCount++;
                    }
                }
                for(uint64_t i=nonEmptyCount;i<this->sub_trees.size(); i++){
                    delete this->sub_trees[i];
                }
                this->sub_trees.resize(nonEmptyCount);
            }
            /*
            void split_big_trees()
            {

                uint64_t xindex = 0;
                while (xindex < this->sub_trees.size())
                {
                    if (this->sub_trees[xindex].children_count() > this->sub_tree_limit_size)
                    {
                        std::cout << "S" << this->sub_trees[xindex].children_count() << "/" << this->sub_tree_limit_size << std::endl;

                        assert(false);
                        this->sub_trees.push_back(STNODE_WTRAVERSER());
                        auto &it = this->sub_trees[sub_trees.size() - 1];
                        it._RLBWTDS = this->_RLBWTDS;
                        this->sub_trees[xindex].split(it);
                    }
                    else
                    {
                        xindex++;
                    }
                }
            }
            */

            uint64_t kk = 0;
#if DEBUG
            uint64_t prev_child_count = 0;
#endif
            void process()
            {
                bool isSingleProcess = false;
                if (current_lcp > 0)
                {
                    
                    if (this->child_count > this->switch_threshold())
                    {
                        kk++;
                        std::cout << "LCP " << this->current_lcp << "/" << (this->child_count) << "/" << this->switch_threshold() << "/" << this->sub_trees.size() << std::endl;
                        std::cout << "Memory: " << this->get_peak_memory() / (1000 * 1000) << "[MB]"
                                  << "/ Optimal: " << this->get_optimal_memory() / (1000 * 1000) << "[MB]" << std::endl;
                        /*
                        if (kk == 3)
                        {
                            std::cout << "STOP" << std::endl;
                            throw -1;
                        }
                        */
                        
                    }
                    
                    

#if DEBUG
                    if (prev_child_count != this->child_count)
                    {
                        std::cout << "LCP = " << current_lcp << "/" << this->child_count << "/" << this->sub_trees.size() << std::endl;
                        prev_child_count = this->child_count;
                    }
                    if (this->_RLBWTDS->str_size() < 100)
                    {
                        this->print();
                    }
#endif


                    isSingleProcess = this->child_count < minimum_child_count || this->thread_count == 1;
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

                uint64_t kp = 0;
                for (auto &it : this->new_trees)
                {
                    while (it.size() > 0)
                    {

                        this->sub_trees.push_back(it[it.size() - 1]);
                        it.pop_back();
                        kp++;
                    }
                }

                this->merge();
                this->remove_empty_trees();

                if ((double)this->sub_trees.size() * 2 < (double)this->sub_trees.capacity())
                {
                    this->sub_trees.shrink_to_fit();
                }

                uint64_t current_child_count = 0;
                uint64_t current_node_count = 0;
                for (auto &it : this->sub_trees)
                {
                    it->build_bits();

                    current_child_count += it->children_count();
                    current_node_count += it->node_count();
                }

                total_counter += (current_child_count - current_node_count);
                this->child_count = current_child_count;
                //this->node_count = current_node_count;
                if (current_child_count > this->peak_child_count)
                {
                    this->peak_child_count = current_child_count;
                }
#if DEBUG

                if (this->_RLBWTDS->str_size() < 100)
                {
                    this->print();
                }

#endif

                assert(total_counter <= strSize);

                current_lcp++;
                assert(this->child_count > 0);
                assert(this->node_count() > 0);
            }
            bool isStop()
            {
                return total_counter == strSize - 1;
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

        private:
            void single_process()
            {
                uint64_t fst_pos = 0;
                for (uint64_t i = 1; i < this->sub_trees.size(); i++)
                {
                    position_stack.push(i);
                }

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

                std::vector<thread> threads;
                for (uint64_t i = 0; i < this->thread_count; i++)
                {
                    threads.push_back(thread(parallel_process_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees), fst_pos_vec[i], ref(position_stack), ref(new_trees[i]), this->sub_tree_limit_size, ref(ems[i])));
                }
                for (thread &t : threads)
                    t.join();

                assert(position_stack.size() == 0);
            }
            /*
            void get_start_indexes(std::vector<ParallelData> &output)
            {
                assert(false);

                throw -1;
                
                for (uint64_t i = 0; i < output.size(); i++)
                {
                    output[i].width = 0;
                }
                uint64_t child_avg_count = (this->child_count / this->sub_tmp_trees.size()) + 1;
                uint64_t start_index = 0;
                uint64_t child_start_index = 0;
                uint64_t p = 0;
                output[p] = ParallelData(start_index, child_start_index);
                uint64_t k = 0;
                for (uint64_t i = 0; i < this->node_count(); i++)
                {

                    output[p].width++;
                    k += this->sub_tree.get_width(i);
                    start_index++;
                    if (k >= child_avg_count)
                    {
                        child_start_index += k;
                        k = 0;
                        p++;
                        if (p < output.size())
                        {
                            output[p] = ParallelData(start_index, child_start_index);
                        }
                        else
                        {
                            std::cout << "Child count = " << this->child_count << ", size = " << this->sub_tmp_trees.size() << ", avg = " << child_avg_count << std::endl;
                            assert(i + 1 == this->sub_tree.node_count());
                        }
                    }
                }
                
            }
            */

            uint64_t get_peak_memory()
            {
                uint64_t k = 0;

                for (auto &it : this->sub_trees)
                {
                    k += it->get_peak_memory();
                }
                /*
                for (auto &it : this->sub_tmp_trees)
                {
                    k += it.get_peak_memory();
                }
                */

                return k;
            }
            uint64_t get_optimal_memory()
            {
                uint64_t k = 0;

                for (auto &it : this->sub_trees)
                {
                    k += it->get_optimal_memory();
                }
                /*
                for (auto &it : this->sub_tmp_trees)
                {
                    k += it.get_optimal_memory();
                }
                */

                return k;
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool