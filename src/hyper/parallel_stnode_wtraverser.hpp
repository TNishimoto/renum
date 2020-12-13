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
        void parallel_process_stnodes(std::vector<STNodeWTraverser<INDEX_SIZE, RLBWTDS>> &trees, uint64_t fst_position, std::stack<uint64_t> &position_stack, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
        {
#if DEBUG
            //std::cout << "Run " << data.start_index << ", " << data.width << ", " <<  data.rank << ", " << tree.node_count() << std::endl;
            //tree.print2();
#endif
            /*
            assert(data.start_index < tree.node_count());
            assert(data.start_index + data.width - 1 < tree.node_count());
            */

            uint64_t pos = fst_position;
            while (pos != UINT64_MAX)
            {

                trees[pos].computeNextLCPIntervalSet(em);

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

            std::vector<STNODE_WTRAVERSER> sub_trees;
            //std::vector<STNODE_WTRAVERSER> sub_tmp_trees;
            std::vector<ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> ems;
            std::vector<LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE>> lightRDs;
            std::vector<SuccinctRangeDistinctDataStructure<INDEX_SIZE>> heavyRDs;
            uint64_t minimum_child_count = 1000;
            uint64_t sub_tree_limit_size = 10000;

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
            bool is_parallel = true;
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
                    k += it.node_count();
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

                sub_trees.push_back(STNODE_WTRAVERSER());
                sub_trees[sub_trees.size() - 1]._RLBWTDS = &_RLBWTDS;

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

                //return this->sub_tree.get_stnode(index, output);
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

            void process()
            {
                if (current_lcp > 0)
                {
                    if (this->child_count > this->switch_threshold())
                    {
                        std::cout << "LCP " << (this->child_count) << "/" << this->switch_threshold() << std::endl;
                    }

#if DEBUG
                    //std::cout << "LCP " << current_lcp << "/" << this->child_count << std::endl;

                    if (this->_RLBWTDS->str_size() < 100)
                    {
                        this->print();
                    }

#endif
                    bool b = (this->child_count < minimum_child_count || this->sub_trees.size() == 1) || !this->is_parallel;
                    //bool b = true;
                    if (b)
                    {
                        this->single_process();
                    }
                    else
                    {
                        //this->allocate_data();
                        //auto start = std::chrono::system_clock::now();

                        this->parallel_process();

                        //auto end = std::chrono::system_clock::now();
                        //double elapsed1 = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
                        //std::cout << "Time, " << elapsed1 << ", " << std::endl;
                    }
                }
                else
                {

                    this->sub_trees[0].first_compute(ems[0]);
                }
                /*
                this->sub_tree.clear();
                for (uint64_t i = 0; i < this->sub_tmp_trees.size(); i++)
                {
                    this->sub_tree.add(this->sub_tmp_trees[i]);
                }
                */

                uint64_t current_child_count = 0;
                uint64_t current_node_count = 0;

                uint64_t xindex = 0;
                while (xindex < this->sub_trees.size())
                {
                    if (this->sub_trees[xindex].children_count() > this->sub_tree_limit_size)
                    {
                        //std::cout << "SPLIT" << this->sub_trees.size() << std::endl;
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
                uint64_t nonEmptyCount = 0;
                for(uint64_t i=0;i<this->sub_trees.size();i++){
                    if(this->sub_trees[i].node_count() > 0){
                        if(i != nonEmptyCount){
                            this->sub_trees[nonEmptyCount].swap(this->sub_trees[i]);
                        }
                        nonEmptyCount++;
                    }
                }
                if(this->sub_trees.size() != nonEmptyCount){
                    std::cout << "POP << " << this->sub_trees.size() << "/" << nonEmptyCount <<std::endl;
                }
                this->sub_trees.resize(nonEmptyCount);

                for (auto &it : this->sub_trees)
                {
                    it.build_bits();

                    current_child_count += it.children_count();
                    current_node_count += it.node_count();
                }

                total_counter += (current_child_count - current_node_count);
                this->child_count = current_child_count;
                //this->node_count = current_node_count;
                if (current_child_count > this->peak_child_count)
                {
                    this->peak_child_count = current_child_count;
                }

                //std::cout << "Memory: " << this->get_peak_memory() / (1000 * 1000) << "[MB]" << "/ Optimal: " << this->get_optimal_memory() / (1000 * 1000) << "[MB]" << std::endl;

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
                /*
                for (uint64_t i = 0; i < this->sub_tmp_trees.size(); i++)
                {
                    this->sub_tmp_trees[i].print_info();
                }
                */
                std::cout << std::endl;
                std::cout << "PRINT PTREE[END]" << std::endl;
            }

        private:
            void single_process()
            {
                /*
                ParallelData pd;
                pd.start_index = 0;
                pd.width = sub_trees.size();
                */
                uint64_t fst_pos = 0;
                for (uint64_t i = 1; i < this->sub_trees.size(); i++)
                {
                    position_stack.push(i);
                }

                parallel_process_stnodes<INDEX_SIZE, RLBWTDS>(ref(sub_trees), fst_pos, ref(position_stack), ref(ems[0]));
                assert(position_stack.size() == 0);
                /*
                std::vector<ParallelData> w_vec;

                w_vec.resize(this->sub_tmp_trees.size());

                this->get_start_indexes(w_vec);
                

                for (uint64_t i = 0; i < this->sub_tmp_trees.size(); i++)
                {
                    if (w_vec[i].width > 0)
                    {
                    }
                }
                */
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
                    if(i < this->thread_count){
                        fst_pos_vec[i] = i;
                    }else{
                        position_stack.push(i);
                    }
                }

                std::vector<thread> threads;
                for (uint64_t i = 0; i < this->thread_count; i++)
                {
                        threads.push_back(thread(parallel_process_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees), fst_pos_vec[i], ref(position_stack), ref(ems[i])));

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
                    k += it.get_peak_memory();
                }

                for (auto &it : this->sub_tmp_trees)
                {
                    k += it.get_peak_memory();
                }

                return k;
            }
            uint64_t get_optimal_memory()
            {
                uint64_t k = 0;

                for (auto &it : this->sub_trees)
                {
                    k += it.get_optimal_memory();
                }

                for (auto &it : this->sub_tmp_trees)
                {
                    k += it.get_optimal_memory();
                }

                return k;
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool