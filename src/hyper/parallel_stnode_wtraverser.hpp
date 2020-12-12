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
        template <typename INDEX_SIZE, typename RLBWTDS>
        void parallel_process_stnodes(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &tree, STNodeWTraverser<INDEX_SIZE, RLBWTDS> &tmp_tree, ParallelData data, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
        {
#if DEBUG
            //std::cout << "Run " << data.start_index << ", " << data.width << ", " <<  data.rank << ", " << tree.node_count() << std::endl;
            //tree.print2();
#endif
            assert(data.start_index < tree.node_count());
            assert(data.start_index + data.width - 1 < tree.node_count());

            if (data.width > 0)
            {
                tmp_tree.computeNextLCPIntervalSet(tree, data.start_index, data.width, data.rank, em);
                //tree.swap(tmp_tree);

                tmp_tree.maximal_repeat_check_vec.resize(tmp_tree.node_count(), false);

                uint64_t size = tmp_tree.w_builder.size();
                RInterval<INDEX_SIZE> it;
                uint64_t L = 0;
                uint64_t x=0;
                for (uint64_t i = 1; i < size; i++)
                {
                    if(tmp_tree.w_builder[i]){
                        uint64_t R = i-1;
                        tmp_tree.get_stnode(L, R, it);
                        bool b = checkMaximalRepeat(it, *(em._RLBWTDS));
                        tmp_tree.maximal_repeat_check_vec[x++] = b;
                        L = i;
                    }
                }
            }
        }

        template <typename INDEX_SIZE, typename RLBWTDS>
        class ParallelSTNodeWTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using STNODE_WTRAVERSER = STNodeWTraverser<INDEX_SIZE, RLBWTDS>;

            /*
            std::vector<RINTERVAL> stnodeVec;
            std::vector<RINTERVAL> childVec;
            std::vector<uint8_t> widthVec;
            std::vector<bool> maximal_repeat_check_vec;
            */


            STNODE_WTRAVERSER sub_tree;
            std::vector<STNODE_WTRAVERSER> sub_tmp_trees;
            std::vector<ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> ems;
            std::vector<LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE>> lightRDs;
            std::vector<SuccinctRangeDistinctDataStructure<INDEX_SIZE>> heavyRDs;
            uint64_t minimum_child_count = 1000;

        public:
            uint64_t current_lcp = 0;
            uint64_t total_counter = 0;
            uint64_t strSize = 0;
            //uint64_t node_count = 0;
            uint64_t child_count = 0;
            uint64_t peak_child_count = 0;
            RLBWTDS *_RLBWTDS;
            bool is_parallel = true;

            STNODE_WTRAVERSER *get_sub_tree()
            {
                return &this->sub_tree;
            }
            uint64_t node_count() const {
                return this->sub_tree.node_count();
            }

            void initialize(uint64_t size, RLBWTDS &_RLBWTDS)
            {
                this->_RLBWTDS = &_RLBWTDS;
                this->strSize = _RLBWTDS.str_size();

                sub_tree._RLBWTDS = &_RLBWTDS;

                for (uint64_t i = 0; i < size; i++)
                {
                    //sub_trees.push_back(STNODE_WTRAVERSER());
                    sub_tmp_trees.push_back(STNODE_WTRAVERSER());
                    sub_tmp_trees[sub_tmp_trees.size() - 1]._RLBWTDS = &_RLBWTDS;

                    lightRDs.push_back(LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE>());
                    lightRDs[lightRDs.size() - 1].preprocess(&_RLBWTDS.bwt);

                    heavyRDs.push_back(SuccinctRangeDistinctDataStructure<INDEX_SIZE>());
                    heavyRDs[heavyRDs.size() - 1].initialize(&_RLBWTDS.wt, &_RLBWTDS.bwt);

                    ems.push_back(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>());
                    ems[ems.size() - 1].initialize(&_RLBWTDS);
                }
                for (uint64_t i = 0; i < this->sub_tmp_trees.size(); i++)
                {

                    ems[i].lightDS = &lightRDs[i];
                    ems[i].heavyDS = &heavyRDs[i];
                }
            }
            void get_stnode(uint64_t index, RINTERVAL &output)
            {
                return this->sub_tree.get_stnode(index, output);

            }
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

            void process()
            {
                if (current_lcp > 0)
                {
                    if(current_lcp % 100 == 0){
                    std::cout << "LCP " << (strSize - total_counter) << std::endl;

                    }

#if DEBUG
                                        std::cout << "LCP " << current_lcp << "/" << this->child_count << std::endl;

                    if (this->_RLBWTDS->str_size() < 100)
                    {
                        this->print();
                    }

#endif
                    bool b = (this->child_count < minimum_child_count || this->sub_tmp_trees.size() == 1) || !this->is_parallel;
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

                    this->sub_tmp_trees[0].first_compute(ems[0]);
                }
                this->sub_tree.clear();
                for (uint64_t i = 0; i < this->sub_tmp_trees.size(); i++)
                {
                    this->sub_tree.add(this->sub_tmp_trees[i]);
                }
                this->sub_tree.build_bits();

                uint64_t current_child_count = sub_tree.children_count();
                uint64_t current_node_count = sub_tree.node_count();

                total_counter += (current_child_count - current_node_count);
                this->child_count = current_child_count;
                //this->node_count = current_node_count;
                if (current_child_count > this->peak_child_count)
                {
                    this->peak_child_count = current_child_count;
                }
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
                this->sub_tree.print_info();
                for (uint64_t i = 0; i < this->sub_tmp_trees.size(); i++)
                {
                    this->sub_tmp_trees[i].print_info();
                }
                std::cout << std::endl;
                std::cout << "PRINT PTREE[END]" << std::endl;

            }
            
        private:
            void single_process()
            {
                std::vector<ParallelData> w_vec;

                w_vec.resize(this->sub_tmp_trees.size());

                this->get_start_indexes(w_vec);
                for (uint64_t i = 0; i < this->sub_tmp_trees.size(); i++)
                {
                    if (w_vec[i].width > 0)
                    {
                        parallel_process_stnodes<INDEX_SIZE, RLBWTDS>(ref(sub_tree), ref(sub_tmp_trees[i]), w_vec[i], ref(ems[i]));
                    }
                }
                /*
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    parallel_process_stnodes<INDEX_SIZE, RLBWTDS>(sub_trees[i], sub_tmp_trees[i], ems[i]);
                }
                */
            }
            void parallel_process()
            {
                #if DEBUG
                std::cout << "PARALLEL PROCESS" << std::endl;
                #endif
                std::vector<ParallelData> w_vec;
                w_vec.resize(this->sub_tmp_trees.size());
                this->get_start_indexes(w_vec);
                std::vector<thread> threads;
                for (uint64_t i = 0; i < this->sub_tmp_trees.size(); i++)
                {
                    if (w_vec[i].width > 0)
                    {
                        threads.push_back(thread(parallel_process_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_tree), ref(sub_tmp_trees[i]), w_vec[i], ref(ems[i])));

                    }

                }
                for (thread &t : threads)
                    t.join();
            }
            void get_start_indexes(std::vector<ParallelData> &output)
            {
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
                for (uint64_t i = 0; i < this->sub_tree.node_count(); i++)
                {

                    output[p].width++;
                    k += this->sub_tree.get_width(i);
                    start_index++;
                    if (k >= child_avg_count)
                    {
                        child_start_index += k;
                        k = 0;
                        p++;
                        if(p < output.size()){
                            output[p] = ParallelData(start_index, child_start_index);
                        }else{
                            std::cout << "Child count = " << this->child_count << ", size = " << this->sub_tmp_trees.size() << ", avg = " << child_avg_count << std::endl;
                            assert(i + 1 == this->sub_tree.node_count());
                        }

                    }
                }
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool