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
        void parallel_process_stnodes(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &tree, STNodeWTraverser<INDEX_SIZE, RLBWTDS> &tmp_tree, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
        {

            if (tree.node_count() > 0)
            {
                tmp_tree.computeNextLCPIntervalSet(tree, em);
                tree.swap(tmp_tree);

                tree.maximal_repeat_check_vec.resize(tree.node_count(), false);
                
                uint64_t size = tree.node_count();
                for (uint64_t i = 0; i < size; i++)
                {
                    const RInterval<INDEX_SIZE> &it = tree.get_stnode(i);
                    bool b = checkMaximalRepeat(it, *(em._RLBWTDS));
                    tree.maximal_repeat_check_vec[i] = b;
                }
            
            }
        }

        template <typename INDEX_SIZE, typename RLBWTDS>
        class ParallelSTNodeWTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using STNODE_WTRAVERSER = STNodeWTraverser<INDEX_SIZE, RLBWTDS>;

            std::vector<STNODE_WTRAVERSER> sub_trees;
            std::vector<STNODE_WTRAVERSER> sub_tmp_trees;
            std::vector<ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> ems;
            std::vector<LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE>> lightRDs;
            std::vector<SuccinctRangeDistinctDataStructure<INDEX_SIZE>> heavyRDs;

        public:
            uint64_t current_lcp = 0;
            uint64_t total_counter = 0;
            uint64_t strSize = 0;
            uint64_t node_count = 0;
            uint64_t child_count = 0;
            uint64_t peak_child_count = 0;
            RLBWTDS *_RLBWTDS;

            std::vector<STNODE_WTRAVERSER> *get_sub_trees()
            {
                return &this->sub_trees;
            }
            uint64_t get_tree_count()
            {
                return this->sub_trees.size();
            }

            void initialize(uint64_t size, RLBWTDS &_RLBWTDS)
            {
                this->_RLBWTDS = &_RLBWTDS;
                this->strSize = _RLBWTDS.str_size();

                for (uint64_t i = 0; i < size; i++)
                {
                    sub_trees.push_back(STNODE_WTRAVERSER());
                    sub_tmp_trees.push_back(STNODE_WTRAVERSER());
                    lightRDs.push_back(LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE>());
                    lightRDs[lightRDs.size() - 1].preprocess(&_RLBWTDS.bwt);

                    heavyRDs.push_back(SuccinctRangeDistinctDataStructure<INDEX_SIZE>());
                    heavyRDs[heavyRDs.size() - 1].initialize(&_RLBWTDS.wt, &_RLBWTDS.bwt);

                    ems.push_back(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>());
                    ems[ems.size() - 1].initialize(&_RLBWTDS);
                }
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {

                    ems[i].lightDS = &lightRDs[i];
                    ems[i].heavyDS = &heavyRDs[i];
                }
            }
            RINTERVAL &get_stnode(uint64_t index)
            {
                uint64_t x = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    uint64_t k = this->sub_trees[i].node_count();
                    if (x <= index && index < x + k)
                    {
                        return this->sub_trees[i].get_stnode(index - x);
                    }
                    else
                    {
                        x += k;
                    }
                }
                std::cout << index << "/" << this->node_count << "/" << x << std::endl;
                assert(false);
                throw -1;
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
                /*
                if (this->child_count > 10000)
                {
                    std::cout << "LCP = " << this->current_lcp << ", node count = " << this->node_count << ", child count = " << this->child_count << std::endl;
                }
                */
                if (current_lcp > 0)
                {
                    bool b = this->child_count < 1000 || this->sub_trees.size() == 1;
                    //bool b = true;
                    if (b)
                    {
                        this->single_process();
                    }
                    else
                    {
                        this->allocate_data();
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
                for(uint64_t i=0;i<this->sub_trees.size();i++){
                    this->sub_trees[i].print2(*this->_RLBWTDS);
                }
                */
                
                uint64_t current_child_count = 0;
                uint64_t current_node_count = 0;

                this->node_count = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    current_child_count += sub_trees[i].children_count();
                    current_node_count += sub_trees[i].node_count();
                }
                total_counter += current_child_count;
                this->child_count = current_child_count;
                this->node_count = current_node_count;
                if (current_child_count > this->peak_child_count)
                {
                    this->peak_child_count = current_child_count;
                }
                assert(total_counter <= strSize);

                current_lcp++;
                assert(this->child_count > 0);
            }
            bool isStop()
            {
                return total_counter == strSize;
            }

            void print()
            {
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    std::cout << "[" << i << ": " << this->sub_trees[i].node_count() << ", " << this->sub_trees[i].children_count() << "]" << std::flush;
                }
                std::cout << std::endl;
            }

            void allocate_data()
            {
                uint64_t child_avg_count = (this->child_count / this->get_tree_count()) + 256;
                std::queue<uint64_t> lower_indexes, upper_indexes;
                for (uint64_t i = 0; i < this->get_tree_count(); i++)
                {
                    if (this->sub_trees[i].children_count() <= child_avg_count)
                    {
                        lower_indexes.push(i);
                    }
                    else
                    {
                        upper_indexes.push(i);
                    }
                }
                if (upper_indexes.size() == 0)
                    return;
                /*
                std::cout << "ALLOCATE, "
                          << "AVERAGE = " << child_avg_count << std::endl;
                this->print();
                */

                while (upper_indexes.size() > 0)
                {
                    uint64_t i = upper_indexes.front();

                    while (this->sub_trees[i].children_count() > child_avg_count)
                    {
                        //this->print();

                        //std::cout << "-" << std::endl;
                        //this->sub_trees[i].print();
                        assert(lower_indexes.size() > 0);
                        uint64_t j = lower_indexes.front();
                        //std::cout << i << " -> " << j << std::endl;

                        this->sub_trees[i].spill(this->sub_trees[j], child_avg_count);
                        if (this->sub_trees[j].children_count() > child_avg_count)
                        {
                            lower_indexes.pop();
                        }
                        //this->print();
                    }
                    upper_indexes.pop();
                }
                //this->print();
            }
            /*
            template <typename FUNC, typename OUTPUT>
            void execute_function_with_multi_thread(FUNC func, int thread_num, std::vector<std::vector<OUTPUT>> &output){
                uint64_t psize = this->node_count / thread_num;
                std::vector<thread> threads;
                uint64_t x = 0;
                while(true){
                    threads.push_back(thread(parallel_process_stnodes<INDEX_SIZE, RLBWTDS>, ref(this), ref(sub_tmp_trees[i]), ref(ems[i])));

                }
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    threads.push_back(thread(parallel_process_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees[i]), ref(sub_tmp_trees[i]), ref(ems[i])));
                }
                for (thread &t : threads)
                    t.join();

            }
            */
        private:
            void single_process()
            {
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    parallel_process_stnodes<INDEX_SIZE, RLBWTDS>(sub_trees[i], sub_tmp_trees[i], ems[i]);
                }
            }
            void parallel_process()
            {
                std::vector<thread> threads;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    threads.push_back(thread(parallel_process_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees[i]), ref(sub_tmp_trees[i]), ref(ems[i])));
                }
                for (thread &t : threads)
                    t.join();
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool