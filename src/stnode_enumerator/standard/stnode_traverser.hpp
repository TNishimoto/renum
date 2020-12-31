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
#include "stnode_sub_traverser.hpp"

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
        void parallel_process_stnodes(std::vector<STNodeSubTraverser<INDEX_SIZE, RLBWTDS> *> &trees, uint64_t fst_position,
                                      std::stack<uint64_t> &position_stack, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
        {

            uint64_t pos = fst_position;
            while (pos != UINT64_MAX)
            {

                trees[pos]->computeNextLCPIntervalSetForParallelProcessing(em);
                //assert(trees[pos]->children_count() <= limit);
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
        void parallel_count_stnodes(std::vector<STNodeSubTraverser<INDEX_SIZE, RLBWTDS> *> &trees, uint64_t fst_position,
                                    std::stack<uint64_t> &position_stack, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em, uint64_t &output_children_count)
        {
            output_children_count = 0;
            uint64_t pos = fst_position;
            while (pos != UINT64_MAX)
            {
                auto pair = trees[pos]->countNextLCPIntervalSet(em);
                output_children_count += pair.second;
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
        /*
        template <typename INDEX_SIZE, typename RLBWTDS>
        void parallel_test_stnodes(std::vector<STNodeSubTraverser<INDEX_SIZE, RLBWTDS> *> &trees, uint64_t fst_position,
                                    std::stack<uint64_t> &position_stack, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em, uint64_t &output_children_count)
        {
            output_children_count = 0;
            uint64_t pos = fst_position;
            while (pos != UINT64_MAX)
            {
                auto pair = trees[pos]->test(em);
                output_children_count += pair.second;
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
        */


        template <typename INDEX_SIZE, typename RLBWTDS>
        class STNodeTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using STNODE_SUB_TRAVERSER = STNodeSubTraverser<INDEX_SIZE, RLBWTDS>;

            std::vector<STNODE_SUB_TRAVERSER *> sub_trees;
            //std::vector<std::vector<STNODE_SUB_TRAVERSER *>> new_trees;

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
            void clear()
            {

                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    this->sub_trees[i]->clear();
                    delete this->sub_trees[i];
                    this->sub_trees[i] = nullptr;
                }
                this->sub_trees.clear();
                this->current_lcp = 0;
                this->_child_count = 0;
                this->_node_count = 0;
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

                auto st = new STNODE_SUB_TRAVERSER(this->_RLBWTDS);
                sub_trees.push_back(st);
                //sub_trees.resize(256);

                //this->new_trees.resize(size);

                for (uint64_t i = 0; i < this->thread_count; i++)
                {
                    /*
                    lightRDs.push_back(LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE>());
                    lightRDs[lightRDs.size() - 1].preprocess(&__RLBWTDS.bwt);

                    heavyRDs.push_back(SuccinctRangeDistinctDataStructure<INDEX_SIZE>());
                    heavyRDs[heavyRDs.size() - 1].initialize(&__RLBWTDS.wt, lastChar);
                    */
                    ems.push_back(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>());
                    ems[ems.size() - 1].initialize(&__RLBWTDS);
                }
                /*
                for (uint64_t i = 0; i < this->thread_count; i++)
                {

                    ems[i].lightDS = &lightRDs[i];
                    ems[i].heavyDS = &heavyRDs[i];
                }
                */
                /*
                double ratio = (double)this->_RLBWTDS->str_size() / (double)this->_RLBWTDS->rle_size();
                double d = std::log2(ratio);
                this->_expected_peak_memory_bits = this->_RLBWTDS->rle_size() * d;
                this->_switch_threshold = this->alpha * (this->expected_peak_memory_bits() / (sizeof(uint64_t) * 8));
                */
            }

            void load(std::ifstream &file)
            {
                this->clear();

                std::cout << "\033[33m";

                std::cout << "LOAD" << std::endl;
                file.seekg(0, std::ios::end);
                uint64_t textSize = (uint64_t)file.tellg() / sizeof(char);
                file.seekg(0, std::ios::beg);
                std::cout << "Data_SIZE " << textSize << std::endl;

                uint64_t sub_tree_count = 0;

                file.read((char *)&this->current_lcp, sizeof(uint64_t));
                file.read((char *)&sub_tree_count, sizeof(uint64_t));

                std::cout << "lcp = " << this->current_lcp << std::endl;
                std::cout << "sub_tree_count = " << sub_tree_count << std::endl;

                for (uint64_t i = 0; i < sub_tree_count; i++)
                {

                    auto st = new STNODE_SUB_TRAVERSER(this->_RLBWTDS);
                    st->load(file);
                    this->sub_trees.push_back(st);
                    this->_child_count += st->children_count();
                    this->_node_count += st->node_count();
                }
                //this->print();
                std::cout << "\033[39m" << std::endl;
            }
            void write(std::ofstream &out)
            {
                uint64_t sub_tree_size = this->sub_trees.size();
                uint64_t lcp = this->current_lcp;
                std::cout << "\033[32m";
                std::cout << "WRITE " << std::endl;

                std::cout << "Subtree count : " << sub_tree_size << std::endl;
                std::cout << "LCP : " << lcp << std::endl;

                out.write((char *)&lcp, sizeof(uint64_t));
                out.write((char *)&sub_tree_size, sizeof(uint64_t));

                for (auto &it : this->sub_trees)
                {
                    it->write(out);
                }
                out.close();
                //this->print();
                std::cout << "\033[39m" << std::endl;
            }

            uint64_t parallel_count()
            {

                std::vector<uint64_t> fst_pos_vec;
                std::vector<uint64_t> children_count_vec;

                fst_pos_vec.resize(this->thread_count, UINT64_MAX);
                children_count_vec.resize(this->thread_count, 0);

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
                    threads.push_back(thread(parallel_count_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees), fst_pos_vec[i], ref(position_stack), ref(ems[i]), ref(children_count_vec[i])));

                }

                for (thread &t : threads)
                    t.join();

                assert(position_stack.size() == 0);
                uint64_t count = 0;
                for (auto &it : children_count_vec)
                {
                    count += it;
                }
                return count;
            }
            uint64_t write_maximal_repeats(std::ofstream &out)
            {
                uint64_t count = 0;
                for (auto &it : this->sub_trees)
                {
                    count += it->write_maximal_repeats(this->current_lcp - 1, out);
                    /*
                    it->get_lcp_intervals(this->current_lcp - 1, tmp);

                    for (auto &intv : tmp)
                    {
                        stool::LCPInterval<INDEX_SIZE> newLCPIntv(intv.i, intv.j, intv.lcp);
                        buffer.push_back(newLCPIntv);
                        tmp.clear();
                        if (buffer.size() >= 8000)
                        {
                            out.write(reinterpret_cast<const char *>(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                            buffer.clear();
                        }
                    }
                    */
                    /*
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
                    */
                }
                return count;
            }
            void get_lcp_intervals(std::vector<stool::LCPInterval<uint64_t>> &output)
            {
                for (auto &it : this->sub_trees)
                {
                    it->get_lcp_intervals(this->current_lcp - 1, output);
                }
            }
            /*
            uint64_t get_stnode(uint64_t L, stool::LCPInterval<uint64_t> &output)
            {
                uint64_t p = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    uint64_t csize = this->sub_trees[i]->children_count();

                    if (p <= L && L < p + csize)
                    {
                        return p + this->sub_trees[i]->get_stnode(L - p, output, this->current_lcp - 1);
                    }
                    else
                    {
                        p += csize;
                    }
                }
                return UINT64_MAX;
            }
            */

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
                std::cout << "PRINT STNODE Traverser" << std::endl;
                std::cout << "[LCP, STNODE_COUNT, CHILDREN_COUNT] = [" << this->current_lcp << ", " << this->_node_count << ", " << this->_child_count << "]" << std::endl;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    this->sub_trees[i]->print();
                }

                std::cout << "[END]" << std::endl;
            }

            uint64_t get_using_memory()
            {
                uint64_t k = 0;

                uint64_t p = this->sub_trees.size() * 32;

                for (auto &it : this->sub_trees)
                {
                    k += it->get_using_memory();
                }

                return p + k;
            }

            void move(std::deque<INDEX_SIZE> &item1, std::deque<bool> &item2, std::deque<bool> &item3)
            {
                for (auto &it : this->sub_trees)
                {
                    it->pop(item1, item2, item3);
                }
                this->remove_empty_trees();
            }

        private:
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
                    threads.push_back(thread(parallel_process_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees), fst_pos_vec[i], ref(position_stack), ref(ems[i])));
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
                /*
                if (this->child_count() == 0)
                    return;
                uint64_t fst_pos = 0;
                */
                for (auto &it : this->sub_trees)
                {
                    it->computeNextLCPIntervalSet(ems[0]);
                }
                /*
                for (uint64_t i = 1; i < this->sub_trees.size(); i++)
                {
                    position_stack.push(i);
                }

                assert(this->child_count() > 0);

                parallel_process_stnodes<INDEX_SIZE, RLBWTDS>(ref(sub_trees), fst_pos, ref(position_stack), ref(ems[0]));
                assert(position_stack.size() == 0);
                */
            }
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
                /*
                for (auto &it : this->new_trees)
                {
                    while (it.size() > 0)
                    {

                        this->sub_trees.push_back(it[it.size() - 1]);
                        it.pop_back();
                    }
                }
                */

                this->merge();
                this->remove_empty_trees();
                this->split();

                if ((double)this->sub_trees.size() * 2 < (double)this->sub_trees.capacity())
                {
                    this->sub_trees.shrink_to_fit();
                }
            }
            void split()
            {
                uint64_t x = 0;
                while (x < this->sub_trees.size())
                {
                    if (this->sub_trees[x]->children_count() > this->sub_tree_limit_size)
                    {
                        auto st = new STNODE_SUB_TRAVERSER(this->_RLBWTDS);
                        this->sub_trees[x]->split(*st);
                        this->sub_trees.push_back(st);
                    }
                    else
                    {
                        x++;
                    }
                }
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