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

#include "standard/stnode_traverser.hpp"
#include "fast/fast_stnode_traverser.hpp"
#include "single/single_stnode_traverser.hpp"

#include "succinct/succinct_sorted_stchildren_builder.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {

        const uint SUCCINCT_MODE = 1;
        const uint STANDARD_MODE = 0;
        const uint FAST_MODE = 2;
        const uint SINGLE_MODE = 3;

        template <typename INDEX_SIZE, typename RLBWTDS>
        class STNodeEnumerator
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            //using STNODE_TRAVERSER = STNodeSubTraverser<INDEX_SIZE, RLBWTDS>;
            using SUCCINCT_STNODE_TRAVERSER = SuccinctSortedSTChildren<INDEX_SIZE, RLBWTDS>;

            STNodeTraverser<INDEX_SIZE, RLBWTDS> standard_st_traverser;
            FastSTNodeTraverser<INDEX_SIZE, RLBWTDS> fast_st_traverser;
            SingleSTNodeTraverser<INDEX_SIZE, RLBWTDS> single_st_traverser;

            //std::vector<SUCCINCT_STNODE_TRAVERSER *> succinct_sub_trees;

            uint64_t print_interval = 100;
            uint64_t print_interval_counter = 0;

        public:
            uint64_t total_counter = 0;
            uint64_t peak_child_count = 0;
            uint64_t alpha = 2;
            uint64_t _expected_peak_memory_bits = 0;
            uint64_t _switch_threshold = 0;
            uint64_t debug_peak_memory = 0;
            uint64_t thread_count = 1;

            uint mode = SINGLE_MODE;
            //uint mode = FAST_MODE;

            RLBWTDS *_RLBWTDS;

            uint64_t node_count() const
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.node_count();
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.node_count();
                }
                else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.node_count();
                }
                else
                {
                    assert(false);
                    throw -1;
                }

                //return this->sub_tree.node_count();
            }
            uint64_t child_count() const
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.child_count();
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.child_count();
                }
                else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.child_count();
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }

            uint64_t expected_peak_memory_bits()
            {
                return this->_expected_peak_memory_bits;
            }

            uint64_t switch_threshold()
            {
                return this->_switch_threshold;
            }

            void initialize(uint64_t _thread_count, RLBWTDS &_RLBWTDS)
            {
                this->_RLBWTDS = &_RLBWTDS;
                //this->strSize = _RLBWTDS.str_size();
                this->standard_st_traverser.initialize(_thread_count, _RLBWTDS);
                this->fast_st_traverser.initialize(_thread_count, _RLBWTDS);
                this->single_st_traverser.initialize(&_RLBWTDS);
                this->thread_count = _thread_count;

#if DEBUG
                if (this->_RLBWTDS->str_size() < 100)
                {
                    this->print_interval = _RLBWTDS.str_size();
                }

#endif

                /*
                double ratio = (double)this->_RLBWTDS->str_size() / (double)this->_RLBWTDS->rle_size();
                double d = std::log2(ratio);
                this->_expected_peak_memory_bits = this->_RLBWTDS->rle_size() * d;
                this->_switch_threshold = this->alpha * (this->expected_peak_memory_bits() / (sizeof(uint64_t) * 8));
                */
            }
            uint64_t write_maximal_repeats(std::ofstream &out)
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.write_maximal_repeats(out);
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.write_maximal_repeats(out);
                }
                else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.write_maximal_repeats(out);
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
            int64_t current_lcp()
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.get_current_lcp();
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.get_current_lcp();
                }
                else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.get_current_lcp();
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }

            void get_lcp_intervals(std::vector<stool::LCPInterval<uint64_t>> &output)
            {
                if (this->mode == STANDARD_MODE)
                {
                    this->standard_st_traverser.get_lcp_intervals(output);
                }
                else if (this->mode == FAST_MODE)
                {
                    assert(false);
                    throw -1;
                }
                else if (this->mode == SINGLE_MODE)
                {
                    this->single_st_traverser.get_lcp_intervals(output);
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
            /*
            uint64_t get_stnode(uint64_t L, stool::LCPInterval<uint64_t> &output)
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.get_stnode(L, output);
                }
                else if (this->mode == FAST_MODE)
                {
                    RINTERVAL output1;
                    uint64_t x = this->fast_st_traverser.get_stnode(L, output1);

                    output.i = this->_RLBWTDS->get_lpos(output1.beginIndex) + output1.beginDiff;
                    output.j = this->_RLBWTDS->get_lpos(output1.endIndex) + output1.endDiff;
                    output.lcp = this->fast_st_traverser.get_current_lcp() - 1;

                    return x;
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
            */

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
                        this->print_info();

                        this->print_interval_counter++;
                    }
                }
#if DEBUG

                if (this->_RLBWTDS->str_size() < 100)
                {
                    std::cout << "Start Enumerate" << std::endl;
                    //this->print();
                }
#endif

                if (this->mode == STANDARD_MODE)
                {
                    this->standard_st_traverser.process();
                }
                else if (this->mode == FAST_MODE)
                {
                    this->fast_st_traverser.process();
                }else if (this->mode == SINGLE_MODE)
                {
                    this->single_st_traverser.process();
                }
                else
                {
                    assert(false);
                    throw -1;
                }

                this->update_info();
                if(this->mode == SINGLE_MODE && this->thread_count > 1){
                    std::cout << "Switch -> STANDARD" << std::endl;
                    this->mode = STANDARD_MODE;
                    STNodeVector<INDEX_SIZE> tmp;
                    this->single_st_traverser.convert_to_vector(tmp);
                    this->standard_st_traverser.import(this->single_st_traverser.get_current_lcp(),tmp);

                }
                /*
                if (this->mode == STANDARD_MODE && (this->child_count() * 10 < this->peak_child_count))
                {
                    std::cout << "Switch" << std::endl;

                    std::deque<INDEX_SIZE> childs_vec;
                    std::deque<bool> first_child_flag_vec;
                    std::deque<bool> maximal_repeat_check_vec;
                    this->standard_st_traverser.move(childs_vec, first_child_flag_vec, maximal_repeat_check_vec);
                    this->fast_st_traverser.import(childs_vec, first_child_flag_vec, maximal_repeat_check_vec, this->standard_st_traverser.get_current_lcp()-1 );
                    this->fast_st_traverser.set_peak(this->peak_child_count);
                    this->mode = FAST_MODE;
                    this->fast_st_traverser.process();

                }
                */

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
                if (this->peak_child_count < this->child_count())
                {
                    this->peak_child_count = this->child_count();
                }
            }

            bool isStop()
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.isStop();
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.isStop();
                }else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.get_current_lcp() >= 0 && this->single_st_traverser.node_count() == 0;
                }
                else
                {
                    assert(false);
                    throw -1;
                }
                //return total_counter == strSize - 1;
            }

            void print()
            {
                if (this->mode == STANDARD_MODE)
                {
                    this->standard_st_traverser.print();
                }
                else if (this->mode == FAST_MODE)
                {
                    this->fast_st_traverser.print();
                }else if (this->mode == SINGLE_MODE)
                {
                    this->single_st_traverser.print();
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
            void print_info()
            {
                std::cout << "[" << (this->print_interval_counter) << "/" << this->print_interval << "] ";
                std::cout << "LCP = " << this->current_lcp();
                std::cout << ", Peak = " << this->peak_child_count;
                std::cout << ", Current = " << this->child_count();
                std::cout << ", current_total = " << (this->_RLBWTDS->str_size() - this->total_counter);
                std::cout << ", Peak = " << this->debug_peak_memory / 1000 << "[KB]" << std::endl;
                //std::cout << "Peak = " << debug_peak_counter << "[KB]" << std::endl;
                /*
                        if(this->print_interval_counter == 3){
                            throw -1;
                        }
                        */
            }

        private:
            uint64_t get_using_memory()
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.get_using_memory();
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.get_using_memory();
                }else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.get_using_memory();
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool