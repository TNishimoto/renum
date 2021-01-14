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
    namespace stnode_on_rlbwt
    {

        const uint SUCCINCT_MODE = 1;
        const uint STANDARD_MODE = 0;
        const uint FAST_MODE = 2;
        const uint SINGLE_MODE = 3;

        template <typename INDEX_SIZE, typename RLBWTDS>
        class SuffixTreeNodes
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            //using STNODE_TRAVERSER = STNodeSubTraverser<INDEX_SIZE, RLBWTDS>;
            using SUCCINCT_STNODE_TRAVERSER = SuccinctSortedSTChildren<INDEX_SIZE, RLBWTDS>;
            using ITERATOR = STNodeIterator<SuffixTreeNodes>;
            using DEPTH_ITERATOR = STDepthIterator<SuffixTreeNodes>;

            STNodeTraverser<INDEX_SIZE, RLBWTDS> standard_st_traverser;
            FastSTNodeTraverser<INDEX_SIZE, RLBWTDS> fast_st_traverser;
            SingleSTNodeTraverser<INDEX_SIZE, ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS>> single_st_traverser;
            std::vector<ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS>> ems;
            bool store_edge_chars = false;

            //std::vector<SUCCINCT_STNODE_TRAVERSER *> succinct_sub_trees;

        public:
            using index_type = INDEX_SIZE;
            using CHAR = typename RLBWTDS::CHAR;
            uint64_t peak_child_count = 0;
            uint64_t thread_count = 1;
            bool use_fast_mode = true;
            uint mode = SINGLE_MODE;
            RLBWTDS *_RLBWTDS;
            stool::stnode_on_rlbwt::RLE<CHAR> *_RLBWT;

            ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS> *get_interval_search_deta_structure() const
            {
                ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS>* r = const_cast<ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS>*>(&ems[0]);
                return r;
            }


            bool has_edge_characters() const
            {
                return this->store_edge_chars;
            }
            /*
            uint64_t print_interval = 100;
            uint64_t print_interval_counter = 0;
            uint64_t total_counter = 0;
            uint64_t alpha = 2;
            uint64_t _expected_peak_memory_bits = 0;
            uint64_t _switch_threshold = 0;
            uint64_t debug_peak_memory = 0;
            */
            //uint mode = FAST_MODE;

            /*
            iterator begin() const
            {
                return iterator(this, true);
            }
            iterator end() const
            {
                return iterator(this, false);
            }
            */

            void clear()
            {
                this->standard_st_traverser.clear();
                this->fast_st_traverser.clear();
                this->single_st_traverser.clear();
                this->peak_child_count = 0;
                //this->total_counter = 0;
                this->mode = SINGLE_MODE;
            }
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

            void initialize(uint64_t _thread_count, RLBWTDS &_RLBWTDS, bool _store_edge_chars)
            {
                this->_RLBWTDS = &_RLBWTDS;
                _RLBWT = this->_RLBWTDS->get_rlbwt();
                this->store_edge_chars = _store_edge_chars;

                for (uint64_t i = 0; i < this->thread_count; i++)
                {
                    ems.push_back(ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS>());
                    ems[ems.size() - 1].initialize(this->_RLBWTDS);
                }

                //this->strSize = _RLBWTDS.str_size();

                this->standard_st_traverser.initialize(_thread_count, _RLBWTDS, this->store_edge_chars);

                this->fast_st_traverser.initialize(_thread_count, _RLBWTDS, this->store_edge_chars);

                this->single_st_traverser.initialize(&ems[0], this->store_edge_chars);
                this->thread_count = _thread_count;

                /*
                double ratio = (double)this->_RLBWTDS->str_size() / (double)this->_RLBWTDS->rle_size();
                double d = std::log2(ratio);
                this->_expected_peak_memory_bits = this->_RLBWTDS->rle_size() * d;
                this->_switch_threshold = this->alpha * (this->expected_peak_memory_bits() / (sizeof(uint64_t) * 8));
                */
            }
            int64_t get_current_lcp() const
            {
                return this->current_lcp();
            }

            int64_t current_lcp() const
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
                    this->fast_st_traverser.get_lcp_intervals(output);
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

            uint64_t get_input_text_length()
            {
                return this->_RLBWT->str_size();
            }

            bool succ()
            {
                if (this->is_finished())
                    return false;

#if DEBUG

                if (this->_RLBWT->str_size() < 100)
                {
                    std::cout << "Start Enumerate" << std::endl;
                    //this->print();
                }
#endif

                if (this->mode == STANDARD_MODE)
                {
                    this->standard_st_traverser.succ();
                }
                else if (this->mode == FAST_MODE)
                {
                    this->fast_st_traverser.succ();
                }
                else if (this->mode == SINGLE_MODE)
                {
                    this->single_st_traverser.succ();
                }
                else
                {
                    assert(false);
                    throw -1;
                }

                this->update_info();
                if (this->mode == SINGLE_MODE && this->thread_count > 1)
                {
                    std::cout << "Switch[SINGLE -> STANDARD]" << std::endl;
                    this->mode = STANDARD_MODE;
                    STNodeVector<INDEX_SIZE> tmp;
                    this->single_st_traverser.convert_to_vector(tmp);
                    this->standard_st_traverser.import(this->single_st_traverser.get_current_lcp(), tmp);
                }

                if (this->use_fast_mode && this->mode == STANDARD_MODE && (this->child_count() * 50 < this->peak_child_count))
                {
                    std::cout << "Switch[STANDARD -> FAST]" << std::endl;
                    STNodeVector<INDEX_SIZE> tmp;
                    uint64_t lcp = this->standard_st_traverser.get_current_lcp();

                    this->standard_st_traverser.to_stnode_vector(tmp);
                    this->standard_st_traverser.clear();

                    this->fast_st_traverser.set_peak(this->peak_child_count / 50);
                    this->fast_st_traverser.import(tmp, lcp);

                    this->mode = FAST_MODE;
                }

#if DEBUG

                if (this->_RLBWT->str_size() < 100)
                {
                    std::cout << "Enumerate END" << std::endl;
                    this->print();
                }

#endif
                return true;
            }
            void update_info()
            {
                //this->total_counter += this->child_count() - this->node_count();
                if (this->peak_child_count < this->child_count())
                {
                    this->peak_child_count = this->child_count();
                }
            }

            bool is_finished() const
            {
                bool b = false;
                if (this->mode == STANDARD_MODE)
                {
                    b = this->standard_st_traverser.is_finished();
                }
                else if (this->mode == FAST_MODE)
                {
                    b = this->fast_st_traverser.is_finished();
                }
                else if (this->mode == SINGLE_MODE)
                {
                    b = this->single_st_traverser.is_finished();
                }
                else
                {
                    assert(false);
                    throw -1;
                }
                return b;

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
                }
                else if (this->mode == SINGLE_MODE)
                {
                    this->single_st_traverser.print();
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
            /*
            void print_info()
            {
                std::cout << "[" << (this->print_interval_counter) << "/" << this->print_interval << "] ";
                std::cout << "LCP = " << this->current_lcp();
                std::cout << ", Peak = " << this->peak_child_count;
                std::cout << ", Current = " << this->child_count();
                std::cout << ", current_total = " << (this->_RLBWT->str_size() - this->total_counter);
                std::cout << ", Peak = " << this->debug_peak_memory / 1000 << "[KB]" << std::endl;
                //std::cout << "Peak = " << debug_peak_counter << "[KB]" << std::endl;
            }
            */

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

            void increment(ITERATOR &iter) const
            {

                if (this->mode == STANDARD_MODE)
                {
                    this->standard_st_traverser.increment(iter.child_index, iter.node_index, iter.array_index);
                }
                else if (this->mode == FAST_MODE)
                {
                    this->fast_st_traverser.increment(iter.child_index, iter.node_index, iter.array_index);
                }
                else if (this->mode == SINGLE_MODE)
                {
                    this->single_st_traverser.increment2(iter.child_index, iter.node_index, iter.array_index);
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
            INDEX_SIZE get_children_count(const ITERATOR &iter) const
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.get_children_count(iter.child_index, iter.array_index);
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.get_children_count(iter.child_index);
                }
                else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.get_children_count(iter.child_index);
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
            
            CHAR get_edge_character(const ITERATOR &iter, uint64_t ith_child) const
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.get_edge_character(iter.child_index + ith_child, iter.array_index);
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.get_edge_character(iter.child_index + ith_child);
                }
                else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.get_edge_character(iter.child_index + ith_child);
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
            INDEX_SIZE get_child_left_boundary(const ITERATOR &iter, uint64_t ith_child) const
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.get_child_left_boundary(iter.child_index + ith_child, iter.array_index);
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.get_child_left_boundary(iter.child_index + ith_child);
                }
                else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.get_child_left_boundary(iter.child_index + ith_child);
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
            INDEX_SIZE get_child_right_boundary(const ITERATOR &iter, uint64_t ith_child) const
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.get_child_right_boundary(iter.child_index + ith_child, iter.array_index);
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.get_child_right_boundary(iter.child_index + ith_child);
                }
                else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.get_child_right_boundary(iter.child_index + ith_child);
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }

            INDEX_SIZE get_left(const ITERATOR &iter) const
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.get_left(iter.child_index, iter.array_index);
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.get_left(iter.child_index);
                }
                else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.get_left(iter.child_index);
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
            INDEX_SIZE get_right(const ITERATOR &iter) const
            {

                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.get_right(iter.child_index, iter.array_index);
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.get_right(iter.child_index);
                }
                else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.get_right(iter.child_index);
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }
            void set_current_first_iterator(ITERATOR &iter) const
            {
                if (this->mode == STANDARD_MODE)
                {
                    this->standard_st_traverser.set_current_first_iterator(iter.child_index, iter.node_index, iter.array_index);
                }
                else if (this->mode == FAST_MODE)
                {
                    this->fast_st_traverser.set_current_first_iterator(iter.child_index, iter.node_index, iter.array_index);
                }
                else if (this->mode == SINGLE_MODE)
                {
                    this->single_st_traverser.set_current_first_iterator(iter.child_index, iter.node_index, iter.array_index);
                }
                else
                {
                    assert(false);
                    throw -1;
                }
            }

            bool check_maximal_repeat(const ITERATOR &iter) const
            {
                if (this->mode == STANDARD_MODE)
                {
                    return this->standard_st_traverser.check_maximal_repeat(iter.node_index, iter.array_index);
                }
                else if (this->mode == FAST_MODE)
                {
                    return this->fast_st_traverser.check_maximal_repeat(iter.node_index);
                }
                else if (this->mode == SINGLE_MODE)
                {
                    return this->single_st_traverser.check_maximal_repeat(iter.node_index);
                }
                else
                {
                    assert(false);
                    throw -1;
                }
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
                }
                else if (this->mode == SINGLE_MODE)
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

    } // namespace stnode_on_rlbwt
} // namespace stool