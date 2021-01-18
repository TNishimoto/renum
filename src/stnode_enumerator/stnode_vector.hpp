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
#include "./node_iterator.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace stnode_on_rlbwt
    {
        /*
        int64_t debug_sum_counter = 0;
        int64_t debug_peak_counter = 0;
        */
        //std::mutex mtx2;

        const uint ARRAYSIZE = 2000;
        const uint CARRAYSIZE = 4000;

        template <typename INDEX_SIZE, typename CHAR = uint8_t>
        class STNodeVector
        {

        public:
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using ITERATOR = STNodeIterator<STNodeVector<INDEX_SIZE, CHAR>>;
            using index_type = INDEX_SIZE;

            std::vector<INDEX_SIZE> childs_vec;
            std::vector<CHAR> edge_char_vec;
            std::vector<bool> first_child_flag_vec;
            std::vector<bool> maximal_repeat_check_vec;
            std::vector<INDEX_SIZE> depth_vec;

            STNodeVector()
            {
            }
            void clear()
            {
                this->childs_vec.clear();
                this->edge_char_vec.clear();
                this->first_child_flag_vec.clear();
                this->maximal_repeat_check_vec.clear();
            }
            uint64_t size() const
            {
                return this->maximal_repeat_check_vec.size();
            }
            uint64_t get_last_width() const
            {
                uint64_t p = UINT64_MAX;
                for (int64_t i = this->first_child_flag_vec.size() - 1; i >= 0; i--)
                {
                    if (this->first_child_flag_vec[i])
                    {
                        p = i;
                        break;
                    }
                }
                return this->first_child_flag_vec.size() - p;
            }
            uint64_t mini_increment(uint64_t L, stool::CharInterval<INDEX_SIZE, uint8_t> &output, bool &leftmost) const
            {
                assert(L < this->first_child_flag_vec.size());
                if (this->first_child_flag_vec[L])
                {
                    return mini_increment(L + 1, output, leftmost);
                }
                else
                {
                    if (this->first_child_flag_vec[L - 1])
                    {
                        output.i = this->childs_vec[L - 1];
                        output.j = this->childs_vec[L];
                        output.c = this->edge_char_vec.size() > 0 ? this->edge_char_vec[L] : 0;
                        leftmost = true;
                    }
                    else
                    {
                        output.i = this->childs_vec[L - 1] + 1;
                        output.j = this->childs_vec[L];
                        output.c = this->edge_char_vec.size() > 0 ? this->edge_char_vec[L] : 0;
                        leftmost = false;
                    }
                    return L + 1;
                }
            }
            void get_last(std::vector<stool::CharInterval<INDEX_SIZE, uint8_t>> &output) const
            {
                uint64_t width = this->get_last_width() - 1;
                stool::CharInterval<INDEX_SIZE, uint8_t> tmp;
                bool b = false;
                uint64_t L = this->first_child_flag_vec.size() - width;
                for (uint64_t i = 0; i < width; i++)
                {
                    L = this->mini_increment(L, tmp, b);
                    output.push_back(tmp);
                }
            }
            uint64_t get_last_left_boundary() const
            {
                assert(this->childs_vec.size() > 0);
                uint64_t i = this->childs_vec.size() - this->get_last_width();
                return this->childs_vec[i];
            }
            uint64_t get_last_right_boundary() const
            {
                assert(this->childs_vec.size() > 0);
                return this->childs_vec[this->childs_vec.size()-1];
            }

            int64_t get_last_depth() const {
                assert(this->depth_vec.size() != 0);
                return this->depth_vec[this->depth_vec.size()-1];
            }
            /*
            std::pair<INDEX_SIZE,INDEX_SIZE> get_last_node() const
            {
                uint64_t width = this->get_last_width() - 1;
                stool::CharInterval<INDEX_SIZE, uint8_t> tmp;
                bool b = false;
                uint64_t L = this->first_child_flag_vec.size() - width;
                for (uint64_t i = 0; i < width; i++)
                {
                    L = this->mini_increment(L, tmp, b);
                    output.push_back(tmp);
                }
            }
            */

            void pop()
            {
                uint64_t width = this->get_last_width();
                for (uint64_t i = 0; i < width; i++)
                {
                    this->first_child_flag_vec.pop_back();
                    this->childs_vec.pop_back();
                    if (this->edge_char_vec.size() > 0)
                    {
                        this->edge_char_vec.pop_back();
                    }
                }
                this->maximal_repeat_check_vec.pop_back();
                if(this->depth_vec.size() > 0){
                    this->depth_vec.pop_back();
                }
            }

            std::pair<uint64_t, uint64_t> compute_import_positions(uint64_t capacity)
            {
                uint64_t node_count = 0;
                uint64_t child_count = 0;

                int64_t i = this->first_child_flag_vec.size() - 1;
                uint64_t width = 0;
                while (capacity > 0 && i >= 0)
                {
                    width++;
                    if (this->first_child_flag_vec[i])
                    {
                        if (capacity >= width)
                        {
                            node_count++;
                            child_count += width;
                            capacity -= width;
                            width = 0;
                        }
                        else
                        {
                            break;
                        }
                    }
                    i--;
                }
                return std::pair<uint64_t, uint64_t>(node_count, child_count - node_count);
            }

            void print() const
            {
                std::cout << "[STNODE_COUNT, CHILDREN_COUNT] = [" << this->maximal_repeat_check_vec.size() << ", "
                          << "XXX"
                          << "]" << std::endl;
                stool::Printer::print("child_vec", this->childs_vec);
                stool::Printer::print_bits("first_child_flag_vec", this->first_child_flag_vec);
                stool::Printer::print_bits("maximal_repeat_check_vec", this->maximal_repeat_check_vec);
            }
            uint64_t increment(uint64_t L) const
            {
                assert(L > 0);
                assert(this->first_child_flag_vec[L - 1]);
                assert(!this->first_child_flag_vec[L]);

                uint64_t R = L + 1;
                while (R < this->first_child_flag_vec.size() && !this->first_child_flag_vec[R])
                {
                    R++;
                }
                return R + 1;
            }
            
            uint64_t increment(uint64_t L, uint64_t &left, uint64_t &right) const
            {
                assert(L > 0);
                assert(this->first_child_flag_vec[L - 1]);
                assert(!this->first_child_flag_vec[L]);

                uint64_t R = L + 1;
                while (R < this->first_child_flag_vec.size() && !this->first_child_flag_vec[R])
                {
                    R++;
                }
                left = this->get_child_left_boundary(L);
                right = this->get_child_right_boundary(R - 1);

                return R + 1;
            }
            void increment(ITERATOR &iter) const
            {
                this->increment2(iter.child_index, iter.node_index, iter.array_index);
            }
            void increment2(INDEX_SIZE &child_index, INDEX_SIZE &node_index, INDEX_SIZE &array_index) const
            {
                child_index = this->increment(child_index);
                if (child_index >= this->first_child_flag_vec.size())
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
            INDEX_SIZE get_children_count(const ITERATOR &iter) const
            {
                return this->get_children_count(iter.child_index);
            }
            INDEX_SIZE get_children_count(INDEX_SIZE child_index) const
            {
                uint64_t R = this->increment(child_index);
                return R - child_index - 1;
            }
            inline uint64_t get_child_left_boundary(uint64_t child_end_pointer) const
            {
                assert(child_end_pointer > 0);
                assert(child_end_pointer < this->childs_vec.size());
                if (this->first_child_flag_vec[child_end_pointer - 1])
                {
                    return this->childs_vec[child_end_pointer - 1];
                }
                else
                {
                    return this->childs_vec[child_end_pointer - 1] + 1;
                }
            }
            inline uint64_t get_child_right_boundary(uint64_t child_end_pointer) const
            {
                assert(child_end_pointer < this->childs_vec.size());
                return this->childs_vec[child_end_pointer];
            }
            CHAR get_edge_character(uint64_t child_pointer) const
            {
                return this->edge_char_vec[child_pointer];
            }
            INDEX_SIZE get_edge_character(const ITERATOR &iter, uint64_t ith_child) const
            {
                return this->get_edge_character(iter.child_index + ith_child);
            }

            INDEX_SIZE get_child_left_boundary(const ITERATOR &iter, uint64_t ith_child) const
            {
                return this->get_child_left_boundary(iter.child_index + ith_child);
            }
            INDEX_SIZE get_child_right_boundary(const ITERATOR &iter, uint64_t ith_child) const
            {
                return this->get_child_right_boundary(iter.child_index + ith_child);
            }

            INDEX_SIZE get_left(const ITERATOR &iter) const
            {
                return this->get_child_left_boundary(iter.child_index);
            }
            INDEX_SIZE get_left(INDEX_SIZE child_index) const
            {
                return this->get_child_left_boundary(child_index);
            }

            INDEX_SIZE get_right(const ITERATOR &iter) const
            {
                uint64_t left = 0, right = 0;
                this->increment(iter.child_index, left, right);
                return right;
            }
            INDEX_SIZE get_right(INDEX_SIZE child_index) const
            {
                uint64_t left = 0, right = 0;
                this->increment(child_index, left, right);
                return right;
            }

            STNodeIterator<STNodeVector> begin() const
            {
                return STNodeIterator<STNodeVector>(this, true);
            }
            STNodeIterator<STNodeVector> end() const
            {
                return STNodeIterator<STNodeVector>(this, false);
            }
            
            void set_current_first_iterator(ITERATOR &it) const
            {
                this->set_current_first_iterator(it.child_index, it.node_index, it.array_index);
            }
            void set_current_first_iterator(INDEX_SIZE &child_index, INDEX_SIZE &node_index, INDEX_SIZE &array_index) const
            {
                if (this->maximal_repeat_check_vec.size() == 0)
                {
                    child_index = std::numeric_limits<INDEX_SIZE>::max();
                    node_index = std::numeric_limits<INDEX_SIZE>::max();
                    array_index = std::numeric_limits<INDEX_SIZE>::max();
                }
                else
                {
                    child_index = 1;
                    node_index = 0;
                    array_index = std::numeric_limits<INDEX_SIZE>::max();
                }
            }
        };
    } // namespace stnode_on_rlbwt
} // namespace stool