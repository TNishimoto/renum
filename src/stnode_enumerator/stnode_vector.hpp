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

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
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
            using RINTERVAL = RInterval<INDEX_SIZE>;

        public:
            std::vector<INDEX_SIZE> childs_vec;
            std::vector<CHAR> edge_char_vec;
            std::vector<bool> first_child_flag_vec;
            std::vector<bool> maximal_repeat_check_vec;

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
            uint64_t mini_increment(uint64_t L, stool::CharInterval<INDEX_SIZE> &output, bool &leftmost) const
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
            void get_last(std::vector<stool::CharInterval<INDEX_SIZE>> &output) const
            {
                uint64_t width = this->get_last_width() - 1;
                stool::CharInterval<INDEX_SIZE> tmp;
                bool b = false;
                uint64_t L = this->first_child_flag_vec.size() - width;
                for (uint64_t i = 0; i < width; i++)
                {
                    L = this->mini_increment(L, tmp, b);
                    output.push_back(tmp);
                }
            }
            void pop()
            {
                uint64_t width = this->get_last_width();
                for (uint64_t i = 0; i < width; i++)
                {
                    this->first_child_flag_vec.pop_back();
                    this->childs_vec.pop_back();
                    if(this->edge_char_vec.size() > 0){
                        this->edge_char_vec.pop_back();

                    }
                }
                this->maximal_repeat_check_vec.pop_back();
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
            /*
            template <typename RLBWTDS>
            void add(uint8_t c, uint64_t count, ExplicitWeinerLinkEmulator<RLBWTDS> &em, bool _store_edge_chars)
            {
                //RINTERVAL copy;
                uint64_t st_left = UINT64_MAX;
                uint64_t st_right = 0;
                std::pair<INDEX_SIZE, INDEX_SIZE> outputInterval;

                for (uint64_t j = 0; j < count; j++)
                {
                    CHAR edgeChar = em.get_child(c, j, outputInterval);
                    //uint64_t left = ds->get_fpos(copy.beginIndex, copy.beginDiff);
                    //uint64_t right = ds->get_fpos(copy.endIndex, copy.endDiff);

                    if (outputInterval.first < st_left)
                    {
                        st_left = outputInterval.first;
                    }
                    if (outputInterval.second > st_right)
                    {
                        st_right = outputInterval.second;
                    }
                    if (j == 0)
                    {
                        this->childs_vec.push_back(outputInterval.first);
                        this->first_child_flag_vec.push_back(true);
                        if (_store_edge_chars)
                        {
                            this->edge_char_vec.push_back(0);
                        }
                    }
                    this->childs_vec.push_back(outputInterval.second);
                    this->first_child_flag_vec.push_back(false);
                    if (_store_edge_chars)
                    {
                        this->edge_char_vec.push_back(edgeChar);
                    }
                }
                bool isMaximalRepeat = em.checkMaximalRepeat(st_left, st_right);

                this->maximal_repeat_check_vec.push_back(isMaximalRepeat);
            }
            */
           /*
            template <typename RLBWTDS>
            void import(ExplicitWeinerLinkEmulator<RLBWTDS> &em, uint8_t c, bool _store_edge_chars)
            {
                auto &currentVec = em.childrenVec[c];
                uint64_t count = currentVec.size();
                this->add(c, count, em, _store_edge_chars);
            }
            */

            void print() const
            {
                std::cout << "[STNODE_COUNT, CHILDREN_COUNT] = [" << this->maximal_repeat_check_vec.size() << ", "
                          << "XXX"
                          << "]" << std::endl;
                stool::Printer::print("child_vec", this->childs_vec);
                stool::Printer::print_bits("first_child_flag_vec", this->first_child_flag_vec);
                stool::Printer::print_bits("maximal_repeat_check_vec", this->maximal_repeat_check_vec);
            }
        };
    } // namespace lcp_on_rlbwt
} // namespace stool