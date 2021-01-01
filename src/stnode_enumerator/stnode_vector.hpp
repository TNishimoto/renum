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

        template <typename INDEX_SIZE>
        class STNodeVector
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;

        public:
            std::vector<INDEX_SIZE> childs_vec;
            std::vector<bool> first_child_flag_vec;
            std::vector<bool> maximal_repeat_check_vec;

            STNodeVector()
            {
            }
            void clear()
            {
                this->childs_vec.clear();
                this->first_child_flag_vec.clear();
                this->maximal_repeat_check_vec.clear();
            }
            uint64_t size()
            {
                return this->maximal_repeat_check_vec.size();
            }
            uint64_t get_last_width()
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
                return std::pair<uint64_t, uint64_t>(node_count, child_count);
            }
            //void move_push(std::vector<INDEX_SIZE> &_childs_vec, std::vector<bool> &_first_child_flag_vec, std::vector<bool> &_maximal_repeat_check_vec){

            /*
            void move(std::deque<INDEX_SIZE> &_childs_vec, std::deque<bool> &_first_child_flag_vec, std::deque<bool> &_maximal_repeat_check_vec)
            {
                uint64_t minSize = std::min(this->first_child_flag_vec.size(), _first_child_flag_vec.size());
                for (uint64_t i = 0; i < minSize; i++)
                {
                    _childs_vec[(i * 2)] = this->childs_vec[(i * 2)];
                    _childs_vec[(i * 2) + 1] = this->childs_vec[(i * 2) + 1];
                    _first_child_flag_vec[i] = this->first_child_flag_vec[i];
                }
                for (uint64_t i = minSize; i < this->first_child_flag_vec.size(); i++)
                {
                    _childs_vec.push_back(this->childs_vec[(i * 2)]);
                    _childs_vec.push_back(this->childs_vec[(i * 2) + 1]);
                    _first_child_flag_vec.push_back(this->first_child_flag_vec[i]);
                }
                while (_first_child_flag_vec.size() > this->first_child_flag_vec.size())
                {
                    _childs_vec.pop_back();
                    _childs_vec.pop_back();
                    _first_child_flag_vec.pop_back();
                }

                uint64_t minSTSize = std::min(this->maximal_repeat_check_vec.size(), _maximal_repeat_check_vec.size());
                for (uint64_t i = 0; i < minSTSize; i++)
                {
                    _maximal_repeat_check_vec[i] = this->maximal_repeat_check_vec[i];
                }
                for (uint64_t i = minSTSize; i < this->maximal_repeat_check_vec.size(); i++)
                {
                    _maximal_repeat_check_vec.push_back(this->maximal_repeat_check_vec[i]);
                }
                while (_maximal_repeat_check_vec.size() > this->maximal_repeat_check_vec.size())
                {
                    _maximal_repeat_check_vec.pop_back();
                }
                this->clear();
            }
            */

            template <typename RLBWTDS>
            void add(uint8_t c, uint64_t count, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em,
                     RLBWTDS *ds)
            {
                RINTERVAL copy;
                uint64_t st_left = UINT64_MAX;
                uint64_t st_right = 0;

                for (uint64_t j = 0; j < count; j++)
                {
                    em.get_child(c, j, copy);
                    uint64_t left = ds->get_fpos(copy.beginIndex, copy.beginDiff);
                    uint64_t right = ds->get_fpos(copy.endIndex, copy.endDiff);

                    if (left < st_left)
                    {
                        st_left = left;
                    }
                    if (right > st_right)
                    {
                        st_right = right;
                    }

                    this->childs_vec.push_back(left);
                    this->childs_vec.push_back(right);
                    this->first_child_flag_vec.push_back(j == 0);
                }
                uint64_t x = ds->get_lindex_containing_the_position(st_left);
                uint64_t d = ds->get_run(x);
                bool isMaximalRepeat = (ds->get_lpos(x) + d) <= st_right;
                this->maximal_repeat_check_vec.push_back(isMaximalRepeat);
            }
            template <typename RLBWTDS>
            void import(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em, RLBWTDS *ds)
            {
                for (uint64_t i = 0; i < em.indexCount; i++)
                {
                    auto c = em.indexVec[i];
                    auto &currentVec = em.childrenVec[c];
                    uint64_t count = currentVec.size();
                    this->add(c, count, em, ds);
                }
            }
        };
    } // namespace lcp_on_rlbwt
} // namespace stool