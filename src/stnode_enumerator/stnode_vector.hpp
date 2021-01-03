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
            uint64_t mini_increment(uint64_t L, uint64_t &left, uint64_t &right, bool &leftmost) const
            {
                assert(L < this->first_child_flag_vec.size());
                if (this->first_child_flag_vec[L])
                {
                    return mini_increment(L+1, left, right, leftmost);
                }
                else
                {
                    if (this->first_child_flag_vec[L - 1])
                    {
                        left = this->childs_vec[L - 1];
                        right = this->childs_vec[L];
                        leftmost = true;
                    }
                    else
                    {
                        left = this->childs_vec[L - 1] + 1;
                        right = this->childs_vec[L];
                        leftmost = false;

                    }
                    return L+1;

                }
            }
            void get_last(uint64_t lcp, std::vector<stool::LCPInterval<INDEX_SIZE>> &output) const {
                uint64_t width = this->get_last_width() - 1;
                uint64_t left = 0; 
                uint64_t right = 0;
                bool b = false;
                uint64_t L = this->first_child_flag_vec.size() - width;
                for(uint64_t i = 0;i<width;i++){
                    L = this->mini_increment(L, left, right, b);
                    output.push_back(stool::LCPInterval<INDEX_SIZE>(left, right, lcp));
                }
            }
            void pop(){
                uint64_t width = this->get_last_width();
                for(uint64_t i = 0;i<width;i++){
                    this->first_child_flag_vec.pop_back();
                    this->childs_vec.pop_back();
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
                    if (j == 0)
                    {
                        this->childs_vec.push_back(left);
                        this->first_child_flag_vec.push_back(true);
                    }
                    this->childs_vec.push_back(right);
                    this->first_child_flag_vec.push_back(false);
                }
                bool isMaximalRepeat = ds->checkMaximalRepeat(st_left, st_right);

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
            
            void print() const
            {
                std::cout << "[STNODE_COUNT, CHILDREN_COUNT] = [" << this->maximal_repeat_check_vec.size() << ", " << "XXX" << "]" << std::endl;
                stool::Printer::print("child_vec", this->childs_vec);
                stool::Printer::print_bits("first_child_flag_vec", this->first_child_flag_vec);
                stool::Printer::print_bits("maximal_repeat_check_vec", this->maximal_repeat_check_vec);
            }
        };
    } // namespace lcp_on_rlbwt
} // namespace stool