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

        template <typename INDEX_SIZE>
        class STNodeVector
        {
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
            uint64_t size(){
                return this->maximal_repeat_check_vec.size();
            }

            //void move_push(std::vector<INDEX_SIZE> &_childs_vec, std::vector<bool> &_first_child_flag_vec, std::vector<bool> &_maximal_repeat_check_vec){
            void move_push(std::deque<INDEX_SIZE> &_childs_vec, std::deque<bool> &_first_child_flag_vec, std::deque<bool> &_maximal_repeat_check_vec){
                uint64_t p = UINT64_MAX;
                for(int64_t i = this->first_child_flag_vec.size() -1; i >= 0;i--){
                    if(this->first_child_flag_vec[i]){
                        p = i;
                        break;
                    }
                }
                uint64_t fsize = this->first_child_flag_vec.size();
                for(uint64_t i = p;i<fsize;i++){
                    _childs_vec.push_back(this->childs_vec[(i*2)]);
                    _childs_vec.push_back(this->childs_vec[(i*2)+1]);
                    _first_child_flag_vec.push_back(this->first_child_flag_vec[i]);
                }
                for(uint64_t i = p;i<fsize;i++){
                    this->childs_vec.pop_back();
                    this->childs_vec.pop_back();
                    this->first_child_flag_vec.pop_back();
                }
                _maximal_repeat_check_vec.push_back(this->maximal_repeat_check_vec[this->maximal_repeat_check_vec.size()-1]);
                this->maximal_repeat_check_vec.pop_back();
            }


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
            
        };
    } // namespace lcp_on_rlbwt
} // namespace stool