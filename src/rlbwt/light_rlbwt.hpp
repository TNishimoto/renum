#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>

namespace stool
{
    namespace lcp_on_rlbwt
    {
        /*
        template <typename CHAR, typename CHARVEC, typename POSVEC>
        class LightRLBWT
        {
        public:
            const CHARVEC &char_vec;
            const POSVEC &lpos_vec;

            LightRLBWT(const CHARVEC &_char_vec, const POSVEC &_lpos_vec) : char_vec(_char_vec), lpos_vec(_lpos_vec){

            }

            uint64_t get_lindex_containing_the_position(uint64_t lposition) const
            {
                auto p = std::upper_bound(this->lpos_vec.begin(), this->lpos_vec.end(), lposition);
                uint64_t pos = std::distance(this->lpos_vec.begin(), p) - 1;
                return pos;
            }
            uint64_t str_size() const
            {
                return lpos_vec[lpos_vec.size() - 1];
            }
            uint8_t get_char_by_run_index(uint64_t _run_index) const
            {
                return char_vec[_run_index];
            }
            uint64_t rle_size() const
            {
                return char_vec.size();
            }
            uint64_t get_run(uint64_t i) const
            {
                return lpos_vec[(i + 1)] - lpos_vec[i];
            }
            uint64_t get_lpos(uint64_t i) const
            {
                return lpos_vec[i];
            }

            uint64_t get_end_rle_lposition() const
            {
                for (uint64_t i = 0; i < char_vec.size(); i++)
                {
                    if (char_vec[i] == 0)
                    {
                        return i;
                    }
                }
                return std::numeric_limits<uint64_t>::max();
            }
        };
        */
    } // namespace lcp_on_rlbwt
} // namespace stool