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

namespace stool
{
    namespace lcp_on_rlbwt
    {

        template <typename INDEX_SIZE, typename RLBWTDS>
        class FastSTNodeSubTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;

            uint64_t _stnode_count = 0;
            std::deque<INDEX_SIZE> childs_vec;
            std::deque<bool> first_child_flag_vec;
            std::deque<bool> maximal_repeat_check_vec;
            RLBWTDS *_RLBWTDS = nullptr;

        public:
            FastSTNodeSubTraverser()
            {
                maximal_repeat_check_vec.resize(0);
            }
            FastSTNodeSubTraverser(RLBWTDS *__RLBWTDS) : _RLBWTDS(__RLBWTDS)
            {
                maximal_repeat_check_vec.resize(0);
            }

        private:
            inline uint64_t get_child_start_position(uint64_t i) const
            {
                return this->childs_vec[(i * 4)];
            }
            inline uint64_t get_child_end_position(uint64_t i) const
            {
                return this->childs_vec[(i * 4) + 3];
            }

            
        public:
        };

    } // namespace lcp_on_rlbwt
} // namespace stool