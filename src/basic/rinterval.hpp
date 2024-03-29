#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>

namespace stool
{
    namespace stnode_on_rlbwt
    {
        template <typename INDEX_SIZE>
        struct RInterval
        {
            INDEX_SIZE beginIndex;
            INDEX_SIZE beginDiff;
            INDEX_SIZE endIndex;
            INDEX_SIZE endDiff;

            void print() const;

            bool is_special() const;

            static RInterval get_special();
        };

    } // namespace stnode_on_rlbwt
} // namespace stool