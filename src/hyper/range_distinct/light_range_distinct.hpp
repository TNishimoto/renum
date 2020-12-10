#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>

#include "stool/src/debug.hpp"
#include "stool/src/elias_fano_vector.hpp"
#include "../../beller/char_interval.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {

        template <typename CHAR_VEC, typename INDEX_SIZE>
        class LightRangeDistinctDataStructure
        {
        public:
            const CHAR_VEC *_char_vec;

            std::vector<int64_t> checker;
            std::vector<uint64_t> beginIndexes;
            std::vector<uint64_t> endIndexes;
            std::vector<uint64_t> charIndexes;
            void preprocess(const CHAR_VEC *__char_vec)
            {
                int32_t charMaxSize = ((int32_t)UINT8_MAX) + 1;
                this->_char_vec = __char_vec;

                checker.resize(charMaxSize, -1);
                beginIndexes.resize(charMaxSize, 0);
                endIndexes.resize(charMaxSize, 0);
                charIndexes.resize(charMaxSize, 0);
            }
            uint64_t range_distinct(INDEX_SIZE i, INDEX_SIZE j, std::vector<CharInterval<INDEX_SIZE>> &output)
            {
                uint64_t count = 0;
                assert(i <= j);
                for (uint64_t x = i; x <= j; x++)
                {
                    uint8_t c = (*this->_char_vec)[x];
                    if (checker[c] == -1)
                    {
                        checker[c] = count;
                        beginIndexes[count] = x;
                        endIndexes[count] = x;
                        charIndexes[count] = c;
                        count++;
                    }
                    else
                    {
                        endIndexes[checker[c]] = x;
                    }
                }
                for (uint64_t x = 0; x < count; x++)
                {
                    output[x] = CharInterval<INDEX_SIZE>(beginIndexes[x], endIndexes[x], charIndexes[x]);
                    checker[charIndexes[x]] = -1;
                }
                assert(count > 0);
                return count;
            }
        };
    } // namespace lcp_on_rlbwt
} // namespace stool