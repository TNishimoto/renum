#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <mutex>
#include <thread>

#include "../../include/rlbwt/range_distinct/light_range_distinct.hpp"

namespace stool
{
    namespace stnode_on_rlbwt
    {
        template <typename CHAR_VEC, typename INDEX_SIZE, typename CHAR>
        void LightRangeDistinctDataStructure<CHAR_VEC, INDEX_SIZE, CHAR>::initialize(const CHAR_VEC *__char_vec)
        {
            int32_t charMaxSize = ((int32_t)UINT8_MAX) + 1;
            this->_char_vec = __char_vec;

            checker.resize(charMaxSize, -1);
            beginIndexes.resize(charMaxSize, 0);
            endIndexes.resize(charMaxSize, 0);
            charIndexes.resize(charMaxSize, 0);
        }

        template <typename CHAR_VEC, typename INDEX_SIZE, typename CHAR>
        uint64_t LightRangeDistinctDataStructure<CHAR_VEC, INDEX_SIZE, CHAR>::range_distinct(INDEX_SIZE i, INDEX_SIZE j, std::vector<CharInterval<INDEX_SIZE, CHAR>> &output)
        {
            //std::lock_guard<std::mutex> lock(std::mutex);

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
                output[x] = CharInterval<INDEX_SIZE, uint8_t>(beginIndexes[x], endIndexes[x], charIndexes[x]);
                checker[charIndexes[x]] = -1;
            }
            assert(count > 0);
            return count;
        }

        template class LightRangeDistinctDataStructure<sdsl::int_vector<>, uint32_t, uint8_t>;
        template class LightRangeDistinctDataStructure<sdsl::int_vector<>, uint64_t, uint8_t>;

    } // namespace stnode_on_rlbwt
} // namespace stool