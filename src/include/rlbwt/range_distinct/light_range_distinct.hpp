#pragma once

#include "../../basic/char_interval.hpp"
#include <sdsl/bit_vectors.hpp>

namespace stool
{
    namespace stnode_on_rlbwt
    {

        template <typename CHAR_VEC, typename INDEX_SIZE, typename CHAR>
        class LightRangeDistinctDataStructure
        {

            const CHAR_VEC *_char_vec;
            //using CHAR = typename CHAR_VEC::index_type;

            std::vector<int64_t> checker;
            std::vector<uint64_t> beginIndexes;
            std::vector<uint64_t> endIndexes;
            std::vector<uint64_t> charIndexes;
        public:
            void initialize(const CHAR_VEC *__char_vec);
            uint64_t range_distinct(INDEX_SIZE i, INDEX_SIZE j, std::vector<CharInterval<INDEX_SIZE, CHAR>> &output);
        };
    } // namespace stnode_on_rlbwt
} // namespace stool