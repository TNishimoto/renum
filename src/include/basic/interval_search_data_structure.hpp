#pragma once
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wt_gmr.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>

#include "./char_interval.hpp"
#include "./type.hpp"

namespace stool
{
    namespace renum
    {

    template <typename CHAR>
    class IntervalSearchDataStructure
    {
    public:
        std::vector<uint64_t> *C;
        stool::WT *wt;
        CHAR lastChar;
        std::vector<CHAR> cs;
        std::vector<uint64_t> cs1;
        std::vector<uint64_t> cs2;


        void initialize(stool::WT *_wt, std::vector<uint64_t> *_C, CHAR _lastChar);
        template <typename INDEX_SIZE>
        uint64_t getIntervals(INDEX_SIZE i, INDEX_SIZE j, std::vector<CharInterval<INDEX_SIZE, CHAR>> &output);
    };
}
} // namespace stool