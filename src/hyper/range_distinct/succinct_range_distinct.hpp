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
#include <sdsl/rmq_support.hpp> // include header for range minimum queries
namespace stool
{
    namespace lcp_on_rlbwt
    {

        template <typename INDEX_SIZE>
        class SuccinctRangeDistinctDataStructure
        {
        public:
            using XCHAR = uint64_t;
            const stool::WT *wt;
            uint8_t lastChar;
            uint64_t size;
            std::vector<XCHAR> cs;
            std::vector<uint64_t> cs1;
            std::vector<uint64_t> cs2;
            std::vector<uint64_t> C;
            void initialize(const stool::WT *_wt, const sdsl::int_vector<> *_bwt)
            {
                std::cout << "BWT = " << _bwt->size() << std::endl;
                wt = _wt;
                this->size = wt->size();
                lastChar = (*_bwt)[_bwt->size() - 1];
                uint64_t CHARMAX = UINT8_MAX + 1;

                cs.resize(CHARMAX, 0);
                cs1.resize(CHARMAX, 0);
                cs2.resize(CHARMAX, 0);

                stool::FMIndex::constructC(*_bwt, C);
            }
            uint64_t range_distinct(INDEX_SIZE i, INDEX_SIZE j, std::vector<CharInterval<INDEX_SIZE>> &output)
            {
                assert(i <= j);
                uint64_t k;
                uint64_t newJ = j + 1 == this->size ? this->size : j + 2;

                assert(i + 1 <= newJ);

                std::cout << i+1 << "/" << newJ << "/" << wt->size() << std::endl;
                sdsl::interval_symbols(*wt, i + 1, newJ, k, cs, cs1, cs2);

                bool b = j + 1 < this->size;
                for (uint64_t x = 0; x < k; x++)
                {
                    uint64_t fstRank = cs1[x] + 1;
                    uint64_t lastRank = cs2[x];

                    std::cout << "x = " << (cs[x]) << " [" << cs1[x] << "/" << cs2[x] << std::endl;

                    uint64_t left = wt->select(fstRank, cs[x]) - 1;

                    uint64_t right = wt->select(lastRank, cs[x]) - 1;
                    std::cout << "LR = " << (cs[x]) << " [" << left << "/" << right << std::endl;

                    if (j + 1 == this->size && cs[x] == lastChar)
                    {
                        //std::cout << "++" << std::endl;
                        right = this->size - 1;
                        b = true;
                    }
                    //assert((*bwt)[left] == (*bwt)[right]);
                    output[x] = CharInterval<INDEX_SIZE>(left, right, cs[x]);
                }
                if (!b)
                {
                    uint64_t left = this->size - 1;
                    uint64_t right = left;
                    
                    output[k++] = CharInterval<INDEX_SIZE>(left, right, lastChar);
                }

                return k;
            }
        };
    } // namespace lcp_on_rlbwt
} // namespace stool