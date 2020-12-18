#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>

//#include "../weiner_interval.hpp"
#include "range_distinct.hpp"
#include "succinct_range_distinct.hpp"
#include "./light_range_distinct.hpp"
#include "../../rlbwt/rinterval.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {   
        /*
        template <typename INDEX_SIZE>
        class RangeDistinctDataStructureOnRLBWT
        {
        public:

            using RINTERVAL = RInterval<INDEX_SIZE>;
            LightRangeDistinctDataStructure<sdsl::int_vector<>, INDEX_SIZE> light_srdds;
            SuccinctRangeDistinctDataStructure<INDEX_SIZE> srdds;
            uint64_t total_cover1 = 0;
            uint64_t num1 = 0;

            uint64_t total_cover2 = 0;
            uint64_t num2 = 0;

            void initialize(const stool::WT *_wt, const sdsl::int_vector<> *_bwt)
            {

                uint64_t CHARMAX = UINT8_MAX + 1;
                light_srdds.preprocess(_bwt);
                srdds.initialize(_wt, _bwt);
            }

            bool check(RInterval<INDEX_SIZE> &range)
            {
                uint64_t CHARMAX = UINT8_MAX + 1;

                std::vector<CharInterval<INDEX_SIZE>> DEBUGcharIntervalTmpVec1;
                std::vector<CharInterval<INDEX_SIZE>> DEBUGcharIntervalTmpVec2;
                DEBUGcharIntervalTmpVec1.resize(CHARMAX);
                DEBUGcharIntervalTmpVec2.resize(CHARMAX);

                assert(range.beginIndex <= range.endIndex);
                uint64_t count1 = srdds.range_distinct(range.beginIndex, range.endIndex, DEBUGcharIntervalTmpVec1);
                uint64_t count2 = light_srdds.range_distinct(range.beginIndex, range.endIndex, DEBUGcharIntervalTmpVec2);
                if(count1 != count2){
                    std::cout << "count distinct" << std::endl;
                    throw -1;
                }else{
                    DEBUGcharIntervalTmpVec1.resize(count1);
                    DEBUGcharIntervalTmpVec2.resize(count2);
                    stool::beller::check(DEBUGcharIntervalTmpVec1, DEBUGcharIntervalTmpVec2);

                }
                return true;

            }
            uint64_t range_distinct(uint64_t l, uint64_t r, std::vector<CharInterval<INDEX_SIZE>> &charIntervalTmpVec)
            {
                uint64_t count = 0;
                assert(l <= r);
                //assert(check(range));

                if (false)
                //if (l - r >= 16)
                {
                    count = srdds.range_distinct(l, r, charIntervalTmpVec);

                    total_cover1 += r - l + 1;
                    num1++;
                }
                else
                {
                    count = light_srdds.range_distinct(l, r, charIntervalTmpVec);
                    //std::cout << "+" << count << std::endl;

                    total_cover2 += r - l + 1;
                    num2++;
                }
                assert(count > 0);
                return count;
            }
        };
        */

    } // namespace lcp_on_rlbwt
} // namespace stool