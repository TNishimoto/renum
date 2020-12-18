#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>

#include "char_interval.hpp"
#include "../hyper/range_distinct/succinct_range_distinct.hpp"

//#include "sa_lcp.hpp"
using namespace std;
using namespace sdsl;

namespace stool
{
    namespace beller
    {
        template <typename INDEX>
        class BellerComponent
        {
        public:
            using INTERVAL = stool::LCPInterval<INDEX>;

            std::vector<std::queue<INTERVAL>> queArr;
            std::vector<bool> checker;
            std::vector<uint64_t> counter;
            std::vector<uint8_t> occurrenceChars;
            std::vector<CharInterval<INDEX>> charIntervalTmpVec;

            uint64_t lcp = 0;
            uint64_t debugCounter = 0;

            void initialize(uint64_t bwtSize)
            {
                uint64_t CHARMAX = UINT8_MAX + 1;

                queArr.resize(CHARMAX);
                counter.resize(CHARMAX, 0);
                checker.resize(bwtSize + 1, false);
                checker[0] = false;
                occurrenceChars.resize(0);
                charIntervalTmpVec.resize(CHARMAX);
            }
        };
        template <typename INDEX>
        std::vector<stool::LCPInterval<INDEX>> computeLCPIntervals(IntervalSearchDataStructure &range, BellerComponent<INDEX> &comp, bool &isEnd)
        {
            /*   
            IntervalSearchDataStructure range2;
            range2.initialize(bwt);
            std::vector<CharInterval<INDEX>> charIntervalTmpVec;
                charIntervalTmpVec.resize(256);
                */
            uint64_t bwtSize = range.wt->size();
            using INTERVAL = stool::LCPInterval<INDEX>;
            //uint64_t n = bwt.size();
            //uint8_t lastChar = bwt[bwt.size() - 1];
            uint64_t max_interval_count = 0;

            std::vector<INTERVAL> output_lcp_intervals;
            //std::cout << "bwt: " << n << std::endl;

            if (comp.lcp == 0)
            {
                INTERVAL fst(0, bwtSize - 1, 0);
                uint64_t charIntvCount = range.getIntervals(fst.i, fst.j, comp.charIntervalTmpVec);

                comp.debugCounter++;

                std::set<uint8_t> nextOccurrenceSet;

                for (uint64_t i = 0; i < charIntvCount; i++)
                {
                    auto &intv = comp.charIntervalTmpVec[i];
                    comp.queArr[intv.c].push(INTERVAL(intv.i, intv.j, 1));
                    nextOccurrenceSet.insert(intv.c);
                }

                for (auto c : nextOccurrenceSet)
                {
                    comp.occurrenceChars.push_back(c);
                }
            }

            //comp.lcp = 0;
            uint64_t interval_count = 0;
            uint64_t last_idx = UINT64_MAX;
            uint64_t last_lb = UINT64_MAX;

            bool occB = false;
            comp.lcp++;
            for (uint64_t x = 0; x < comp.queArr.size(); x++)
            {
                auto &que = comp.queArr[x];
                uint64_t queSize = que.size();
                comp.counter[x] = queSize;
                interval_count += queSize;
            }
            if (max_interval_count < interval_count)
                max_interval_count = interval_count;
            //std::cout << "lcp = " << comp.lcp << ", count = " << interval_count << ", max = " << max_interval_count << std::endl;

            std::set<uint8_t> nextOccurrenceSet;
            for (auto &c : comp.occurrenceChars)
            {
                auto &que = comp.queArr[c];
                uint64_t queSize = comp.counter[c];

                while (queSize > 0)
                {
                    occB = true;

                    auto top = que.front();
                    que.pop();
                    queSize--;

                    //std::cout << top.to_string() << ",checker[" << (top.j + 1) << "]=" << checker[top.j + 1] << ", last_idx = " << last_idx << std::endl;

                    if (!comp.checker[top.j + 1])
                    {
                        if (last_lb == UINT64_MAX)
                        {
                            last_lb = top.i;
                        }

                        comp.checker[top.j + 1] = true;
                        last_idx = top.j + 1;
                        //auto tmp = getIntervals(top.i, top.j, bwt, C, wt);

                        uint64_t charIntvCount = range.getIntervals(top.i, top.j, comp.charIntervalTmpVec);
                        /*
                        uint64_t charIntvCount2 = range2.getIntervals(top.i, top.j, charIntervalTmpVec);

                        for (uint64_t i = 0; i < charIntvCount; i++)
                        {
                            auto &intv = comp.charIntervalTmpVec[i];
                            INTERVAL test1(intv.i, intv.j, top.lcp + 1);
                            std::cout << test1.to_string();
                        }
                        std::cout << std::endl;
                        for (uint64_t i = 0; i < charIntvCount2; i++)
                        {
                            auto &intv = charIntervalTmpVec[i];
                            INTERVAL test2(intv.i, intv.j, top.lcp + 1);
                            std::cout << test2.to_string();
                        }
                        std::cout << std::endl;
                        throw -1;
                        */

                        comp.debugCounter++;

                        //check(tmp, tmpx);
                        for (uint64_t i = 0; i < charIntvCount; i++)
                        {
                            auto &intv = comp.charIntervalTmpVec[i];
                            comp.queArr[intv.c].push(INTERVAL(intv.i, intv.j, top.lcp + 1));
                            nextOccurrenceSet.insert(intv.c);
                        }
                    }
                    else
                    {
                        comp.checker[top.j + 1] = true;
                        if (top.i == last_idx)
                        {
                            output_lcp_intervals.push_back(INTERVAL(last_lb, top.j, top.lcp - 1));

                            last_lb = UINT64_MAX;
                            last_idx = UINT64_MAX;

                            //auto tmp = getIntervals(top.i, top.j, bwt, C, wt);
                            uint64_t charIntvCount = range.getIntervals(top.i, top.j, comp.charIntervalTmpVec);
                            comp.debugCounter++;

                            //check(tmp, tmpx);

                            for (uint64_t i = 0; i < charIntvCount; i++)
                            {
                                auto &intv = comp.charIntervalTmpVec[i];
                                comp.queArr[intv.c].push(INTERVAL(intv.i, intv.j, top.lcp + 1));
                                nextOccurrenceSet.insert(intv.c);
                            }
                        }
                    }
                }
            }
            comp.occurrenceChars.resize(0);
            for (auto c : nextOccurrenceSet)
            {
                comp.occurrenceChars.push_back(c);
            }
            if (!occB)
            {
                comp.checker.pop_back();
                isEnd = true;
            }

            return output_lcp_intervals;

        } // namespace beller
        template <typename INDEX>
        std::vector<stool::LCPInterval<INDEX>> computeLCPIntervals(IntervalSearchDataStructure &range)
        {
            stool::beller::BellerComponent<INDEX> comp;
            comp.initialize(range.wt->size());
            bool isEnd = false;
            std::vector<stool::LCPInterval<INDEX>> r;

            while (!isEnd)
            {
                auto r2 = computeLCPIntervals(range, comp, isEnd);
                for (auto it : r2)
                {
                    r.push_back(it);
                }
            }
            return r;
        }
        /*
        template <typename INDEX>
        std::vector<stool::LCPInterval<INDEX>> computeLCPIntervals2(wt_huff<> &wt)
        {
            IntervalSearchDataStructure range;
            range.initialize(bwt);
            stool::beller::BellerComponent<INDEX> comp;
            comp.initialize(bwt);
            bool isEnd = false;
            std::vector<stool::LCPInterval<INDEX>> r;

            uint64_t size = bwt.size();
            while (!isEnd)
            {
                auto r2 = computeLCPIntervals(size, range, comp, isEnd);
                for (auto it : r2)
                {
                    r.push_back(it);
                }
            }
            return r;
        }
        */

        template <typename INDEX>
        std::vector<stool::LCPInterval<INDEX>> computeMaximalSubstrings(IntervalSearchDataStructure &range, BellerComponent<INDEX> &comp, bool &isEnd, uint8_t lastChar, sdsl::bit_vector::rank_1_type &bwt_bit_rank1)
        {
            std::vector<stool::LCPInterval<INDEX>> r;
            //uint64_t size = range.wt->size();

            auto r2 = computeLCPIntervals(range, comp, isEnd);
            for (auto it : r2)
            {
                uint64_t k1 = it.i == 0 ? 0 : bwt_bit_rank1(it.i);
                uint64_t k2 = bwt_bit_rank1(it.j + 1);
                bool b = !((k2 - k1 == 0) || ((k2 - k1) == (it.j - it.i + 1)));
                if (b)
                {
                    r.push_back(it);
                }
                /*
                uint8_t fstChar = it.i == size - 1 ? lastChar : (*range.wt)[it.i + 1];
                uint8_t lstChar = it.j == size - 1 ? lastChar : (*range.wt)[it.j + 1];
                if (fstChar == lstChar)
                {
                    uint64_t p1 = range.wt->rank(it.i + 1, fstChar);
                    uint64_t p2 = range.wt->rank(it.j + 1, fstChar);
                    

                    //std::cout << it.to_string() << std::endl;
                    //std::cout << p1 << "/" << p2 << std::endl;
                }
                else
                {
                    if (!b)
                    {
                        throw -1;
                    }
                    r.push_back(it);
                }
                */
            }
            return r;
        }
        template <typename INDEX>
        uint64_t outputMaximalSubstrings(IntervalSearchDataStructure &range, std::ofstream &out, uint8_t lastChar, sdsl::bit_vector::rank_1_type &bwt_bit_rank1)
        {
            stool::beller::BellerComponent<INDEX> comp;
            uint64_t size = range.wt->size();
            comp.initialize(size);
            bool isEnd = false;
            uint64_t count = 0;
            while (!isEnd)
            {

                auto r2 = computeMaximalSubstrings(range, comp, isEnd, lastChar, bwt_bit_rank1);
                for (auto &it : r2)
                {
                    out.write(reinterpret_cast<const char *>(&it), sizeof(stool::LCPInterval<INDEX>));
                }

                count += r2.size();
            }

            std::cout << "Range Distinct Count = " << comp.debugCounter << "/" << size << std::endl;

            uint64_t dx = range.wt->select(1, 0) - 1;
            auto last = stool::LCPInterval<INDEX>(dx, dx, size);
            out.write(reinterpret_cast<const char *>(&last), sizeof(stool::LCPInterval<INDEX>));

            count += 1;
            return count;
        }

    } // namespace beller
} // namespace stool