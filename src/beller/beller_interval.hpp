#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>

#include "char_interval.hpp"

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

            void initialize(int_vector<> &bwt)
            {
                uint64_t CHARMAX = UINT8_MAX + 1;

                queArr.resize(CHARMAX);
                counter.resize(CHARMAX, 0);
                checker.resize(bwt.size() + 1, false);
                checker[0] = false;
                occurrenceChars.resize(0);
                charIntervalTmpVec.resize(CHARMAX);
            }
        };
        template <typename INDEX>
        std::vector<stool::LCPInterval<INDEX>> computeLCPIntervals(int_vector<> &bwt, IntervalSearchDataStructure &range, BellerComponent<INDEX> &comp, bool &isEnd)
        {
            using INTERVAL = stool::LCPInterval<INDEX>;
            uint64_t n = bwt.size();
            //uint8_t lastChar = bwt[bwt.size() - 1];
            uint64_t max_interval_count = 0;

            std::vector<INTERVAL> r;
            //std::cout << "bwt: " << n << std::endl;

            if (comp.lcp == 0)
            {
                INTERVAL fst(0, n - 1, 0);
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
                            r.push_back(INTERVAL(last_lb, top.j, top.lcp - 1));

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

            return r;

        } // namespace beller
        template <typename INDEX>
        std::vector<stool::LCPInterval<INDEX>> computeLCPIntervals(int_vector<> &bwt)
        {
            IntervalSearchDataStructure range;
            range.initialize(bwt);
            stool::beller::BellerComponent<INDEX> comp;
            comp.initialize(bwt);
            bool isEnd = false;
            std::vector<stool::LCPInterval<INDEX>> r;

            while (!isEnd)
            {
                auto r2 = computeLCPIntervals(bwt, range, comp, isEnd);
                for (auto it : r2)
                {
                    r.push_back(it);
                }
            }
            return r;
        }
        template <typename INDEX>
        std::vector<stool::LCPInterval<INDEX>> computeMaximalSubstrings(int_vector<> &bwt, IntervalSearchDataStructure &range, BellerComponent<INDEX> &comp, bool &isEnd)
        {
            std::vector<stool::LCPInterval<INDEX>> r;
            auto r2 = computeLCPIntervals(bwt, range, comp, isEnd);
            for (auto it : r2)
            {
                uint8_t fstChar = bwt[it.i];
                uint8_t lstChar = bwt[it.j];
                if (fstChar == lstChar)
                {
                    uint64_t p1 = range.wt.rank(it.i + 1, fstChar);
                    uint64_t p2 = range.wt.rank(it.j + 1, fstChar);
                    //std::cout << it.to_string() << std::endl;
                    //std::cout << p1 << "/" << p2 << std::endl;

                    if (p2 - p1 != it.j - it.i)
                    {
                        r.push_back(it);
                    }
                }
                else
                {
                    r.push_back(it);
                }
            }
            return r;
        }
        template <typename INDEX>
        uint64_t outputMaximalSubstrings(int_vector<> &bwt, std::ofstream &out)
        {
            IntervalSearchDataStructure range;
            range.initialize(bwt);

            stool::beller::BellerComponent<INDEX> comp;
            comp.initialize(bwt);
            bool isEnd = false;
            uint64_t count = 0;
            while (!isEnd)
            {

                auto r2 = computeMaximalSubstrings(bwt, range, comp, isEnd);
                for (auto &it : r2)
                {
                    out.write(reinterpret_cast<const char *>(&it), sizeof(stool::LCPInterval<INDEX>));
                }

                count += r2.size();
            }

                std::cout << "Range Distinct Count = " << comp.debugCounter << "/" << bwt.size() << std::endl;
            uint64_t dx = range.wt.select(1, 0) - 1;
            auto last = stool::LCPInterval<INDEX>(dx, dx, bwt.size());
            out.write(reinterpret_cast<const char *>(&last), sizeof(stool::LCPInterval<INDEX>));
            count += 1;
            return count;
        }

    } // namespace beller
} // namespace stool