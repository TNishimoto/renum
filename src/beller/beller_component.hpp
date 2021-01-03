#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>

#include "char_interval.hpp"
#include "../rlbwt/range_distinct/succinct_range_distinct.hpp"

//#include "sa_lcp.hpp"
using namespace std;
using namespace sdsl;

namespace stool
{
    namespace beller
    {

        template <typename INDEX>
        class OutputStructure
        {
        public:
            std::vector<stool::LCPInterval<INDEX>> outputIntervals;
            std::ofstream *out = nullptr;
            sdsl::bit_vector::rank_1_type *bwt_bit_rank1;
            bool is_output_maximal_substrings;
            uint64_t count = 0;
            uint64_t lcp_interval_count = 0;
            uint64_t peak = 0;

            void push(stool::LCPInterval<INDEX> &interval)
            {
                lcp_interval_count++;
                if (this->is_output_maximal_substrings)
                {
                    uint64_t k1 = interval.i == 0 ? 0 : (*bwt_bit_rank1)(interval.i);
                    uint64_t k2 = (*bwt_bit_rank1)(interval.j + 1);


                    bool isNotMaximalRepeat = ((k2 - k1 == 0) || ((k2 - k1) == (interval.j - interval.i + 1)));

                    if (!isNotMaximalRepeat)
                    {
                        this->outputIntervals.push_back(interval);
                        count++;
                    }
                }
                else
                {
                    this->outputIntervals.push_back(interval);
                        count++;

                }

                if (this->out != nullptr && this->outputIntervals.size() > 8192)
                {
                    out->write(reinterpret_cast<const char *>(&this->outputIntervals[0]), sizeof(stool::LCPInterval<INDEX>) * this->outputIntervals.size());
                    this->outputIntervals.clear();
                }
            }
            void finish()
            {
                if (this->out != nullptr && this->outputIntervals.size() > 0)
                {
                    out->write(reinterpret_cast<const char *>(&this->outputIntervals[0]), sizeof(stool::LCPInterval<INDEX>) * this->outputIntervals.size());
                    this->outputIntervals.clear();
                }

            }
        };

        struct BellerSmallComponent
        {
            uint64_t last_idx;
            uint64_t last_lb;
            bool occB;
            std::set<uint8_t> nextOccurrenceSet;

            BellerSmallComponent()
            {
            }
            void initialize()
            {
                last_idx = UINT64_MAX;
                last_lb = UINT64_MAX;
                occB = false;
            }
        };

        template <typename INDEX>
        class BellerComponent
        {
        public:
            using INTERVAL = stool::LCPInterval<INDEX>;

            std::vector<std::queue<INTERVAL>> intervalQueues;
            std::vector<bool> checker;
            std::vector<uint64_t> counter;
            std::vector<uint8_t> occurrenceChars;
            std::vector<CharInterval<INDEX>> charIntervalTmpVec;

            uint64_t lcp = 0;
            uint64_t debugCounter = 0;

            void initialize(IntervalSearchDataStructure &range)
            {
                uint64_t bwtSize = range.wt->size();
                uint64_t CHARMAX = UINT8_MAX + 1;

                intervalQueues.resize(CHARMAX);
                counter.resize(CHARMAX, 0);
                checker.resize(bwtSize + 1, false);
                checker[0] = false;
                occurrenceChars.resize(0);
                charIntervalTmpVec.resize(CHARMAX);

                this->first_process(range);
            }

            void process(IntervalSearchDataStructure &range, BellerSmallComponent &bsc, stool::beller::OutputStructure<INDEX> &os)
            {
                uint64_t k = 0;
                for (auto &c : this->occurrenceChars)
                {
                    auto &que = this->intervalQueues[c];
                    uint64_t queSize = this->counter[c];

                    while (queSize > 0)
                    {
                        bsc.occB = true;

                        auto top = que.front();
                        que.pop();
                        queSize--;

                        if (!this->checker[top.j + 1])
                        {
                            if (bsc.last_lb == UINT64_MAX)
                            {
                                bsc.last_lb = top.i;
                            }

                            this->checker[top.j + 1] = true;
                            bsc.last_idx = top.j + 1;

                            uint64_t charIntvCount = range.getIntervals(top.i, top.j, this->charIntervalTmpVec);

                            this->debugCounter++;

                            for (uint64_t i = 0; i < charIntvCount; i++)
                            {
                                auto &intv = this->charIntervalTmpVec[i];
                                this->intervalQueues[intv.c].push(INTERVAL(intv.i, intv.j, top.lcp + 1));
                                k++;
                                bsc.nextOccurrenceSet.insert(intv.c);
                            }
                        }
                        else
                        {
                            this->checker[top.j + 1] = true;
                            if (top.i == bsc.last_idx)
                            {
                                if (os.out != nullptr)
                                {
                                    INTERVAL iv(bsc.last_lb, top.j, top.lcp - 1);
                                    os.push(iv);
                                }

                                bsc.last_lb = UINT64_MAX;
                                bsc.last_idx = UINT64_MAX;

                                uint64_t charIntvCount = range.getIntervals(top.i, top.j, this->charIntervalTmpVec);
                                this->debugCounter++;

                                for (uint64_t i = 0; i < charIntvCount; i++)
                                {
                                    auto &intv = this->charIntervalTmpVec[i];
                                    this->intervalQueues[intv.c].push(INTERVAL(intv.i, intv.j, top.lcp + 1));
                                    k++;

                                    bsc.nextOccurrenceSet.insert(intv.c);
                                }
                            }
                        }
                    }
                }
                if(os.peak < k){
                    os.peak = k;
                }
            }
            void computeLCPIntervals(IntervalSearchDataStructure &range, bool &isEnd, stool::beller::OutputStructure<INDEX> &os)
            {
                using INTERVAL = stool::LCPInterval<INDEX>;
                uint64_t max_interval_count = 0;

                std::vector<INTERVAL> output_lcp_intervals;

                BellerSmallComponent bsc;
                bsc.initialize();

                uint64_t interval_count = 0;

                this->lcp++;
                for (uint64_t x = 0; x < this->intervalQueues.size(); x++)
                {
                    auto &que = this->intervalQueues[x];
                    uint64_t queSize = que.size();
                    this->counter[x] = queSize;
                    interval_count += queSize;
                }
                if (max_interval_count < interval_count)
                    max_interval_count = interval_count;

                this->process(range, bsc, os);

                this->occurrenceChars.resize(0);
                for (auto c : bsc.nextOccurrenceSet)
                {
                    this->occurrenceChars.push_back(c);
                }
                if (!bsc.occB)
                {
                    this->checker.pop_back();
                    isEnd = true;
                    os.finish();
                }

                //return output_lcp_intervals;
            }

        private:
            void first_process(IntervalSearchDataStructure &range)
            {
                uint64_t bwtSize = range.wt->size();
                INTERVAL fst(0, bwtSize - 1, 0);
                uint64_t charIntvCount = range.getIntervals(fst.i, fst.j, this->charIntervalTmpVec);

                this->debugCounter++;

                std::set<uint8_t> nextOccurrenceSet;

                for (uint64_t i = 0; i < charIntvCount; i++)
                {
                    auto &intv = this->charIntervalTmpVec[i];
                    this->intervalQueues[intv.c].push(INTERVAL(intv.i, intv.j, 1));
                    nextOccurrenceSet.insert(intv.c);
                }

                for (auto c : nextOccurrenceSet)
                {
                    this->occurrenceChars.push_back(c);
                }
            }
        };

    } // namespace beller
} // namespace stool