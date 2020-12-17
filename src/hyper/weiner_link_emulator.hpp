#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>

namespace stool
{
    namespace lcp_on_rlbwt
    {
        uint64_t DEBUG_LIMIT = 1000;
        /*
            This is a data structure to ...
        */
        template <typename INDEX_SIZE, typename RLBWTDS>
        class ExplicitWeinerLinkEmulator
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;

        public:
            std::vector<std::vector<RINTERVAL>> childrenVec;
            std::vector<RINTERVAL> stnodeVec;
            std::vector<uint64_t> indexVec;
            std::vector<bool> stnodeOccFlagArray;

            uint64_t indexCount = 0;
            uint64_t explicitChildCount = 0;
            uint64_t range_distinct_threshold = 16;

            // For range distinct
            std::vector<uint8_t> charTmpVec;
            vector<RINTERVAL> rIntervalTmpVec;
            std::vector<CharInterval<INDEX_SIZE>> charIntervalTmpVec;

            RLBWTDS *_RLBWTDS;
            LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE> *lightDS = nullptr;
            SuccinctRangeDistinctDataStructure<INDEX_SIZE> *heavyDS = nullptr;

            void initialize(RLBWTDS *_rlbwtds)
            {
                _RLBWTDS = _rlbwtds;
                uint64_t CHARMAX = UINT8_MAX + 1;
                childrenVec.resize(CHARMAX);
                indexVec.resize(CHARMAX);

                stnodeOccFlagArray.resize(CHARMAX, false);
                stnodeVec.resize(CHARMAX);

                rIntervalTmpVec.resize(CHARMAX);
                charTmpVec.resize(CHARMAX);

                charIntervalTmpVec.resize(CHARMAX);
            }
            uint64_t get_explicit_stnode_count()
            {
                return this->indexCount;
            }
            uint64_t get_explicit_children_count()
            {
                uint64_t k = 0;
                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto &it = this->indexVec[i];
                    k += childrenVec[it].size();
                }
                return k;
            }

            void clear()
            {
                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto &it = this->indexVec[i];
                    childrenVec[it].clear();

                    stnodeOccFlagArray[it] = false;
                }
                indexCount = 0;
                explicitChildCount = 0;
            }
            /*
            void move_st_internal_nodes(std::deque<uint64_t> &outputExplicitChildrenVec, std::deque<bool> &leftmost_child_bits)
            {
                for (uint64_t i = 0; i < this->indexCount; i++)
                {

                    auto character = this->indexVec[i];

                    auto &currentVec = this->childrenVec[character];
                    uint64_t count = this->childrenVec[character].size();
                    leftmost_child_bits.push_back(true);
                    for (uint64_t i = 1; i < count; i++)
                    {
                        leftmost_child_bits.push_back(false);
                    }

                    for (uint64_t j = 0; j < count; j++)
                    {
                    outputExplicitChildrenVec.push_back(currentVec[j].beginIndex);
                    outputExplicitChildrenVec.push_back(currentVec[j].beginDiff);
                    outputExplicitChildrenVec.push_back(currentVec[j].endIndex);
                    outputExplicitChildrenVec.push_back(currentVec[j].endDiff);
                        //output.push_weiner(currentVec[j]);
                    }
                }
            }
            */
            void move_st_internal_nodes(std::deque<uint64_t> &outputExplicitChildrenVec, std::deque<bool> &leftmost_child_bits, std::deque<bool> &maximal_repeat_check_vec, uint8_t c)
            {
                auto &currentVec = this->childrenVec[c];
                uint64_t count = currentVec.size();
                leftmost_child_bits.push_back(true);
                for (uint64_t i = 1; i < count; i++)
                {
                    leftmost_child_bits.push_back(false);
                }
                uint64_t st_left = UINT64_MAX;
                uint64_t st_right = 0;

                for (uint64_t j = 0; j < count; j++)
                {

                    uint64_t left = this->_RLBWTDS->get_fpos(currentVec[j].beginIndex, currentVec[j].beginDiff);
                    uint64_t right = this->_RLBWTDS->get_fpos(currentVec[j].endIndex, currentVec[j].endDiff);

                    if (left < st_left)
                    {
                        st_left = left;
                    }
                    if (right > st_right)
                    {
                        st_right = right;
                    }

                    outputExplicitChildrenVec.push_back(left);
                    outputExplicitChildrenVec.push_back(right);
                }

                uint64_t x = this->_RLBWTDS->get_lindex_containing_the_position(st_left);
                uint64_t d = this->_RLBWTDS->get_run(x);
                bool isMaximalRepeat = (this->_RLBWTDS->get_lpos(x) + d) <= st_right;
                maximal_repeat_check_vec.push_back(isMaximalRepeat);
            }

            bool pushExplicitWeinerInterval(const RINTERVAL &w, uint8_t c)
            {
                auto &lcpIntv = this->stnodeVec[c];
                bool isLastChild = lcpIntv.endIndex == w.endIndex && lcpIntv.endDiff == w.endDiff;
                bool isFirstChild = lcpIntv.beginIndex == w.beginIndex && lcpIntv.beginDiff == w.beginDiff;
                bool b = !(isFirstChild && isLastChild);

                if (b)
                {

                    if (this->childrenVec[c].size() == 0)
                    {

                        this->indexVec[this->indexCount] = c;
                        this->indexCount++;
                    }
                    this->childrenVec[c].push_back(w);
                    explicitChildCount++;
                }
                return b;
            }
            void pushLCPInterval(const RINTERVAL &w, uint8_t c)
            {
                this->stnodeVec[c] = w;
                this->stnodeOccFlagArray[c] = true;
            }

            //template <typename LPOSDS, typename RANGEDS>
            void computeSTNodeCandidates(const RINTERVAL &w)
            {
                RINTERVAL frontL = this->_RLBWTDS->getIntervalOnL(w);
                uint64_t resultCount = this->range_distinct(frontL);
                //this->_RLBWTDS->rangeOnRLBWT.range_distinct(frontL);
                for (uint64_t i = 0; i < resultCount; i++)
                {
                    typename RLBWTDS::UCHAR c = this->charTmpVec[i];
                    auto &it = this->rIntervalTmpVec[i];

                    this->pushLCPInterval(it, c);
                }
            }
            void computeSTNodeCandidates2(const RINTERVAL &w)
            {
                uint64_t resultCount = this->range_distinct(w);
                for (uint64_t i = 0; i < resultCount; i++)
                {
                    typename RLBWTDS::UCHAR c = this->charTmpVec[i];
                    auto &it = this->rIntervalTmpVec[i];

                    this->pushLCPInterval(it, c);
                }
            }

        public:
            void computeSTChildren(const RINTERVAL &w)
            {
                RINTERVAL frontL = this->_RLBWTDS->getIntervalOnL(w);

                assert(frontL.beginIndex <= frontL.endIndex);
                uint64_t resultCount = this->range_distinct(frontL);

                for (uint64_t i = 0; i < resultCount; i++)
                {
                    typename RLBWTDS::UCHAR c = this->charTmpVec[i];
                    auto &it = this->rIntervalTmpVec[i];

                    //auto &lcpIntv = this->stnodeVec[c];
                    //bool isLastChild = lcpIntv.endIndex == it.endIndex && lcpIntv.endDiff == it.endDiff;
                    //if(!isLastChild){
                    this->pushExplicitWeinerInterval(it, c);
                    //}
                }
            }
            void computeSTChildren2(const RINTERVAL &w)
            {

                assert(w.beginIndex <= w.endIndex);
                uint64_t resultCount = this->range_distinct(w);

                for (uint64_t i = 0; i < resultCount; i++)
                {
                    typename RLBWTDS::UCHAR c = this->charTmpVec[i];
                    auto &it = this->rIntervalTmpVec[i];

                    //auto &lcpIntv = this->stnodeVec[c];
                    //bool isLastChild = lcpIntv.endIndex == it.endIndex && lcpIntv.endDiff == it.endDiff;
                    //if(!isLastChild){
                    this->pushExplicitWeinerInterval(it, c);
                    //}
                }
            }
#if DEBUG
            bool check(const RInterval<INDEX_SIZE> &range)
            {
                uint64_t CHARMAX = UINT8_MAX + 1;

                std::vector<CharInterval<INDEX_SIZE>> DEBUGcharIntervalTmpVec1;
                std::vector<CharInterval<INDEX_SIZE>> DEBUGcharIntervalTmpVec2;
                DEBUGcharIntervalTmpVec1.resize(CHARMAX);
                DEBUGcharIntervalTmpVec2.resize(CHARMAX);

                assert(range.beginIndex <= range.endIndex);
                uint64_t count1 = heavyDS->range_distinct(range.beginIndex, range.endIndex, DEBUGcharIntervalTmpVec1);
                uint64_t count2 = lightDS->range_distinct(range.beginIndex, range.endIndex, DEBUGcharIntervalTmpVec2);
                if (count1 != count2)
                {
                    std::cout << "count distinct" << std::endl;
                    throw -1;
                }
                else
                {
                    DEBUGcharIntervalTmpVec1.resize(count1);
                    DEBUGcharIntervalTmpVec2.resize(count2);
                    stool::beller::check(DEBUGcharIntervalTmpVec1, DEBUGcharIntervalTmpVec2);
                }
                return true;
            }
#endif
            uint64_t range_distinct(const RInterval<INDEX_SIZE> &range)
            {
                assert(this->lightDS != nullptr);
                uint64_t count = 0;
                if (range.endIndex - range.beginIndex <= range_distinct_threshold)
                {
                    count = this->lightDS->range_distinct(range.beginIndex, range.endIndex, this->charIntervalTmpVec);
                }
                else
                {
                    count = this->heavyDS->range_distinct(range.beginIndex, range.endIndex, this->charIntervalTmpVec);
                }
                //check(range);
                //uint64_t count = this->lightDS->range_distinct(range.beginIndex, range.endIndex, this->charIntervalTmpVec);

                assert(count > 0);

                for (uint64_t x = 0; x < count; x++)
                {
                    auto &it = this->charIntervalTmpVec[x];
                    INDEX_SIZE cBeginIndex = it.i;
                    INDEX_SIZE cEndIndex = it.j;
                    INDEX_SIZE cBeginDiff = cBeginIndex == range.beginIndex ? range.beginDiff : 0;
                    uint64_t end_run = (this->_RLBWTDS->lpos_vec)[cEndIndex + 1] - (this->_RLBWTDS->lpos_vec)[cEndIndex];
                    INDEX_SIZE cEndDiff = cEndIndex == range.endIndex ? range.endDiff : end_run - 1;

                    RInterval<INDEX_SIZE> cInterval;
                    cInterval.beginIndex = cBeginIndex;
                    cInterval.beginDiff = cBeginDiff;
                    cInterval.endIndex = cEndIndex;
                    cInterval.endDiff = cEndDiff;

                    charTmpVec[x] = it.c;
                    rIntervalTmpVec[x] = cInterval;
                }

                return count;
            }
            /*
            uint64_t computeNextLCPIntervalSet(const RINTERVAL &lcpIntv, const std::deque<RINTERVAL> &weinerVec, uint64_t rank, uint64_t width)
            {

                this->clear();

                this->computeSTNodeCandidates(lcpIntv);

                for (uint64_t j = 0; j < width; j++)
                {
                    this->computeSTChildren(weinerVec[rank + j]);
                }
                this->fit();

                this->check();

                return this->indexCount;
            }
            */
#if DEBUG
            /*
            void check2(const RINTERVAL &lcpIntv)
            {
                stool::LCPInterval<uint64_t> intv2 = lcpIntv.get_lcp_interval(this->_RLBWTDS->stnc->current_lcp - 1, this->_RLBWTDS->_fposDS);
                this->_RLBWTDS->checkWeinerLink(intv2, this->stnodeVec, this->indexVec, this->indexCount);
            }
            */
            void check2(const LCPInterval<uint64_t> &lcpIntv)
            {
                this->_RLBWTDS->checkWeinerLink(lcpIntv, this->stnodeVec, this->indexVec, this->indexCount);
            }
            void check3(uint64_t left, uint64_t right)
            {
                LCPInterval<uint64_t> lcpIntv;
                lcpIntv.i = left;
                lcpIntv.j = right;
                lcpIntv.lcp = this->_RLBWTDS->stnc->current_lcp - 1;
                this->_RLBWTDS->checkWeinerLink(lcpIntv, this->stnodeVec, this->indexVec, this->indexCount);
            }

#endif

            uint64_t computeFirstLCPIntervalSet()
            {
                this->clear();
                RINTERVAL lcpIntv;
                lcpIntv.beginIndex = this->_RLBWTDS->get_end_rle_lposition();
                lcpIntv.beginDiff = 0;
                lcpIntv.endIndex = this->_RLBWTDS->get_start_rle_lposition();
                lcpIntv.endDiff = this->_RLBWTDS->get_run(lcpIntv.endIndex) - 1;

                uint64_t begin_lindex = 0;
                uint64_t begin_diff = 0;
                uint64_t end_lindex = this->_RLBWTDS->rle_size() - 1;
                uint64_t end_diff = this->_RLBWTDS->get_run(end_lindex) - 1;

                RINTERVAL tmpArg;
                tmpArg.beginIndex = begin_lindex;
                tmpArg.beginDiff = begin_diff;
                tmpArg.endIndex = end_lindex;
                tmpArg.endDiff = end_diff;
                //std::vector<CHAR> charOutputVec;

                uint64_t resultCount = this->range_distinct(tmpArg);
                uint64_t counter = 0;
                this->pushLCPInterval(lcpIntv, 0);

                for (uint64_t x = 0; x < resultCount; x++)
                {
                    auto &it = this->rIntervalTmpVec[x];
                    assert(this->_RLBWTDS->get_char_by_run_index(it.beginIndex) == this->_RLBWTDS->get_char_by_run_index(it.endIndex));
                    //output.push_weiner(it);
                    this->pushExplicitWeinerInterval(it, 0);

                    counter++;
                }
                this->fit(true);

#if DEBUG
                stool::LCPInterval<uint64_t> lcpIntv2 = lcpIntv.get_lcp_interval(0, this->_RLBWTDS->_fposDS);
                this->check3(lcpIntv2.i, lcpIntv2.j);
#endif

                return 1;
                //return counter;
            }
            void fit(bool first)
            {
                uint64_t k = 0;

                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto c = this->indexVec[i];
                    uint64_t explicitChildrenCount = this->childrenVec[c].size();

#if DEBUG
                    if (this->_RLBWTDS->str_size() < 100)
                    {
                        std::cout << "FOUND ";
                        this->stnodeVec[c].print2(this->_RLBWTDS->_fposDS);
                        std::cout << std::endl;
                    }
#endif
                    assert(this->_RLBWTDS->checkLCPInterval(this->stnodeVec[c]) == (explicitChildrenCount > 0));
                    if (explicitChildrenCount > 0)
                    {

                        this->indexVec[k] = this->indexVec[i];

                        k++;
                    }
                }
                this->indexCount = k;

                if (!first)
                {

                    for (uint64_t i = 0; i < this->indexCount; i++)
                    {
                        auto c = this->indexVec[i];
                        uint64_t minCharIndex = 0;
                        uint64_t maxCharIndex = 0;

                        auto &currentVec = this->childrenVec[c];
                        uint64_t count = this->childrenVec[c].size();
                        for (uint64_t i = 0; i < count; i++)
                        {
                            auto &it = currentVec[i];
                            bool isLeft = it.beginIndex < currentVec[minCharIndex].beginIndex || (it.beginIndex == currentVec[minCharIndex].beginIndex && it.beginDiff < currentVec[minCharIndex].beginDiff);
                            if (isLeft)
                            {
                                minCharIndex = i;
                            }
                            bool isRight = it.endIndex > currentVec[maxCharIndex].endIndex || (it.endIndex == currentVec[maxCharIndex].endIndex && it.endDiff > currentVec[maxCharIndex].endDiff);
                            if (isRight)
                            {
                                maxCharIndex = i;
                            }
                        }
                        if (minCharIndex != 0)
                        {
                            auto tmp = currentVec[0];
                            currentVec[0] = currentVec[minCharIndex];
                            currentVec[minCharIndex] = tmp;
                        }

                        if (maxCharIndex == 0)
                            maxCharIndex = minCharIndex;
                        if (maxCharIndex != count - 1)
                        {
                            auto tmp2 = currentVec[count - 1];
                            currentVec[count - 1] = currentVec[maxCharIndex];
                            currentVec[maxCharIndex] = tmp2;
                        }
                    }
                }
                else
                {
                    for (uint64_t i = 0; i < this->indexCount; i++)
                    {
                        auto c = this->indexVec[i];
                        uint64_t minCharIndex = 0;
                        uint64_t maxCharIndex = 0;

                        auto &currentVec = this->childrenVec[c];
                        uint64_t count = this->childrenVec[c].size();
                        for (uint64_t i = 0; i < count; i++)
                        {
                            auto &it = currentVec[i];
                            bool isLeft = this->_RLBWTDS->bwt[it.beginIndex] < this->_RLBWTDS->bwt[currentVec[minCharIndex].beginIndex];
                            if (isLeft)
                            {
                                minCharIndex = i;
                            }
                            bool isRight = this->_RLBWTDS->bwt[it.beginIndex] > this->_RLBWTDS->bwt[currentVec[maxCharIndex].beginIndex];
                            if (isRight)
                            {
                                maxCharIndex = i;
                            }
                        }
                        if (minCharIndex != 0)
                        {
                            auto tmp = currentVec[0];
                            currentVec[0] = currentVec[minCharIndex];
                            currentVec[minCharIndex] = tmp;
                        }

                        if (maxCharIndex == 0)
                            maxCharIndex = minCharIndex;
                        if (maxCharIndex != count - 1)
                        {
                            auto tmp2 = currentVec[count - 1];
                            currentVec[count - 1] = currentVec[maxCharIndex];
                            currentVec[maxCharIndex] = tmp2;
                        }
                    }
                }
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool