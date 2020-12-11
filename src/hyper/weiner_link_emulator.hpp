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
        /*
            This is a data structure to ...
        */
        template <typename INDEX_SIZE, typename RLBWTDS>
        class ExplicitWeinerLinkEmulator
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;

            std::vector<std::vector<RINTERVAL>> childrenVec;
            std::vector<RINTERVAL> stnodeVec;
            std::vector<uint64_t> indexVec;
            std::vector<bool> stnodeOccFlagArray;

            uint64_t indexCount = 0;
            uint64_t explicitChildCount = 0;

        public:
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

            void move_st_internal_nodes(std::vector<RINTERVAL> &outputSTVec, std::vector<RINTERVAL> &outputExplicitChildrenVec, std::vector<uint8_t> &outputWidthVec)
            {
                for (uint64_t i = 0; i < this->indexCount; i++)
                {

                    auto character = this->indexVec[i];


                    auto &currentVec = this->childrenVec[character];
                    uint64_t count = this->childrenVec[character].size();
                    //uint64_t count = this->explicitChildrenCountVec[character];
                    outputSTVec.push_back(this->stnodeVec[character]);
                    outputWidthVec.push_back(count);
                    for (uint64_t j = 0; j < count; j++)
                    {
                        outputExplicitChildrenVec.push_back(currentVec[j]);
                        //output.push_weiner(currentVec[j]);
                    }
                }
            }

        private:
            void pushExplicitWeinerInterval(const RINTERVAL &w, uint8_t c)
            {

                if (this->childrenVec[c].size() == 0)
                {

                    this->indexVec[this->indexCount] = c;
                    this->indexCount++;
                }
                this->childrenVec[c].push_back(w);

                explicitChildCount++;
            }
            void pushLCPInterval(const RINTERVAL &w, uint8_t c)
            {
                this->stnodeVec[c] = w;
                this->stnodeOccFlagArray[c] = true;
            }

            //template <typename LPOSDS, typename RANGEDS>
            void getNextSAIntervalsForLCPIntervals(const RINTERVAL &w)
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

        public:
            void computeNextExplicitWeinerIntervals(const RINTERVAL &w)
            {
                RINTERVAL frontL = this->_RLBWTDS->getIntervalOnL(w);

                assert(frontL.beginIndex <= frontL.endIndex);
                uint64_t resultCount = this->range_distinct(frontL);

                for (uint64_t i = 0; i < resultCount; i++)
                {
                    typename RLBWTDS::UCHAR c = this->charTmpVec[i];
                    auto &it = this->rIntervalTmpVec[i];

                    auto &lcpIntv = this->stnodeVec[c];
                    bool isLastChild = lcpIntv.endIndex == it.endIndex && lcpIntv.endDiff == it.endDiff;
                    if (!isLastChild)
                    {
                        this->pushExplicitWeinerInterval(it, c);
                    }
                }
            }
            uint64_t thshold = 16;
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
            uint64_t range_distinct(const RInterval<INDEX_SIZE> &range)
            {
                assert(this->lightDS != nullptr);
                uint64_t count = 0;
                if (range.endIndex - range.beginIndex > thshold)
                {
                    count = this->heavyDS->range_distinct(range.beginIndex, range.endIndex, this->charIntervalTmpVec);
                }
                else
                {
                    count = this->lightDS->range_distinct(range.beginIndex, range.endIndex, this->charIntervalTmpVec);
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
            void computeNextLCPIntervalSet(const RINTERVAL &lcpIntv, const std::vector<RINTERVAL> &weinerVec, uint64_t weinerVecSize, uint64_t rank)
            {
                this->clear();
                uint64_t i = rank;

                INDEX_SIZE lcpIntvBeginPos = this->_RLBWTDS->get_fpos(lcpIntv.beginIndex, lcpIntv.beginDiff);
                INDEX_SIZE lcpIntvEndPos = this->_RLBWTDS->get_fpos(lcpIntv.endIndex, lcpIntv.endDiff);

                this->getNextSAIntervalsForLCPIntervals(lcpIntv);

                while (i < weinerVecSize)
                {

                    INDEX_SIZE weinerBeginPos = this->_RLBWTDS->get_fpos(weinerVec[i].beginIndex, weinerVec[i].beginDiff);
                    INDEX_SIZE weinerEndPos = this->_RLBWTDS->get_fpos(weinerVec[i].endIndex, weinerVec[i].endDiff);
                    if (lcpIntvBeginPos <= weinerBeginPos && weinerEndPos <= lcpIntvEndPos)
                    {
                        this->computeNextExplicitWeinerIntervals(weinerVec[i]);
                    }
                    else
                    {

                        break;
                    }

                    i++;
                }
                this->fit();

                //intervalTemporary.move(outputSet);
            }
            void computeFirstLCPIntervalSet()
            {
                this->clear();
                RINTERVAL lcpIntv;
                size_t minIndex = this->_RLBWTDS->get_end_rle_lposition();

                uint64_t maxIndex = 0;
                typename RLBWTDS::UCHAR c = 0;
                for (uint64_t i = 0; i < this->_RLBWTDS->rle_size(); i++)
                {
                    typename RLBWTDS::UCHAR c2 = this->_RLBWTDS->get_char_by_run_index(i);
                    if (c < c2)
                    {
                        c = c2;
                        maxIndex = i;
                    }
                    else if (c == c2)
                    {
                        maxIndex = i;
                    }
                }
                uint64_t diff = this->_RLBWTDS->get_run(maxIndex) - 1;
                lcpIntv.beginIndex = minIndex;
                lcpIntv.beginDiff = 0;
                lcpIntv.endIndex = maxIndex;
                lcpIntv.endDiff = diff;

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

                //return counter;
            }
            void fit()
            {
                uint64_t k = 0;
                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto character = this->indexVec[i];
                    uint64_t explicitChildrenCount = this->childrenVec[character].size();
                    if (explicitChildrenCount > 0)
                    {
                        this->indexVec[k++] = this->indexVec[i];
                    }
                }
                this->indexCount = k;
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool