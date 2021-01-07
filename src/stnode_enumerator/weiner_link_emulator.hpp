#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>
#include "../rlbwt/range_distinct/light_range_distinct.hpp"

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

            bool checker_on = false;

            uint64_t indexCount = 0;
            uint64_t explicitChildCount = 0;
            uint64_t range_distinct_threshold = 16;

            // For range distinct
            std::vector<uint8_t> charTmpVec;
            vector<RINTERVAL> rIntervalTmpVec;
            std::vector<CharInterval<INDEX_SIZE>> charIntervalTmpVec;

            RLBWTDS *_RLBWTDS;
            //LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE> *lightDS = nullptr;
            //SuccinctRangeDistinctDataStructure<INDEX_SIZE> *heavyDS = nullptr;

            LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE> lightRangeSearcher;
            SuccinctRangeDistinctDataStructure<INDEX_SIZE> heavyRangeSearcher;
            uint64_t get_input_text_length()
            {
                return this->_RLBWTDS->str_size();
            }

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

                uint8_t lastChar = _RLBWTDS->bwt[_RLBWTDS->bwt.size() - 1];

                lightRangeSearcher.preprocess(&_RLBWTDS->bwt);
                heavyRangeSearcher.initialize(&_RLBWTDS->wt, lastChar);
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
            void get_child(uint8_t c, uint64_t index, RINTERVAL &output)
            {
                auto &it = this->childrenVec[c][index];
                uint64_t left = this->_RLBWTDS->get_fpos(it.beginIndex, it.beginDiff);
                uint64_t right = this->_RLBWTDS->get_fpos(it.endIndex, it.endDiff);
                this->_RLBWTDS->to_rinterval(left, right, output);
            }
            void get_child(uint8_t c, uint64_t index, std::pair<INDEX_SIZE, INDEX_SIZE> &output)
            {
                auto &it = this->childrenVec[c][index];
                uint64_t left = this->_RLBWTDS->get_fpos(it.beginIndex, it.beginDiff);
                uint64_t right = this->_RLBWTDS->get_fpos(it.endIndex, it.endDiff);
                output.first = left;
                output.second = right;
            }
            bool checkMaximalRepeat(uint64_t left, uint64_t right)
            {
                return this->_RLBWTDS->checkMaximalRepeat(left, right);
            }

            void executeWeinerLinkSearch(std::pair<INDEX_SIZE, INDEX_SIZE> &node, std::vector<std::pair<INDEX_SIZE, INDEX_SIZE>> &children, std::vector<uint8_t> &output_chars)
            {
                this->clear();
                this->computeSTNodeCandidates(node.first, node.second);
                for (auto &it : children)
                {
                    this->computeSTChildren(it.first, it.second);
                }
                this->fit(false);
#if DEBUG
                if (this->_RLBWTDS->stnc != nullptr)
                {
                    this->verify_next_lcp_interval(node.first, node.second);
                }
#endif
                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto c = this->indexVec[i];
                    output_chars.push_back(c);
                }
            }
            void executeWeinerLinkSearch(const RINTERVAL &node, std::vector<RINTERVAL> &children, std::vector<uint8_t> &output_chars)
            {
                this->clear();
                this->computeSTNodeCandidates(node);
                for (auto &it : children)
                {
                    this->computeSTChildren(it);
                }
                this->fit(false);

#if DEBUG
                if (this->_RLBWTDS->stnc != nullptr && this->checker_on)
                {
                    uint64_t left = this->_RLBWTDS->get_lpos(node.beginIndex) + node.beginDiff;
                    uint64_t right = this->_RLBWTDS->get_lpos(node.endIndex) + node.endDiff;

                    this->verify_next_lcp_interval(left, right);
                }
#endif

                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto c = this->indexVec[i];
                    output_chars.push_back(c);
                }
            }

            std::vector<std::pair<uint64_t, uint64_t>> getFirstChildren()
            {
                std::vector<std::pair<uint64_t, uint64_t>> r;
                uint64_t count = this->heavyRangeSearcher.range_distinct(0, this->_RLBWTDS->rle_size() - 1, this->charIntervalTmpVec);
                //r.resize(count - 1);

                for (uint64_t x = 0; x < count; x++)
                {
                    auto &it = this->charIntervalTmpVec[x];
                    //uint8_t c = it.c;
                    uint64_t run = this->_RLBWTDS->get_run(it.j) - 1;
                    uint64_t i = this->_RLBWTDS->get_fpos(it.i, 0);
                    uint64_t j = this->_RLBWTDS->get_fpos(it.j, run);
                    r.push_back(std::pair<uint64_t, uint64_t>(i, j));
                }
                sort(r.begin(), r.end(), [&](const std::pair<uint64_t, uint64_t> &lhs, const std::pair<uint64_t, uint64_t> &rhs) {
                    return lhs.first < rhs.first;
                });
                return r;
            }
            uint64_t get_width(uint8_t c)
            {
                auto &currentVec = this->childrenVec[c];
                return currentVec.size();
            }

        private:

            bool pushExplicitWeinerInterval(const RINTERVAL &w, uint8_t c)
            {
                bool isValid = this->stnodeOccFlagArray[c];
                if (!isValid)
                    return false;
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
            void computeSTNodeCandidates(const RINTERVAL &w)
            {
                uint64_t resultCount = this->range_distinct(w);
                for (uint64_t i = 0; i < resultCount; i++)
                {
                    typename RLBWTDS::UCHAR c = this->charTmpVec[i];
                    auto &it = this->rIntervalTmpVec[i];

                    //this->pushLCPInterval(it, c);
                    if (c != 0)
                    {
                        this->stnodeVec[c] = it;
                        this->stnodeOccFlagArray[c] = true;
                    }
                }
            }

            void computeSTNodeCandidates(INDEX_SIZE left, INDEX_SIZE right)
            {
                RINTERVAL intv;
                this->_RLBWTDS->to_rinterval(left, right, intv);
                this->computeSTNodeCandidates(intv);
            }
            void computeSTChildren(const RINTERVAL &w)
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
            void computeSTChildren(INDEX_SIZE left, INDEX_SIZE right)
            {
                RINTERVAL child;
                this->_RLBWTDS->to_rinterval(left, right, child);
                this->computeSTChildren(child);
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
                uint64_t count1 = heavyRangeSearcher.range_distinct(range.beginIndex, range.endIndex, DEBUGcharIntervalTmpVec1);
                uint64_t count2 = lightRangeSearcher.range_distinct(range.beginIndex, range.endIndex, DEBUGcharIntervalTmpVec2);
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
                uint64_t count = 0;
                if (range.endIndex - range.beginIndex <= range_distinct_threshold)
                {
                    count = this->lightRangeSearcher.range_distinct(range.beginIndex, range.endIndex, this->charIntervalTmpVec);
                }
                else
                {
                    count = this->heavyRangeSearcher.range_distinct(range.beginIndex, range.endIndex, this->charIntervalTmpVec);
                }

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
#if DEBUG

            void verify_next_lcp_interval(uint64_t left, uint64_t right)
            {
                LCPInterval<uint64_t> lcpIntv;
                lcpIntv.i = left;
                lcpIntv.j = right;
                assert(this->_RLBWTDS->stnc != nullptr);
                lcpIntv.lcp = this->_RLBWTDS->stnc->get_lcp() - 1;
                this->_RLBWTDS->checkWeinerLink(lcpIntv, this->stnodeVec, this->indexVec, this->indexCount);
            }

#endif

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
                    if (this->_RLBWTDS->stnc != nullptr && this->checker_on)
                    {
                        assert(this->_RLBWTDS->checkLCPInterval(this->stnodeVec[c]) == (explicitChildrenCount > 0));
                    }
#endif
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