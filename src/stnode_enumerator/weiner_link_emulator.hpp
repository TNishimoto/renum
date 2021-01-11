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
#include "../debug/stnode_checker.hpp"
#include "../rlbwt/rle.hpp"
#include "./stnode_vector.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        /*
            This is a data structure to ...
        */
        template <typename RLBWTDS>
        class ExplicitWeinerLinkEmulator
        {

            stool::lcp_on_rlbwt::STNodeChecker *stnc;
        public:
            using CHAR = typename RLBWTDS::CHAR;
            using INDEX = typename RLBWTDS::INDEX;
            using RINTERVAL = RInterval<INDEX>;

            std::vector<std::vector<RINTERVAL>> childrenVec;
            std::vector<std::vector<uint8_t>> edgeCharVec;

            std::vector<RINTERVAL> stnodeVec;
            std::vector<uint64_t> indexVec;
            std::vector<bool> stnodeOccFlagArray;

            bool checker_on = false;

            uint64_t indexCount = 0;
            uint64_t explicitChildCount = 0;
            uint64_t range_distinct_threshold = 16;
            CHAR special_character = 0;

            // For range distinct
            std::vector<uint8_t> charTmpVec;
            vector<RINTERVAL> rIntervalTmpVec;
            std::vector<CharInterval<INDEX>> charIntervalTmpVec;

            RLBWTDS *_RLBWTDS;
            stool::lcp_on_rlbwt::RLE<CHAR> *rlbwt;
            LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX> lightRangeSearcher;
            SuccinctRangeDistinctDataStructure<INDEX> heavyRangeSearcher;
            uint64_t get_input_text_length() const
            {
                return this->rlbwt->str_size();
            }


            void initialize(RLBWTDS *_rlbwtds)
            {
                _RLBWTDS = _rlbwtds;
                this->rlbwt = _RLBWTDS->get_rlbwt();
                this->special_character = rlbwt->get_smallest_character();
                uint64_t CHARMAX = UINT8_MAX + 1;
                childrenVec.resize(CHARMAX);
                edgeCharVec.resize(CHARMAX);
                indexVec.resize(CHARMAX);

                stnodeOccFlagArray.resize(CHARMAX, false);
                stnodeVec.resize(CHARMAX);

                rIntervalTmpVec.resize(CHARMAX);
                charTmpVec.resize(CHARMAX);

                charIntervalTmpVec.resize(CHARMAX);

                uint8_t lastChar = rlbwt->get_char_by_run_index(rlbwt->rle_size() - 1);
                //rlbwt->bwt[rlbwt->bwt.size() - 1];

                auto head_chars = _RLBWTDS->get_head_chars_pointer();
                auto wt  = _RLBWTDS->get_wavelet_tree_pointer();

                lightRangeSearcher.preprocess(head_chars);
                heavyRangeSearcher.initialize(wt, lastChar);
            }
            uint64_t get_explicit_stnode_count() const 
            {
                return this->indexCount;
            }
            uint64_t get_explicit_children_count() const
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
                    edgeCharVec[it].clear();

                    stnodeOccFlagArray[it] = false;
                }
                indexCount = 0;
                explicitChildCount = 0;
            }
            uint8_t get_child(uint8_t c, uint64_t index, RINTERVAL &output) const
            {
                auto &it = this->childrenVec[c][index];
                uint64_t left = this->_RLBWTDS->get_fpos(it.beginIndex, it.beginDiff);
                uint64_t right = this->_RLBWTDS->get_fpos(it.endIndex, it.endDiff);
                this->_RLBWTDS->to_rinterval(left, right, output);
                return this->edgeCharVec[c][index];
            }
            uint8_t get_child(uint8_t c, uint64_t index, std::pair<INDEX, INDEX> &output) const
            {
                auto &it = this->childrenVec[c][index];
                uint64_t left = this->_RLBWTDS->get_fpos(it.beginIndex, it.beginDiff);
                uint64_t right = this->_RLBWTDS->get_fpos(it.endIndex, it.endDiff);
                output.first = left;
                output.second = right;
                return this->edgeCharVec[c][index];
            }
            bool checkMaximalRepeat(uint64_t left, uint64_t right) const
            {
                return this->_RLBWTDS->checkMaximalRepeat(left, right);
            }

            void output(uint8_t c, bool _store_edge_chars, stool::lcp_on_rlbwt::STNodeVector<INDEX, CHAR> &output_vec)
            {
                //RINTERVAL copy;
                
                auto &currentVec = this->childrenVec[c];
                uint64_t count = currentVec.size();

                uint64_t st_left = UINT64_MAX;
                uint64_t st_right = 0;
                std::pair<INDEX, INDEX> outputInterval;

                for (uint64_t j = 0; j < count; j++)
                {
                    CHAR edgeChar = this->get_child(c, j, outputInterval);
                    //uint64_t left = ds->get_fpos(copy.beginIndex, copy.beginDiff);
                    //uint64_t right = ds->get_fpos(copy.endIndex, copy.endDiff);

                    if (outputInterval.first < st_left)
                    {
                        st_left = outputInterval.first;
                    }
                    if (outputInterval.second > st_right)
                    {
                        st_right = outputInterval.second;
                    }
                    if (j == 0)
                    {
                        output_vec.childs_vec.push_back(outputInterval.first);
                        output_vec.first_child_flag_vec.push_back(true);
                        if (_store_edge_chars)
                        {
                            output_vec.edge_char_vec.push_back(0);
                        }
                    }
                    output_vec.childs_vec.push_back(outputInterval.second);
                    output_vec.first_child_flag_vec.push_back(false);
                    if (_store_edge_chars)
                    {
                        output_vec.edge_char_vec.push_back(edgeChar);
                    }
                }
                bool isMaximalRepeat = this->checkMaximalRepeat(st_left, st_right);

                output_vec.maximal_repeat_check_vec.push_back(isMaximalRepeat);
            }
            void executeWeinerLinkSearch(std::pair<INDEX, INDEX> &node,
                                         std::vector<std::pair<INDEX, INDEX>> &children, std::vector<uint8_t> *edgeChars, std::vector<uint8_t> &output_chars)
            {
                this->clear();
                RINTERVAL intv;
                this->_RLBWTDS->to_rinterval(node.first, node.second, intv);

                this->computeSTNodeCandidates(intv);
                RINTERVAL child;
                if (edgeChars == nullptr)
                {
                    for (auto &it : children)
                    {
                        this->_RLBWTDS->to_rinterval(it.first, it.second, child);
                        this->computeSTChildren(child, 0);
                    }
                }
                else
                {
                    for (uint64_t i = 0; i < children.size(); i++)
                    {
                        this->_RLBWTDS->to_rinterval(children[i].first, children[i].second, child);
                        this->computeSTChildren(child, (*edgeChars)[i]);
                    }
                }
                this->fit();
#if DEBUG
                if (this->stnc != nullptr)
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
            void executeWeinerLinkSearch(const RINTERVAL &node, std::vector<RINTERVAL> &children, std::vector<uint8_t> *edgeChars, std::vector<uint8_t> &output_chars)
            {

                this->clear();

                this->computeSTNodeCandidates(node);
                if (edgeChars == nullptr)
                {
                    for (auto &it : children)
                    {
                        this->computeSTChildren(it, 0);

                    }
                }
                else
                {
                    for (uint64_t i = 0; i < children.size(); i++)
                    {
                        this->computeSTChildren(children[i], (*edgeChars)[i]);
                    }
                }
                this->fit();

#if DEBUG
                if (this->stnc != nullptr && this->checker_on)
                {
                    uint64_t left = this->rlbwt->get_lpos(node.beginIndex) + node.beginDiff;
                    uint64_t right = this->rlbwt->get_lpos(node.endIndex) + node.endDiff;

                    this->verify_next_lcp_interval(left, right);
                }
#endif

                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto c = this->indexVec[i];
                    output_chars.push_back(c);
                }
            }

            std::vector<CharInterval<INDEX>> getFirstChildren()
            {
                std::vector<CharInterval<INDEX>> r;
                uint64_t count = this->heavyRangeSearcher.range_distinct(0, this->rlbwt->rle_size() - 1, this->charIntervalTmpVec);
                //r.resize(count - 1);

                for (uint64_t x = 0; x < count; x++)
                {
                    auto &it = this->charIntervalTmpVec[x];
                    //uint8_t c = it.c;
                    uint64_t run = this->rlbwt->get_run(it.j) - 1;
                    uint64_t i = this->_RLBWTDS->get_fpos(it.i, 0);
                    uint64_t j = this->_RLBWTDS->get_fpos(it.j, run);
                    r.push_back(CharInterval<INDEX>(i, j, it.c));
                }
                sort(r.begin(), r.end(), [&](const CharInterval<INDEX> &lhs, const CharInterval<INDEX> &rhs) {
                    return lhs.c < rhs.c;
                });
                return r;
            }
            uint64_t get_width(uint8_t c)
            {
                auto &currentVec = this->childrenVec[c];
                return currentVec.size();
            }

            bool pushExplicitWeinerInterval(const RINTERVAL &w, uint8_t c, uint8_t edgeChar)
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
                    this->edgeCharVec[c].push_back(edgeChar);
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
                    if (c != this->special_character)
                    {
                        this->stnodeVec[c] = it;
                        this->stnodeOccFlagArray[c] = true;
                    }
                }
            }
            void computeSTChildren(const RINTERVAL &w, uint8_t edgeChar)
            {

                assert(w.beginIndex <= w.endIndex);
                uint64_t resultCount = this->range_distinct(w);

                for (uint64_t i = 0; i < resultCount; i++)
                {
                    typename RLBWTDS::UCHAR c = this->charTmpVec[i];
                    auto &it = this->rIntervalTmpVec[i];
                    this->pushExplicitWeinerInterval(it, c, edgeChar);
                }
            }

#if DEBUG
            bool check(const RInterval<INDEX> &range)
            {
                uint64_t CHARMAX = UINT8_MAX + 1;

                std::vector<CharInterval<INDEX>> DEBUGcharIntervalTmpVec1;
                std::vector<CharInterval<INDEX>> DEBUGcharIntervalTmpVec2;
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
            uint64_t range_distinct(const RInterval<INDEX> &range)
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
                    INDEX cBeginIndex = it.i;
                    INDEX cEndIndex = it.j;
                    INDEX cBeginDiff = cBeginIndex == range.beginIndex ? range.beginDiff : 0;
                    uint64_t end_run = this->rlbwt->get_run(cEndIndex);
                    INDEX cEndDiff = cEndIndex == range.endIndex ? range.endDiff : end_run - 1;

                    RInterval<INDEX> cInterval;
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
            std::vector<stool::LCPInterval<uint64_t>> createNextWeinerLinkNodes(uint64_t lcp)
            {
                std::vector<stool::LCPInterval<uint64_t>> wlinks;
                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    stool::LCPInterval<uint64_t> intv = this->stnodeVec[this->indexVec[i]].get_lcp_interval(lcp, this->_RLBWTDS->_fposDS);
                    wlinks.push_back(intv);
                }
                return wlinks;
            }
            void verify_next_lcp_interval(uint64_t left, uint64_t right)
            {
                LCPInterval<uint64_t> lcpIntv;
                lcpIntv.i = left;
                lcpIntv.j = right;
                assert(this->stnc != nullptr);
                lcpIntv.lcp = this->stnc->get_lcp() - 1;

                if (this->stnc != nullptr)
                {
                    uint64_t lcp = this->stnc->get_lcp();
                    auto wlinks = createNextWeinerLinkNodes(lcp);

                    this->stnc->check_weiner_links(left, right, wlinks);
                }
            }
#endif

            void fit()
            {
                uint64_t k = 0;

                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto c = this->indexVec[i];
                    uint64_t explicitChildrenCount = this->childrenVec[c].size();

#if DEBUG
                    if (this->rlbwt->str_size() < 100)
                    {
                        std::cout << "FOUND ";
                        this->stnodeVec[c].print2(this->_RLBWTDS->_fposDS);
                        std::cout << std::endl;
                    }
                    /*
                    if (this->stnc != nullptr && this->checker_on)
                    {
                        assert(this->_RLBWTDS->checkLCPInterval(this->stnodeVec[c]) == (explicitChildrenCount > 0));
                    }
                    */
#endif
                    if (explicitChildrenCount > 0)
                    {

                        this->indexVec[k] = this->indexVec[i];

                        k++;
                    }
                }
                this->indexCount = k;

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
        };

    } // namespace lcp_on_rlbwt
} // namespace stool