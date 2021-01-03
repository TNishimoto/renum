#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>
#include "../rlbwt/rinterval.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        /*
            This is a data structure to ...
        */
        template <typename INDEX_SIZE>
        class ExplicitWeinerLinkSearch
        {
            using RINTERVAL = std::pair<INDEX_SIZE, INDEX_SIZE>;

        public:
            std::vector<std::vector<RINTERVAL>> childrenVec;
            std::vector<RINTERVAL> stnodeVec;
            std::vector<uint64_t> indexVec;
            std::vector<bool> stnodeOccFlagArray;

            uint64_t indexCount = 0;
            uint64_t explicitChildCount = 0;
            uint64_t range_distinct_threshold = 16;

            std::vector<uint8_t> charTmpVec;
            vector<RINTERVAL> rIntervalTmpVec;
            std::vector<CharInterval<INDEX_SIZE>> charIntervalTmpVec;

            stool::IntervalSearchDataStructure *searcher;
            sdsl::bit_vector::rank_1_type *bwt_bit_rank1;
            uint64_t strSize = 0;
            //LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE> lightRangeSearcher;
            //SuccinctRangeDistinctDataStructure<INDEX_SIZE> heavyRangeSearcher;

            void initialize(IntervalSearchDataStructure *_searcher, sdsl::bit_vector::rank_1_type *_bwt_bit_rank1, uint64_t _strSize)
            {
                this->bwt_bit_rank1 = _bwt_bit_rank1;
                uint64_t CHARMAX = UINT8_MAX + 1;
                childrenVec.resize(CHARMAX);
                indexVec.resize(CHARMAX);

                stnodeOccFlagArray.resize(CHARMAX, false);
                stnodeVec.resize(CHARMAX);

                rIntervalTmpVec.resize(CHARMAX);
                charTmpVec.resize(CHARMAX);

                charIntervalTmpVec.resize(CHARMAX);
                this->searcher = _searcher;
                this->strSize = _strSize;
            }
            bool checkMaximalRepeat(uint64_t left, uint64_t right)
            {

                uint64_t k1 = left == 0 ? 0 : (*bwt_bit_rank1)(left);
                uint64_t k2 = (*bwt_bit_rank1)(right + 1);
                bool b = !((k2 - k1 == 0) || ((k2 - k1) == (right - left + 1)));
                return b;
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
            void get_child(uint8_t c, uint64_t index, std::pair<INDEX_SIZE, INDEX_SIZE> &output)
            {
                auto &it = this->childrenVec[c][index];
                output.first = it.first;
                output.second = it.second;
            }
            uint64_t get_width(uint8_t c)
            {
                auto &currentVec = this->childrenVec[c];
                return currentVec.size();
            }

            bool pushExplicitWeinerInterval(INDEX_SIZE left, INDEX_SIZE right, uint8_t c)
            {
                bool isValid = this->stnodeOccFlagArray[c];
                if (!isValid)
                    return false;
                auto &lcpIntv = this->stnodeVec[c];
                bool isLastChild = lcpIntv.second == right;
                bool isFirstChild = lcpIntv.first == left;
                bool b = !(isFirstChild && isLastChild);

                if (b)
                {

                    if (this->childrenVec[c].size() == 0)
                    {

                        this->indexVec[this->indexCount] = c;
                        this->indexCount++;
                    }
                    this->childrenVec[c].push_back(RINTERVAL(left, right));
                    explicitChildCount++;
                }
                return b;
            }
            void computeSTNodeCandidates(INDEX_SIZE left, INDEX_SIZE right)
            {

                uint64_t resultCount = this->searcher->getIntervals(left, right, this->charIntervalTmpVec);
                for (uint64_t i = 0; i < resultCount; i++)
                {
                    auto &it = this->charIntervalTmpVec[i];

                    if (it.c != 0)
                    {
                        this->stnodeVec[it.c] = RINTERVAL(it.i, it.j);
                        this->stnodeOccFlagArray[it.c] = true;
                    }
                }
            }
            void executeWeinerLinkSearch(std::pair<INDEX_SIZE, INDEX_SIZE> &node, std::vector<std::pair<INDEX_SIZE, INDEX_SIZE>> &children, std::vector<uint8_t> &output_chars)
            {
                this->clear();
                this->computeSTNodeCandidates(node.first, node.second);
                for (auto &it : children)
                {
                    this->computeSTChildren(it.first, it.second);
                }
                this->fit();

                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto c = this->indexVec[i];
                    output_chars.push_back(c);
                }
            }

            void computeSTChildren(INDEX_SIZE left, INDEX_SIZE right)
            {

                assert(left <= right);
                uint64_t resultCount = this->searcher->getIntervals(left, right, this->charIntervalTmpVec);

                for (uint64_t i = 0; i < resultCount; i++)
                {
                    auto &it = this->charIntervalTmpVec[i];

                    this->pushExplicitWeinerInterval(it.i, it.j, it.c);
                }
            }

            std::vector<RINTERVAL> getFirstChildren()
            {

                std::vector<RINTERVAL> r;
                INDEX_SIZE left = 0;
                INDEX_SIZE right = this->strSize - 1;
                uint64_t count = this->searcher->getIntervals(left, right, this->charIntervalTmpVec);

                //r.resize(count - 1);

                for (uint64_t x = 0; x < count; x++)
                {
                    auto &it = this->charIntervalTmpVec[x];
                    r.push_back(RINTERVAL(it.i, it.j));
                }
                sort(r.begin(), r.end(), [&](const RINTERVAL &lhs, const RINTERVAL &rhs) {
                    return lhs.first < rhs.first;
                });
                return r;
            }
            void fit()
            {
                uint64_t k = 0;

                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto c = this->indexVec[i];
                    uint64_t explicitChildrenCount = this->childrenVec[c].size();

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
                        bool isLeft = it.first < currentVec[minCharIndex].first;
                        if (isLeft)
                        {
                            minCharIndex = i;
                        }
                        bool isRight = it.second > currentVec[maxCharIndex].second;
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