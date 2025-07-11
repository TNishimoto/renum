#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>
#include "../basic/rinterval.hpp"
#include "./stnode_vector.hpp"
#include "../basic/interval_search_data_structure.hpp"
#include "../basic/sdsl_functions.hpp"

namespace stool
{
    namespace renum
    {

        /*
            This is a data structure to ...
        */
        template <typename INDEX_SIZE, typename CHAR = uint8_t>
        class ExplicitWeinerLinkComputer
        {
            using RINTERVAL = std::pair<INDEX_SIZE, INDEX_SIZE>;
            using CHAR_INTERVAL = CharInterval<INDEX_SIZE, CHAR>;
            std::vector<std::vector<CHAR_INTERVAL>> childrenVec;
            //std::vector<RINTERVAL> stnodeVec;
            std::vector<uint64_t> indexVec;
            //std::vector<std::vector<uint8_t>> edgeCharVec;
            //using CHAR = uint8_t;
            using INDEX = INDEX_SIZE;

            std::vector<bool> stnodeOccFlagArray;

            uint64_t indexCount = 0;
            //uint64_t explicitChildCount = 0;
            uint64_t range_distinct_threshold = 16;

            std::vector<uint8_t> charTmpVec;
            std::vector<RINTERVAL> rIntervalTmpVec;
            std::vector<CHAR_INTERVAL> charIntervalTmpVec;

            stool::renum::IntervalSearchDataStructure<CHAR> *searcher;
            //sdsl::bit_vector::rank_1_type *bwt_bit_rank1;
            uint64_t strSize = 0;

        public:
            //LightRangeDistinctDataStructure<typename RLBWTDS::CHAR_VEC, INDEX_SIZE> lightRangeSearcher;
            //SuccinctRangeDistinctDataStructure<INDEX_SIZE> heavyRangeSearcher;

            void initialize(IntervalSearchDataStructure<CHAR> *_searcher, uint64_t _strSize)
            {
                //this->bwt_bit_rank1 = _bwt_bit_rank1;
                uint64_t CHARMAX = UINT8_MAX + 1;
                childrenVec.resize(CHARMAX);
                indexVec.resize(CHARMAX);
                //edgeCharVec.resize(CHARMAX);

                stnodeOccFlagArray.resize(CHARMAX, false);
                //stnodeVec.resize(CHARMAX);

                rIntervalTmpVec.resize(CHARMAX);
                charTmpVec.resize(CHARMAX);

                charIntervalTmpVec.resize(CHARMAX);
                this->searcher = _searcher;
                this->strSize = _strSize;
            }
            void executeWeinerLinkSearch(std::pair<INDEX_SIZE, INDEX_SIZE> &node, std::vector<std::pair<INDEX_SIZE, INDEX_SIZE>> &children, std::vector<uint8_t> *edgeChars, stool::renum::STNodeVector<INDEX_SIZE, CHAR> &output_vec)
            {
                (void) node;
                bool _store_edge_chars = edgeChars != nullptr;

                this->clear();
                //this->computeSTNodeCandidates(node.first, node.second);
                if (_store_edge_chars)
                {
                    for (uint64_t i = 0; i < children.size(); i++)
                    {
                        this->computeSTChildren(children[i].first, children[i].second, (*edgeChars)[i]);
                    }
                }
                else
                {
                    for (auto &it : children)
                    {
                        this->computeSTChildren(it.first, it.second, 0);
                    }
                }
                this->fit();
                this->output(_store_edge_chars, output_vec);
            }
            void executeWeinerLinkSearch(stool::renum::STNodeVector<INDEX_SIZE, CHAR> &stack)
            {
                //uint64_t left = stack.get_last_left_boundary();
                //uint64_t right = stack.get_last_right_boundary();

                this->clear();
                //this->computeSTNodeCandidates(left, right);

                stool::renum::CharInterval<INDEX_SIZE, uint8_t> tmp;
                uint64_t L = stack.childs_vec.size() - stack.get_last_width();
                bool b = false;
                while (L < stack.childs_vec.size())
                {
                    L = stack.mini_increment(L, tmp, b);
                    this->computeSTChildren(tmp.i, tmp.j, tmp.c);
                }
                this->fit();

                bool _store_edge_chars = stack.edge_char_vec.size() > 0;
                int64_t lcp = stack.get_last_depth();
                stack.pop();

                this->sortedOutput(_store_edge_chars, lcp, stack);
            }

            std::vector<CharInterval<INDEX_SIZE, uint8_t>> getFirstChildren()
            {

                std::vector<CharInterval<INDEX_SIZE, uint8_t>> r;
                INDEX_SIZE left = 0;
                INDEX_SIZE right = this->strSize - 1;
                uint64_t count = this->searcher->getIntervals(left, right, this->charIntervalTmpVec);

                //r.resize(count - 1);

                for (uint64_t x = 0; x < count; x++)
                {
                    auto &it = this->charIntervalTmpVec[x];
                    r.push_back(CharInterval<INDEX_SIZE, uint8_t>(it.i, it.j, it.c));
                }
                sort(r.begin(), r.end(), [&](const CharInterval<INDEX_SIZE, uint8_t> &lhs, const CharInterval<INDEX_SIZE, uint8_t> &rhs) {
                    return lhs.c < rhs.c;
                });
                return r;
            }
            uint64_t get_input_text_length() const
            {
                return this->searcher->wt->size();
            }

        private:
            /*
            bool checkMaximalRepeat(uint64_t left, uint64_t right) const
            {

                uint64_t k1 = left == 0 ? 0 : (*bwt_bit_rank1)(left);
                uint64_t k2 = (*bwt_bit_rank1)(right + 1);
                bool isNotMaximalRepeat = ((k2 - k1 == 0) || ((k2 - k1) == (right - left + 1)));
            
                return !isNotMaximalRepeat;
            }
            */
            bool checkMaximalRepeat(uint64_t left, uint64_t right) const
            {
                uint8_t right_c = SDSLFunction::get_Char((*this->searcher->wt), right, this->searcher->lastChar);
                uint64_t right_c_num = SDSLFunction::get_rank((*this->searcher->wt), right + 1, right_c, this->searcher->lastChar);
                uint64_t left_c_num = SDSLFunction::get_rank((*this->searcher->wt), left, right_c, this->searcher->lastChar);

                bool isNotMaximalRepeat = ((right_c_num - left_c_num == 0) || ((right_c_num - left_c_num) == (right - left + 1)));
                return !isNotMaximalRepeat;

            }

            /*
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
            */

            void clear(CHAR c)
            {
                childrenVec[c].clear();
                stnodeOccFlagArray[c] = false;
            }
            void clear()
            {

                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto &it = this->indexVec[i];
                    childrenVec[it].clear();
                    //edgeCharVec[it].clear();

                    stnodeOccFlagArray[it] = false;
                }
                indexCount = 0;
                //explicitChildCount = 0;
            }
            uint8_t get_child(uint8_t c, uint64_t index, std::pair<INDEX_SIZE, INDEX_SIZE> &output) const
            {
                auto &it = this->childrenVec[c][index];
                output.first = it.i;
                output.second = it.j;
                return it.c;
            }
            /*
            uint64_t get_width(uint8_t c) const
            {
                auto &currentVec = this->childrenVec[c];
                return currentVec.size();
            }
            */

            void output(bool _store_edge_chars, stool::renum::STNodeVector<INDEX_SIZE, CHAR> &output_vec)
            {
                //RINTERVAL copy;

                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    CHAR c = this->indexVec[i];
                    this->output(c, _store_edge_chars, -1, output_vec);
                }
            }
            std::pair<INDEX_SIZE, INDEX_SIZE> get_node_ianterval(CHAR c)
            {
                auto &currentVec = this->childrenVec[c];
                uint64_t count = currentVec.size();
                uint64_t st_left = UINT64_MAX;
                uint64_t st_right = 0;
                std::pair<INDEX_SIZE, INDEX_SIZE> outputInterval;

                for (uint64_t j = 0; j < count; j++)
                {
                    this->get_child(c, j, outputInterval);
                    if (outputInterval.first < st_left)
                    {
                        st_left = outputInterval.first;
                    }
                    if (outputInterval.second > st_right)
                    {
                        st_right = outputInterval.second;
                    }
                }
                return std::pair<INDEX_SIZE, INDEX_SIZE>(st_left, st_right);
            }
            int64_t get_largest_next_interval_index()
            {
                int64_t largest_index = -1;
                int64_t largest_width = 0;

                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    CHAR c = this->indexVec[i];
                    std::pair<INDEX_SIZE, INDEX_SIZE> p = this->get_node_ianterval(c);
                    uint64_t w = p.second - p.first + 1;
                    if ((uint64_t)largest_width < w)
                    {
                        largest_index = c;
                        largest_width = w;
                    }
                }
                return largest_index;
            }
            void output(CHAR c, bool _store_edge_chars, int64_t parent_lcp, stool::renum::STNodeVector<INDEX_SIZE, CHAR> &output_vec)
            {
                auto &currentVec = this->childrenVec[c];
                uint64_t count = currentVec.size();

                uint64_t st_left = UINT64_MAX;
                uint64_t st_right = 0;
                std::pair<INDEX_SIZE, INDEX_SIZE> outputInterval;

                for (uint64_t j = 0; j < count; j++)
                {
                    CHAR edgeChar = this->get_child(c, j, outputInterval);
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
                if (parent_lcp != -1)
                {
                    output_vec.depth_vec.push_back(parent_lcp + 1);
                    output_vec.width_vec.push_back(count);
                }
            }
            void sortedOutput(bool _store_edge_chars, int64_t parent_lcp, stool::renum::STNodeVector<INDEX_SIZE, CHAR> &output_vec)
            {
                //RINTERVAL copy;
                if (this->indexCount > 0)
                {
                    int64_t largest_index = this->get_largest_next_interval_index();
                    this->output(largest_index, _store_edge_chars, parent_lcp, output_vec);

                    for (uint64_t i = 0; i < this->indexCount; i++)
                    {
                        int64_t c = this->indexVec[i];
                        if (c != largest_index)
                        {
                            this->output(c, _store_edge_chars, parent_lcp, output_vec);
                        }
                    }
                }
            }

        private:
            void pushExplicitWeinerInterval(INDEX_SIZE left, INDEX_SIZE right, uint8_t c, uint8_t edgeChar)
            {
                /*
                bool isValid = this->stnodeOccFlagArray[c];
                if (!isValid)
                    return false;
                auto &lcpIntv = this->stnodeVec[c];
                bool isLastChild = lcpIntv.second == right;
                bool isFirstChild = lcpIntv.first == left;
                bool b = !(isFirstChild && isLastChild);
                */
                bool firstOccurance = !this->stnodeOccFlagArray[c];
                if (firstOccurance)
                {
                    this->stnodeOccFlagArray[c] = true;
                    this->indexVec[this->indexCount++] = c;
                }
                this->childrenVec[c].push_back(CHAR_INTERVAL(left, right, edgeChar));
                //explicitChildCount++;

                //this->edgeCharVec[c].push_back(edgeChar);

                /*
                if (this->childrenVec[c].size() == 0)
                {

                    this->indexVec[this->indexCount] = c;
                    this->indexCount++;
                }
                */

                //return b;
            }
            /*
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
            */
            void computeSTChildren(INDEX_SIZE left, INDEX_SIZE right, uint8_t edgeChar)
            {

                assert(left <= right);
                uint64_t resultCount = this->searcher->getIntervals(left, right, this->charIntervalTmpVec);

                for (uint64_t i = 0; i < resultCount; i++)
                {
                    auto &it = this->charIntervalTmpVec[i];

                    this->pushExplicitWeinerInterval(it.i, it.j, it.c, edgeChar);
                }
            }

            void fit()
            {
                uint64_t k = 0;

                for (uint64_t i = 0; i < this->indexCount; i++)
                {
                    auto c = this->indexVec[i];
                    uint64_t explicitChildrenCount = this->childrenVec[c].size();

                    if (explicitChildrenCount >= 2)
                    {
                        this->indexVec[k] = this->indexVec[i];

                        k++;
                    }
                    else
                    {
                        this->clear(this->indexVec[i]);
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
                        bool isLeft = it.i < currentVec[minCharIndex].i;
                        if (isLeft)
                        {
                            minCharIndex = i;
                        }
                        bool isRight = it.j > currentVec[maxCharIndex].j;
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

        class WeinerLinkCommonFunctions
        {
        public:
            /*
            template <typename EM>
            static void output(const EM em, uint8_t c, bool _store_edge_chars, stool::renum::STNodeVector<typename EM::INDEX, typename EM::CHAR> &output_vec)
            {
                //RINTERVAL copy;

                auto &currentVec = em.childrenVec[c];
                uint64_t count = currentVec.size();

                uint64_t st_left = UINT64_MAX;
                uint64_t st_right = 0;
                std::pair<typename EM::INDEX, typename EM::INDEX> outputInterval;

                for (uint64_t j = 0; j < count; j++)
                {
                    typename EM::CHAR edgeChar = em.get_child(c, j, outputInterval);
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
                bool isMaximalRepeat = em.checkMaximalRepeat(st_left, st_right);

                output_vec.maximal_repeat_check_vec.push_back(isMaximalRepeat);
            }
            template <typename EM>
            static void output(const EM em, bool _store_edge_chars, stool::renum::STNodeVector<typename EM::INDEX, typename EM::CHAR> &output_vec)
            {
                for (uint64_t i = 0; i < em.indexCount; i++)
                {
                    typename EM::CHAR c = em.indexVec[i];
                    output(em, c, _store_edge_chars, output_vec);
                }
            }
            */
            template <typename EM, typename ITERATOR>
            static void compute_weiner_links(EM &em, const ITERATOR &it, stool::renum::STNodeVector<typename EM::INDEX, typename EM::CHAR> &output)
            {
                using CHAR = typename EM::CHAR;
                using INDEX = typename EM::INDEX;

                std::pair<INDEX, INDEX> node;
                node.first = it.get_left();
                node.second = it.get_right();

                std::vector<std::pair<INDEX, INDEX>> children;
                std::vector<CHAR> edgeChars;
                std::vector<CHAR> output_chars;

                uint64_t child_count = it.get_children_count();
                //uint64_t depth = it.get_depth();
                for (uint64_t i = 0; i < child_count; i++)
                {
                    uint64_t left = it.get_child_left_boundary(i);
                    uint64_t right = it.get_child_right_boundary(i);
                    children.push_back(std::pair<INDEX, INDEX>(left, right));
                    char c = it.get_edge_character(i);
                    edgeChars.push_back(c);
                }
                em.executeWeinerLinkSearch(node, children, it.has_edge_characters() ? &edgeChars : nullptr, output);
            }
        };

    } // namespace renum
} // namespace stool