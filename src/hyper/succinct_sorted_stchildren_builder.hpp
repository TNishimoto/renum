#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <deque>
#include <vector>
#include <type_traits>
#include "weiner_link_emulator.hpp"
#include "stool/src/elias_fano_vector.hpp"
#include "succinct_sorted_stchildren.hpp"
//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        template <typename INDEX_SIZE, typename RLBWTDS>
        class SuccinctSortedSTChildrenBuilder
        {
        public:
            using RINTERVAL = RInterval<INDEX_SIZE>;
            uint64_t builderRangeStart;
            uint64_t builderRangeEnd;
            RLBWTDS *_RLBWTDS = nullptr;
            SuccinctSortedSTChildren<INDEX_SIZE, RLBWTDS> *prevSortedSTChildren;

            std::vector<stool::EliasFanoVectorBuilder> next_tmp;
            std::vector<std::vector<bool>> next_tmp_bits;
            std::vector<uint64_t> occurrence_counter;
            std::vector<uint64_t> minimum_counter;
            std::vector<uint64_t> maximum_counter;

            std::stack<uint8_t> foundCharacters;
            std::vector<bool> foundCharacterFlags;

            std::vector<bool> rangeDistinctFlags;

            uint64_t next_child_count = 0;
            uint64_t next_node_count = 0;
            uint64_t nextSTChildrenCount = 0;

            uint64_t debug = 0;

            void initialize(uint64_t start, uint64_t end, RLBWTDS *__RLBWTDS, SuccinctSortedSTChildren<INDEX_SIZE, RLBWTDS> *prev)
            {
                this->builderRangeStart = start;
                this->builderRangeEnd = end;
                this->_RLBWTDS = __RLBWTDS;
                this->next_tmp.resize(256);
                this->next_tmp_bits.resize(256);
                this->foundCharacterFlags.resize(256, false);
                this->occurrence_counter.resize(256, 0);
                this->minimum_counter.resize(256, UINT64_MAX);
                this->maximum_counter.resize(256, 0);

                prevSortedSTChildren = prev;
            }
            void set(uint64_t universe, uint64_t total_element_count)
            {
                for (uint64_t i = 0; i < this->next_tmp.size(); i++)
                {
                    this->occurrence_counter[i] = 0;
                    this->next_tmp[i].initialize(universe, total_element_count * 2);
                }
            }
            void set()
            {

                for (uint64_t i = 0; i < this->next_tmp.size(); i++)
                {
                    if (this->occurrence_counter[i] > 0)
                    {
                        uint64_t universe = this->maximum_counter[i] - this->minimum_counter[i] + 1;
                        this->next_tmp[i].initialize(universe, this->occurrence_counter[i] * 2);
                    }
                }
            }

            void clear()
            {
                throw -1;
            }
            static void merge(SuccinctSortedSTChildren<INDEX_SIZE, RLBWTDS> &next, std::vector<SuccinctSortedSTChildrenBuilder *> &elements)
            {
                uint64_t child_count_sum = 0;
                uint64_t node_count_sum = 0;
                for (auto &it : elements)
                {
                    child_count_sum += it->next_child_count;
                    node_count_sum += it->next_node_count;
                }

                next._RLBWTDS = elements[0]->_RLBWTDS;
                stool::EliasFanoVectorBuilder builder;
                builder.initialize(next._RLBWTDS->str_size(), child_count_sum * 2);
                assert(builder.universe <= UINT32_MAX);
                sdsl::bit_vector leftmost_child_bits;
                leftmost_child_bits.resize(child_count_sum + 1);

                uint64_t x = 0;
                for (uint64_t c = 0; c < 256; c++)
                {
                    for (auto &it : elements)
                    {

                        if (it->occurrence_counter[c] > 0)
                        {
                            builder.merge(it->next_tmp[c], it->minimum_counter[c]);

                            for (uint64_t j = 0; j < it->next_tmp_bits[c].size(); j++)
                            {
                                assert(x < leftmost_child_bits.size());

                                leftmost_child_bits[x++] = it->next_tmp_bits[c][j];

                            }
                        }
                    }
                }
                assert(x + 1 ==  leftmost_child_bits.size());
                leftmost_child_bits[x] = true;
                for (uint64_t x = 0; x < leftmost_child_bits.size(); x++)
                {
                    std::cout << (leftmost_child_bits[x] ? "1" : "0");
                }
                std::cout << std::endl;
                builder.finish();
                next.build(builder, leftmost_child_bits, node_count_sum, child_count_sum);
            }
            /*
            void buildSuccinctSortedSTChildren(SuccinctSortedSTChildren<INDEX_SIZE, RLBWTDS> &next)
            {
                next._RLBWTDS = this->_RLBWTDS;
                stool::EliasFanoVectorBuilder builder;
                builder.initialize(this->_RLBWTDS->str_size(), next_child_count * 2);

                sdsl::bit_vector leftmost_child_bits;
                leftmost_child_bits.resize(next_child_count + 1);

                uint64_t x = 0;

                for (uint64_t c = 0; c < this->next_tmp_bits.size(); c++)
                {
                    if (this->next_tmp_bits[c].size() > 0)
                    {
                        std::vector<uint64_t> ptmp;
                        this->next_tmp[c].to_vector(ptmp);
                        for (uint64_t i = 0; i < this->next_tmp_bits[c].size(); i++)
                        {
                            builder.push(ptmp[i * 2]);
                            builder.push(ptmp[(i * 2) + 1]);
                            leftmost_child_bits[x++] = this->next_tmp_bits[c][i];
                        }
                    }
                }
                leftmost_child_bits[x] = true;
                builder.finish();

                next.build(builder, leftmost_child_bits, this->next_node_count, this->next_child_count);
            }
            */
            void push(LightweightInterval intv, char c)
            {
                next_tmp[c].push(intv.left);
                next_tmp[c].push(intv.right);
                next_tmp_bits[c].push_back(intv.is_leftmost);

                this->occurrence_counter[c] += 1;
                this->minimum_counter[c] = 0;

                //this->next_tmp[c].push_back(intv);
                this->next_child_count++;
                if (!this->foundCharacterFlags[c])
                {
                    this->foundCharacterFlags[c] = true;
                    this->foundCharacters.push(c);
                }
                assert(this->next_node_count > 0 || intv.is_leftmost);
                if (intv.is_leftmost)
                {
                    next_node_count++;
                }
            }

            void computeNextSTNodes(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em, uint64_t st_index)
            {
                uint64_t st_left_index = prevSortedSTChildren->get_stnode_left_index(st_index);
                uint64_t st_right_index = prevSortedSTChildren->get_stnode_right_index(st_index);
                uint64_t leftmost_child_left = prevSortedSTChildren->get_child_node_left(st_left_index);
                uint64_t rightmost_child_right = prevSortedSTChildren->get_child_node_right(st_right_index);

                RINTERVAL w;
                this->_RLBWTDS->to_rinterval(leftmost_child_left, rightmost_child_right, w);

                em.clear();
                em.computeSTNodeCandidates2(w);

                //uint64_t x = 0;
                for (uint64_t i = st_left_index; i <= st_right_index; i++)
                {
                    uint64_t child_left = prevSortedSTChildren->get_child_node_left(i);
                    uint64_t child_right = prevSortedSTChildren->get_child_node_right(i);
                    this->_RLBWTDS->to_rinterval(child_left, child_right, w);
                    em.computeSTChildren2(w);
                }
                em.fit();
#if DEBUG
                LCPInterval<uint64_t> test;
                test.i = leftmost_child_left;
                test.j = rightmost_child_right;
                test.lcp = this->_RLBWTDS->stnc->current_lcp;
                em.check2(test);
#endif
                /*
                assert(this->debug >= em.indexCount);
                std::cout << "DEB = " << this->debug << "/" << em.indexCount << std::endl;
                this->debug -= em.indexCount;
                */

                this->next_node_count += em.indexCount;
                for (uint64_t i = 0; i < em.indexCount; i++)
                {
                    auto c = em.indexVec[i];
                    auto &currentVec = em.childrenVec[c];
                    uint64_t count = currentVec.size();

                    if (!this->foundCharacterFlags[c])
                    {
                        this->foundCharacterFlags[c] = true;
                        this->foundCharacters.push(c);
                    }

                    for (uint64_t x = 0; x < count; x++)
                    {
                        RINTERVAL &w1 = currentVec[x];
                        uint64_t w_left = this->_RLBWTDS->get_fpos(w1.beginIndex, w1.beginDiff);
                        uint64_t w_right = this->_RLBWTDS->get_fpos(w1.endIndex, w1.endDiff);
#if DEBUG
                        if (this->_RLBWTDS->str_size() < 100)
                        {
                            std::cout << "FOUND: [" << w_left << ", " << w_right << "]";
                        }
#endif
                        next_tmp[c].push(w_left - this->minimum_counter[c]);
                        next_tmp[c].push(w_right - this->minimum_counter[c]);
                        next_tmp_bits[c].push_back(x == 0);

                        //next_tmp[c].push_back(LightweightInterval(w_left, w_right, x == 0));
                        next_child_count++;
                    }
#if DEBUG
                    if (this->_RLBWTDS->str_size() < 100)
                    {
                        std::cout << std::endl;
                    }

#endif
                }
            }
            uint64_t countNextSTNodes(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em, uint64_t st_index)
            {
                uint64_t st_left_index = prevSortedSTChildren->get_stnode_left_index(st_index);
                uint64_t st_right_index = prevSortedSTChildren->get_stnode_right_index(st_index);
                uint64_t leftmost_child_left = prevSortedSTChildren->get_child_node_left(st_left_index);
                uint64_t rightmost_child_right = prevSortedSTChildren->get_child_node_right(st_right_index);

                RINTERVAL w;
                this->_RLBWTDS->to_rinterval(leftmost_child_left, rightmost_child_right, w);

                em.clear();
                em.computeSTNodeCandidates2(w);
                uint64_t k = 0;
                for (uint64_t i = st_left_index; i <= st_right_index; i++)
                {
                    uint64_t child_left = prevSortedSTChildren->get_child_node_left(i);
                    uint64_t child_right = prevSortedSTChildren->get_child_node_right(i);
                    this->_RLBWTDS->to_rinterval(child_left, child_right, w);
                    em.computeSTChildren2(w);
                }
                em.fit();

                for (uint64_t i = 0; i < em.indexCount; i++)
                {
                    auto c = em.indexVec[i];
                    auto &currentVec = em.childrenVec[c];
                    uint64_t count = currentVec.size();
                    this->occurrence_counter[c] += count;
                    k += count;

                    for (uint64_t x = 0; x < count; x++)
                    {
                        RINTERVAL &w1 = currentVec[x];
                        uint64_t w_left = this->_RLBWTDS->get_fpos(w1.beginIndex, w1.beginDiff);
                        uint64_t w_right = this->_RLBWTDS->get_fpos(w1.endIndex, w1.endDiff);
                        if (w_left < this->minimum_counter[c])
                        {
                            this->minimum_counter[c] = w_left;
                        }
                        if (w_right > this->maximum_counter[c])
                        {
                            this->maximum_counter[c] = w_right;
                        }
                    }
                }

                return k;
            }
            uint64_t countNextSTNodes(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                assert(this->next_node_count == 0);
                assert(this->next_child_count == 0);

                uint64_t k = 0;
                for (uint64_t i = this->builderRangeStart; i <= this->builderRangeEnd; i++)
                {
                    k += this->countNextSTNodes(em, i);
                }
                return k;
            }

            void computeNextSTNodes(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                assert(this->next_node_count == 0);
                assert(this->next_child_count == 0);

                for (uint64_t i = this->builderRangeStart; i <= this->builderRangeEnd; i++)
                {
                    this->computeNextSTNodes(em, i);
                }
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool