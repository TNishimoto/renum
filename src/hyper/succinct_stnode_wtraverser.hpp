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
//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        template <typename INDEX_SIZE, typename RLBWTDS>
        class SuccinctSTNodeWTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            stool::EliasFanoVector children;
            sdsl::bit_vector leftmost_child_bits;
            sdsl::bit_vector::select_1_type leftmost_child_bits_selecter;

            bool builded = false;
            uint64_t _stnode_count = 0;
            uint64_t _children_count = 0;

        public:
            RLBWTDS *_RLBWTDS = nullptr;

            SuccinctSTNodeWTraverser()
            {
            }
            SuccinctSTNodeWTraverser(RLBWTDS *__RLBWTDS) : _RLBWTDS(__RLBWTDS)
            {
            }
            void create_next_wtraverser()
            {
            }
            uint64_t children_count()
            {
                return this->_children_count;
            }
            uint64_t node_count()
            {
                return this->_stnode_count;
            }

            uint64_t get_stnode_left_index(uint64_t i)
            {
                assert(i < this->_stnode_count);
                uint64_t L = this->leftmost_child_bits_selecter(i + 1);
                return L;
            }
            uint64_t get_stnode_right_index(uint64_t i)
            {
                assert(i < this->_stnode_count);

                uint64_t R = this->leftmost_child_bits_selecter(i + 2);
                return R - 1;
            }

            uint64_t get_child_node_left(uint64_t i)
            {
                return this->children[i * 2];
            }
            uint64_t get_child_node_right(uint64_t i)
            {
                return this->children[(i * 2) + 1];
            }
            uint64_t get_stnode2(uint64_t L, stool::LCPInterval<uint64_t> &output, uint64_t lcp)
            {
                assert(this->leftmost_child_bits[L]);
                uint64_t R = L + 1;
                while (!this->leftmost_child_bits[R])
                {
                    R++;
                }
                R--;
                output.i = this->get_child_node_left(L);
                output.j = this->get_child_node_right(R);
                output.lcp = lcp;

                return R + 1;
            }

            void build(stool::EliasFanoVectorBuilder &_children, sdsl::bit_vector &_leftmost_child_bits, uint64_t _stnode_count, uint64_t __children_count)
            {
                this->children.build_from_builder(_children);
                this->_children_count = __children_count;
                this->_stnode_count = _stnode_count;
                this->leftmost_child_bits.swap(_leftmost_child_bits);
                /*
                std::cout << "LEFTMOST BITS: ";
                for (uint64_t i = 0; i < this->leftmost_child_bits.size(); i++)
                {
                    std::cout << (this->leftmost_child_bits[i] ? "1" : "0");
                }
                std::cout << std::endl;
                */

                assert(this->leftmost_child_bits[0]);

                sdsl::bit_vector::select_1_type b_sel(&this->leftmost_child_bits);
                leftmost_child_bits_selecter.set_vector(&this->leftmost_child_bits);
                leftmost_child_bits_selecter.swap(b_sel);
                assert(leftmost_child_bits_selecter(1) == 0);

            }

            std::vector<std::pair<uint64_t, uint64_t>> to_plain()
            {
                std::vector<std::pair<uint64_t, uint64_t>> r;
                for (uint64_t i = 0; i < this->_children_count; i++)
                {
                    uint64_t left = this->children[i * 2];
                    uint64_t right = this->children[(i * 2) + 1];

                    r.push_back(std::pair<uint64_t, uint64_t>(left, right));
                }
                return r;
            }
            uint64_t get_using_memory() const
            {
                uint64_t x1 = sdsl::size_in_bytes(leftmost_child_bits);
                uint64_t x2 = this->children.get_using_memory();
                return x1 + x2;
            }
        };
        struct LightweightInterval
        {
            uint64_t left;
            uint64_t right;
            bool is_leftmost;

            LightweightInterval(){

            };
            LightweightInterval(uint64_t _left, uint64_t _right, bool _is_leftmost) : left(_left), right(_right), is_leftmost(_is_leftmost)
            {
            }
        };
        template <typename INDEX_SIZE, typename RLBWTDS>
        class SuccinctSTNodeWTraverserBuilder
        {
        public:
            using RINTERVAL = RInterval<INDEX_SIZE>;
            uint64_t builderRangeStart;
            uint64_t builderRangeEnd;
            RLBWTDS *_RLBWTDS = nullptr;
            std::vector<std::vector<LightweightInterval>> next_tmp;
            std::stack<uint8_t> foundCharacters;
            std::vector<bool> foundCharacterFlags;
            uint64_t next_child_count = 0;
            uint64_t next_node_count = 0;

            void initialize(uint64_t start, uint64_t end, RLBWTDS *__RLBWTDS)
            {
                this->builderRangeStart = start;
                this->builderRangeEnd = end;
                this->_RLBWTDS = __RLBWTDS;
                this->next_tmp.resize(256);
                this->foundCharacterFlags.resize(256, false);
            }
            void clear()
            {
                throw -1;
            }
            void buildSuccinctSTNodeWTraverser(SuccinctSTNodeWTraverser<INDEX_SIZE, RLBWTDS> &next)
            {
                next._RLBWTDS = this->_RLBWTDS;
                stool::EliasFanoVectorBuilder builder;
                builder.initialize(this->_RLBWTDS->str_size(), next_child_count * 2);

                sdsl::bit_vector leftmost_child_bits;
                leftmost_child_bits.resize(next_child_count + 1);

                uint64_t x = 0;
                for (auto &it : this->next_tmp)
                {
                    for (auto &it2 : it)
                    {
                        builder.push(it2.left);
                        builder.push(it2.right);

                        leftmost_child_bits[x++] = it2.is_leftmost;
                    }
                }
                leftmost_child_bits[x] = true;
                builder.finish();

                next.build(builder, leftmost_child_bits, this->next_node_count, this->next_child_count);
            }
            void push(LightweightInterval intv, char c)
            {
                this->next_tmp[c].push_back(intv);
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

            void computeNextSTNodes(SuccinctSTNodeWTraverser<INDEX_SIZE, RLBWTDS> &prev, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em, uint64_t st_index)
            {
                uint64_t st_left_index = prev.get_stnode_left_index(st_index);

                uint64_t st_right_index = prev.get_stnode_right_index(st_index);

                uint64_t leftmost_child_left = prev.get_child_node_left(st_left_index);

                uint64_t rightmost_child_right = prev.get_child_node_right(st_right_index);

                RINTERVAL w;
                this->_RLBWTDS->to_rinterval(leftmost_child_left, rightmost_child_right, w);

                em.clear();
                em.computeSTNodeCandidates2(w);

                //uint64_t x = 0;
                for (uint64_t i = st_left_index; i <= st_right_index; i++)
                {
                    uint64_t child_left = prev.get_child_node_left(i);
                    uint64_t child_right = prev.get_child_node_right(i);
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

                        next_tmp[c].push_back(LightweightInterval(w_left, w_right, x == 0));
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
            void computeNextSTNodes(SuccinctSTNodeWTraverser<INDEX_SIZE, RLBWTDS> &prev, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                for (uint64_t i = this->builderRangeStart; i <= this->builderRangeEnd; i++)
                {
                    this->computeNextSTNodes(prev, em, i);
                }
            }

            /*
            bool computeNextLCPIntervalSet(SuccinctSTNodeWTraverser<INDEX_SIZE, RLBWTDS> &prev, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                RINTERVAL intv;

                assert(this->childVec.size() + 1 == this->w_builder.size());
                uint64_t size = this->w_builder.size();
                uint64_t x = 1;
                for (uint64_t i = 1; i < size; i++)
                {
                    if (this->w_builder[x])
                    {
                        this->get_stnode(0, x - 1, intv);
                        em.clear();
                        em.computeSTNodeCandidates(intv);

                        for (uint64_t i = 0; i < x; i++)
                        {
                            auto &child = this->childVec[0];
                            em.computeSTChildren(child);
                            this->childVec.pop_front();
                            this->w_builder.pop_front();
                        }
                        em.fit();

                        //std::lock_guard<std::mutex> lock(test_mtx);
                        for (uint64_t i = 0; i < em.indexCount; i++)
                        {
                            auto c = em.indexVec[i];
                            auto &currentVec = em.childrenVec[c];
                            uint64_t count = currentVec.size();

                            bool limitOver = ((this->childVec.size() + count) > limit);
                            if (limitOver)
                            {
                                if (tmp.size() == 0 || tmp[tmp.size() - 1]->childVec.size() + count > limit)
                                {
                                    isSplit = true;
                                    tmp.push_back(new STNodeWTraverser(this->_RLBWTDS));
                                }
                            }

                            STNodeWTraverser *it = limitOver ? tmp[tmp.size() - 1] : this;

                            em.move_st_internal_nodes(it->childVec, it->w_builder, c);
                            assert(it->childVec.size() <= limit);
                            it->_stnode_count++;
                        }

                        this->_stnode_count--;

                        x = 1;
                    }
                    else
                    {
                        x++;
                    }
                }

                this->w_builder.pop_front();
                //this->_stnode_count = tmp_count;
                assert(this->childVec.size() == this->w_builder.size());

                return isSplit;
            }
            */
        };

    } // namespace lcp_on_rlbwt
} // namespace stool