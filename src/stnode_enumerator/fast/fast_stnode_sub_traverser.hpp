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
#include "../weiner_link_emulator.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {

        template <typename INDEX_SIZE, typename RLBWTDS>
        class FastSTNodeSubTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;

            uint64_t _stnode_count = 0;
            std::deque<INDEX_SIZE> childs_vec;
            std::deque<bool> first_child_flag_vec;
            std::deque<bool> maximal_repeat_check_vec;
            RLBWTDS *_RLBWTDS = nullptr;

        public:
            FastSTNodeSubTraverser()
            {
            }
            FastSTNodeSubTraverser(RLBWTDS *__RLBWTDS) : _RLBWTDS(__RLBWTDS)
            {
            }

        private:
            inline void get_child_start_position(uint64_t i, RINTERVAL &output) const
            {
                output.beginIndex = this->childs_vec[(i * 4)];
                output.beginDiff = this->childs_vec[(i * 4) + 1];
            }
            inline void get_child_end_position(uint64_t i, RINTERVAL &output) const
            {
                output.endIndex = this->childs_vec[(i * 4) + 2];
                output.endDiff = this->childs_vec[(i * 4) + 3];
            }

            void add(uint8_t c, uint64_t count, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                RINTERVAL copy;
                uint64_t st_left = UINT64_MAX;
                uint64_t st_right = 0;

                for (uint64_t j = 0; j < count; j++)
                {
                    em.get_child(c, j, copy);
                    uint64_t left = this->_RLBWTDS->get_fpos(copy.beginIndex, copy.beginDiff);
                    uint64_t right = this->_RLBWTDS->get_fpos(copy.endIndex, copy.endDiff);
                    this->_RLBWTDS->to_rinterval(left, right, copy);

                    if (left < st_left)
                    {
                        st_left = left;
                    }
                    if (right > st_right)
                    {
                        st_right = right;
                    }

                    this->childs_vec.push_back(copy.beginIndex);
                    this->childs_vec.push_back(copy.beginDiff);
                    this->childs_vec.push_back(copy.endIndex);
                    this->childs_vec.push_back(copy.endDiff);

                    this->first_child_flag_vec.push_back(j == 0);
                }
                uint64_t x = this->_RLBWTDS->get_lindex_containing_the_position(st_left);
                uint64_t d = this->_RLBWTDS->get_run(x);
                bool isMaximalRepeat = (this->_RLBWTDS->get_lpos(x) + d) <= st_right;
                this->maximal_repeat_check_vec.push_back(isMaximalRepeat);
                this->_stnode_count++;
            }

        public:
            bool check_maximal_repeat(uint64_t st_index) const
            {
                return this->maximal_repeat_check_vec[st_index];
            }
            uint64_t read_st_node(uint64_t L, RINTERVAL &output) const
            {
                assert(this->first_child_flag_vec[L]);
                assert(!this->first_child_flag_vec[L + 1]);

                uint64_t R = L + 1;
                while (R < this->first_child_flag_vec.size() && !this->first_child_flag_vec[R])
                {
                    R++;
                }

                R--;

                this->get_child_start_position(L, output);
                this->get_child_end_position(R, output);

                return R + 1;
            }
            uint64_t get_stnode(uint64_t L, stool::LCPInterval<uint64_t> &output, uint64_t lcp)
            {
                assert(this->first_child_flag_vec[L]);
                assert(!this->first_child_flag_vec[L + 1]);

                RINTERVAL intv;
                uint64_t newL = read_st_node(L, intv);

                output.i = this->_RLBWTDS->get_lpos(intv.beginIndex) + intv.beginDiff;
                output.j = this->_RLBWTDS->get_lpos(intv.endIndex) + intv.endDiff;
                output.lcp = lcp;

                return newL;
            }

            uint64_t children_count() const
            {
                return this->childs_vec.size() / 4;
            }
            uint64_t node_count() const
            {
                return this->_stnode_count;
            }

            void clear()
            {
                this->_stnode_count = 0;
                this->childs_vec.clear();
                this->maximal_repeat_check_vec.clear();
                this->first_child_flag_vec.clear();
            }
            void swap(FastSTNodeSubTraverser<INDEX_SIZE, RLBWTDS> &item)
            {
                this->childs_vec.swap(item.childs_vec);
                this->first_child_flag_vec.swap(item.first_child_flag_vec);

                this->maximal_repeat_check_vec.swap(item.maximal_repeat_check_vec);
                uint64_t tmp2 = this->_stnode_count;
                this->_stnode_count = item._stnode_count;
                item._stnode_count = tmp2;

                auto tmp3 = this->_RLBWTDS;
                this->_RLBWTDS = item._RLBWTDS;
                item._RLBWTDS = tmp3;
            }
            void first_compute(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                this->first_child_flag_vec.clear();

                //em.computeFirstLCPIntervalSet();
                auto r = em.getFirstChildren();
                RINTERVAL intv;
                for (uint64_t i = 0; i < r.size(); i++)
                {
                    auto &it = r[i];

                    this->_RLBWTDS->to_rinterval(it.first, it.second, intv);

                    this->childs_vec.push_back(intv.beginIndex);
                    this->childs_vec.push_back(intv.beginDiff);
                    this->childs_vec.push_back(intv.endIndex);
                    this->childs_vec.push_back(intv.endDiff);
                    this->first_child_flag_vec.push_back(i == 0);
                }
                this->maximal_repeat_check_vec.push_back(true);
                this->_stnode_count++;
            }
            bool computeNextLCPIntervalSet(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {

                bool isSplit = false;
                RINTERVAL intv;

                assert(this->children_count() == this->first_child_flag_vec.size());
                uint64_t size = this->_stnode_count;
                uint64_t L = 0;

                for (uint64_t i = 0; i < size; i++)
                {
                    assert(this->first_child_flag_vec[L]);

                    uint64_t nextL = this->read_st_node(L, intv);
                    em.clear();
#if DEBUG
                    if (this->_RLBWTDS->str_size() < 100)
                    {
                        std::cout << "Search input = ";
                        std::cout << std::endl;
                    }
#endif
                    em.computeSTNodeCandidates(intv);

                    RINTERVAL child;
                    for (uint64_t i = 0; i < nextL; i++)
                    {
                        this->get_child_start_position(0, child);
                        this->get_child_end_position(0, child);
                        em.computeSTChildren(child);
                        this->childs_vec.pop_front();
                        this->childs_vec.pop_front();
                        this->childs_vec.pop_front();
                        this->childs_vec.pop_front();

                        this->first_child_flag_vec.pop_front();
                    }

                    em.fit(false);
#if DEBUG
                    if (this->_RLBWTDS->stnc != nullptr)
                    {
                        uint64_t _left = this->_RLBWTDS->get_lpos(intv.beginIndex) + intv.beginDiff;
                        uint64_t _right = this->_RLBWTDS->get_lpos(intv.endIndex) + intv.endDiff;

                        em.verify_next_lcp_interval(_left, _right);
                    }
#endif
                    for (uint64_t i = 0; i < em.indexCount; i++)
                    {
                        auto c = em.indexVec[i];
                        auto &currentVec = em.childrenVec[c];
                        uint64_t count = currentVec.size();

                        //bool limitOver = false;
                        this->add(c, count, em);
                    }
                    this->maximal_repeat_check_vec.pop_front();
                    this->_stnode_count--;
                }

                assert(this->children_count() == this->first_child_flag_vec.size());

                return isSplit;
            }

            void print()
            {
                std::cout << "[" << this->node_count() << ", " << this->children_count() << "]" << std::endl;
            }
            void print_info()
            {
                uint64_t L = 0;
                RINTERVAL intv;
                for (uint64_t i = 0; i < this->node_count(); i++)
                {
                    L = read_st_node(L, intv);
                    intv.print();
                }

                std::cout << std::endl;

                for (uint64_t i = 0; i < this->childs_vec.size(); i += 2)
                {
                    assert(i + 1 < this->childs_vec.size());
                    std::cout << "[" << this->childs_vec[i] << ", " << this->childs_vec[i + 1] << "]";
                }

                std::cout << std::endl;
            }
            uint64_t get_using_memory()
            {
                uint64_t x1 = this->childs_vec.size() * sizeof(INDEX_SIZE);
                uint64_t x2 = (this->first_child_flag_vec.size() * 1) / 8;
                //uint64_t x3 = (this->leftmost_child_bits.size() * 1) / 8;
                uint64_t x4 = (maximal_repeat_check_vec.size() * 1) / 8;
                return x1 + x2 + x4;
            }

#if DEBUG
            void bit_check()
            {
                uint64_t k = 0;
                for (uint64_t i = 0; i < first_child_flag_vec.size(); i++)
                {
                    if (first_child_flag_vec[i])
                    {
                        k++;
                    }
                }
                assert(this->_stnode_count == k);
            }
#endif
            void merge(FastSTNodeSubTraverser<INDEX_SIZE, RLBWTDS> &item)
            {
                while (item.childs_vec.size() > 0)
                {
                    this->childs_vec.push_back(item.childs_vec[0]);
                    item.childs_vec.pop_front();
                }
                while (item.first_child_flag_vec.size() > 0)
                {
                    this->first_child_flag_vec.push_back(item.first_child_flag_vec[0]);
                    item.first_child_flag_vec.pop_front();
                }
                while (item.maximal_repeat_check_vec.size() > 0)
                {
                    this->maximal_repeat_check_vec.push_back(item.maximal_repeat_check_vec[0]);
                    item.maximal_repeat_check_vec.pop_front();
                }
                this->_stnode_count += item._stnode_count;
                item.clear();
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool