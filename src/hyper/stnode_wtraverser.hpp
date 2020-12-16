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
//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        template <typename INDEX_SIZE, typename RLBWTDS>
        class STNodeWTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;

            uint64_t _stnode_count = 0;
            std::deque<uint64_t> childs_vec;
            std::deque<bool> first_child_flag_vec;
            std::vector<bool> maximal_repeat_check_vec;
            RLBWTDS *_RLBWTDS = nullptr;

        public:

            RINTERVAL get_child(uint64_t i){
                RINTERVAL intv;
                intv.beginIndex = this->childs_vec[(i*4)];
                intv.beginDiff = this->childs_vec[(i*4)+1];
                intv.endIndex = this->childs_vec[(i*4)+2];
                intv.endDiff = this->childs_vec[(i*4)+3];
                return intv;
            }
            


            STNodeWTraverser()
            {
            }
            STNodeWTraverser(RLBWTDS *__RLBWTDS) : _RLBWTDS(__RLBWTDS)
            {
            }

            /*
            void get_stnode(uint64_t i, RINTERVAL &output)
            {
                assert(builded);
                assert(i <= this->_stnode_count);
                uint64_t L = this->leftmost_child_bits_selecter(i + 1);
                uint64_t R = this->leftmost_child_bits_selecter(i + 2) - 1;

                this->get_stnode(L, R, output);
            }
            */
            
            void get_stnode(uint64_t L, uint64_t R, RINTERVAL &output)
            {
                RINTERVAL intv1 = this->get_child(L);
                RINTERVAL intv2 = this->get_child(R);

                uint64_t beg_index = intv1.beginIndex, end_index = intv2.endIndex,
                         beg_diff = intv1.beginDiff, end_diff = intv2.endDiff;

                output.beginIndex = beg_index;
                output.beginDiff = beg_diff;
                output.endIndex = end_index;
                output.endDiff = end_diff;
                if (this->_RLBWTDS->bwt[beg_index] != this->_RLBWTDS->bwt[end_index])
                {
                    output.beginIndex = this->_RLBWTDS->get_end_rle_lposition();
                    output.beginDiff = 0;
                    output.endIndex = this->_RLBWTDS->get_start_rle_lposition();
                    output.endDiff = this->_RLBWTDS->get_run(output.endIndex) - 1;
                }
            }
            
            uint64_t get_stnode2(uint64_t L, RINTERVAL &output)
            {
                assert(this->first_child_flag_vec[L]);
                assert(!this->first_child_flag_vec[L+1]);

                uint64_t R = L + 1;
                while (!this->first_child_flag_vec[R] )
                {
                    R++;
                }
                R--;

                
                RINTERVAL intv1 = this->get_child(L);
                RINTERVAL intv2 = this->get_child(R);
                uint64_t beg_index = intv1.beginIndex, end_index = intv2.endIndex,
                         beg_diff = intv1.beginDiff, end_diff = intv2.endDiff;

                output.beginIndex = beg_index;
                output.beginDiff = beg_diff;
                output.endIndex = end_index;
                output.endDiff = end_diff;
                if (this->_RLBWTDS->bwt[beg_index] != this->_RLBWTDS->bwt[end_index])
                {
                    output.beginIndex = this->_RLBWTDS->get_end_rle_lposition();
                    output.beginDiff = 0;
                    output.endIndex = this->_RLBWTDS->get_start_rle_lposition();
                    output.endDiff = this->_RLBWTDS->get_run(output.endIndex) - 1;
                }
                return R + 1;
            }
            uint64_t get_stnode2(uint64_t L, stool::LCPInterval<uint64_t> &output, uint64_t lcp)
            {
                assert(this->first_child_flag_vec[L]);
                assert(!this->first_child_flag_vec[L+1]);

                RINTERVAL tmp;
                uint64_t newL = get_stnode2(L, tmp);
                uint64_t left = this->_RLBWTDS->get_fpos(tmp.beginIndex, tmp.beginDiff);
                uint64_t right = this->_RLBWTDS->get_fpos(tmp.endIndex, tmp.endDiff);
                output.i = left;
                output.j = right;
                output.lcp = lcp;
                return newL;

            }

            /*
            std::vector<RINTERVAL> *get_stnode_vec()
            {
                return &this->stnodeVec;
            }
            */

            std::vector<RINTERVAL> *get_child_vec()
            {
                return &this->child_vec;
            }

            /*
            uint64_t size(){
                return this->stnodeVec.size();
            }
            */
            uint64_t children_count() const
            {
                return this->childs_vec.size() / 4;
            }
            uint64_t node_count() const
            {
                return this->_stnode_count;
            }
            /*
            uint64_t get_width(uint64_t i) const
            {
                if (builded)
                {
                    uint64_t q = this->leftmost_child_bits_selecter(i + 2) - this->leftmost_child_bits_selecter(i + 1);
                    return q;
                }
                else
                {
                    std::cout << "Not supported get_width" << std::endl;
                    throw -1;
                }
            }
            */

            void clear()
            {
                this->_stnode_count = 0;
                this->childs_vec.clear();
                this->maximal_repeat_check_vec.clear();
                //this->first_child_flag_vec.clear();

                //this->leftmost_child_bits.clear();
            }
            void swap(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &item)
            {
                this->childs_vec.swap(item.childs_vec);
                this->first_child_flag_vec.swap(item.first_child_flag_vec);
                //bool tmp = this->builded;
                //this->builded = item.builded;
                //item.builded = tmp;

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
                first_child_flag_vec.clear();

                this->_stnode_count += em.computeFirstLCPIntervalSet();
                em.move_st_internal_nodes(this->childs_vec, first_child_flag_vec);
                /*
#if DEBUG                
                this->_RLBWTDS->checkLCPInterval(this->stnodeVec[0]);
#endif
*/

                this->maximal_repeat_check_vec.resize(1);
                this->maximal_repeat_check_vec[0] = true;

            }
            /*
            uint64_t get_child_rank(uint64_t i)
            {
                return this->leftmost_child_bits_selecter(i + 1);
            }
            */
            bool computeNextLCPIntervalSet(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em, std::vector<STNodeWTraverser *> &tmp, uint64_t limit)
            {
                //this->first_child_flag_vec.clear();
                bool isSplit = false;
                RINTERVAL intv;
                //uint64_t tmp_count = 0;

                assert(this->children_count() + 1 == this->first_child_flag_vec.size());
                uint64_t size = this->first_child_flag_vec.size();
                uint64_t x = 1;
                for (uint64_t i = 1; i < size; i++)
                {
                    if (this->first_child_flag_vec[x])
                    {
                        this->get_stnode(0, x - 1, intv);
                        em.clear();
                        #if DEBUG
                            if(this->_RLBWTDS->str_size() < 100){
                                std::cout << "Search input = ";
                                intv.print2(this->_RLBWTDS->_fposDS);
                                std::cout << std::endl;

                            }
                        #endif
                        em.computeSTNodeCandidates(intv);

                        for (uint64_t i = 0; i < x; i++)
                        {
                            RINTERVAL child = this->get_child(0);
                            //auto &child = this->child_vec[0];
                            em.computeSTChildren(child);
                            this->childs_vec.pop_front();
                            this->childs_vec.pop_front();
                            this->childs_vec.pop_front();
                            this->childs_vec.pop_front();

                            this->first_child_flag_vec.pop_front();
                        }
                        em.fit();
#if DEBUG
                        em.check2(intv);
#endif

                        //std::lock_guard<std::mutex> lock(test_mtx);
                        for (uint64_t i = 0; i < em.indexCount; i++)
                        {
                            auto c = em.indexVec[i];
                            auto &currentVec = em.childrenVec[c];
                            uint64_t count = currentVec.size();

                            bool limitOver = ((this->children_count() + count) > limit);
                            if (limitOver)
                            {
                                if (tmp.size() == 0 || tmp[tmp.size() - 1]->children_count() + count > limit)
                                {
                                    isSplit = true;
                                    tmp.push_back(new STNodeWTraverser(this->_RLBWTDS));
                                }
                            }

                            STNodeWTraverser *it = limitOver ? tmp[tmp.size() - 1] : this;

                            em.move_st_internal_nodes(it->childs_vec, it->first_child_flag_vec, c);

                            assert(it->children_count() <= limit);
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

                this->first_child_flag_vec.pop_front();
                //this->_stnode_count = tmp_count;
                assert(this->children_count() == this->first_child_flag_vec.size());

                return isSplit;
            }
            
            void add_the_last_bit_into_bit_array()
            {
                first_child_flag_vec.push_back(true);
            }
            

            void print()
            {
                std::cout << "[" << this->node_count() << ", " << this->children_count() << "]" << std::endl;
            }
            void print_info()
            {
                RINTERVAL it;
                it.beginIndex = 0;
                it.beginDiff = 0;
                it.endIndex = 0;
                it.endDiff = 0;

                uint64_t L = 0;
                for (uint64_t i = 0; i < this->node_count(); i++)
                {

                    L = this->get_stnode2(L, it);
                    it.print2(this->_RLBWTDS->_fposDS);
                }
                std::cout << std::endl;
            }
            uint64_t get_peak_memory()
            {
                uint64_t x1 = this->child_vec.size() * sizeof(uint64_t);
                uint64_t x2 = (this->first_child_flag_vec.size() * 1) / 8;
                //uint64_t x3 = (this->leftmost_child_bits.size() * 1) / 8;
                uint64_t x4 = (maximal_repeat_check_vec.capacity() * 1) / 8;

                return x1 + x2 + x4;
            }
            uint64_t get_optimal_memory()
            {
                uint64_t x1 = this->child_vec.size() * sizeof(uint64_t);
                uint64_t x2 = (this->first_child_flag_vec.size() * 1) / 8;
                //uint64_t x3 = (this->leftmost_child_bits.size() * 1) / 8;
                uint64_t x4 = (maximal_repeat_check_vec.size() * 1) / 8;
                return x1 + x2 + x4;
            }
            /*
            void split(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &item)
            {
                uint64_t k = this->child_vec.size() / 2;

                while (!this->first_child_flag_vec[k])
                {
                    k--;
                }
                uint64_t num = 0;
                for (uint64_t i = k; i < this->child_vec.size(); i++)
                {
                    if (this->first_child_flag_vec[i])
                    {
                        num++;
                    }
                    item.child_vec.push_back(this->child_vec[i]);
                    item.first_child_flag_vec.push_back(this->first_child_flag_vec[i]);
                }
                item._stnode_count += num;
                this->_stnode_count -= num;
                this->child_vec.resize(k);
                this->first_child_flag_vec.resize(k);

                this->child_vec.shrink_to_fit();
                this->first_child_flag_vec.shrink_to_fit();
                this->bit_check();
            }
            */
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
            void merge(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &item)
            {
                while (item.childs_vec.size() > 0)
                {
                    this->childs_vec.push_back(item.childs_vec[0]);
                    item.childs_vec.pop_front();
                }
                while(item.first_child_flag_vec.size() > 0){
                    this->first_child_flag_vec.push_back(item.first_child_flag_vec[0]);
                    item.first_child_flag_vec.pop_front();
                }
                this->_stnode_count += item._stnode_count;
                item.clear();
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool