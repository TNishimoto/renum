#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
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
            //std::vector<RINTERVAL> stnodeVec;
            std::vector<RINTERVAL> childVec;
            //std::vector<uint8_t> _widthVec;

            sdsl::bit_vector leftmost_child_bits;
            sdsl::bit_vector::select_1_type leftmost_child_bits_selecter;
            bool builded = false;
            uint64_t _stnode_count = 0;

        public:
            std::vector<bool> w_builder;

            RLBWTDS *_RLBWTDS;

            std::vector<bool> maximal_repeat_check_vec;

            STNodeWTraverser()
            {
            }
            void get_stnode(uint64_t i, RINTERVAL &output)
            {
                assert(builded);
                uint64_t L = this->leftmost_child_bits_selecter(i + 1);
                uint64_t R = this->leftmost_child_bits_selecter(i + 2) - 1;

                this->get_stnode(L, R, output);
            }
            void get_stnode(uint64_t L, uint64_t R, RINTERVAL &output)
            {
                uint64_t beg_index = childVec[L].beginIndex, end_index = childVec[R].endIndex,
                         beg_diff = childVec[L].beginDiff, end_diff = childVec[R].endDiff;
                /*
                output.beginIndex = beg_index;
                output.beginDiff = beg_diff;
                output.endIndex = end_index;
                output.endDiff = end_diff;

                for (uint64_t x = L + 1; x <= R; x++)
                {
                    auto &it = childVec[x];

                    bool isLeft = it.beginIndex < beg_index || (it.beginIndex == beg_index && it.beginDiff < beg_diff);
                    assert(!isLeft);
                    if (isLeft)
                    {
                        beg_index = it.beginIndex;
                        beg_diff = it.beginDiff;
                    }

                    bool isRight = it.endIndex > end_index || (it.endIndex == end_index && it.endDiff > end_diff);
                    assert(!isRight);

                    if (isRight)
                    {
                        end_index = it.endIndex;
                        end_diff = it.endDiff;
                    }
                }
                */
                
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

                //return this->stnodeVec[i];
            }
            /*
            std::vector<RINTERVAL> *get_stnode_vec()
            {
                return &this->stnodeVec;
            }
            */

            std::vector<RINTERVAL> *get_child_vec()
            {
                return &this->childVec;
            }

            /*
            uint64_t size(){
                return this->stnodeVec.size();
            }
            */
            uint64_t children_count() const
            {
                return this->childVec.size();
            }
            uint64_t node_count() const
            {
                return this->_stnode_count;
            }

            uint64_t get_width(uint64_t i) const
            {
                //uint64_t p = this->_widthVec[i];
                if (builded)
                {
                    uint64_t q = this->leftmost_child_bits_selecter(i + 2) - this->leftmost_child_bits_selecter(i + 1);
                    return q;
                    //assert(p == q);
                }
                else
                {
                    std::cout << "Not supported get_width" << std::endl;
                    throw -1;
                }
            }
            void clear()
            {
                this->_stnode_count = 0;
                //this->stnodeVec.clear();
                this->childVec.clear();
                //this->_widthVec.clear();
                this->maximal_repeat_check_vec.clear();
                this->w_builder.clear();
                this->builded = false;

                //this->leftmost_child_bits.clear();
            }
            void shrink(uint64_t node_capacity, uint64_t child_capacity){
                if(this->childVec.capacity() > child_capacity){
                    this->childVec.reserve(child_capacity);
                    this->childVec.shrink_to_fit();
                    this->w_builder.reserve(child_capacity);
                    this->w_builder.shrink_to_fit();
                }
                if(this->maximal_repeat_check_vec.capacity() > node_capacity){
                    this->maximal_repeat_check_vec.reserve(node_capacity);
                }
            }
            void first_compute(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                w_builder.clear();

                this->_stnode_count += em.computeFirstLCPIntervalSet();
                em.move_st_internal_nodes(this->childVec, w_builder);
                /*
#if DEBUG                
                this->_RLBWTDS->checkLCPInterval(this->stnodeVec[0]);
#endif
*/

                this->maximal_repeat_check_vec.resize(1);
                this->maximal_repeat_check_vec[0] = true;
            }
            uint64_t get_child_rank(uint64_t i){
                return this->leftmost_child_bits_selecter(i + 1);
            }
            void computeNextLCPIntervalSet(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &inputSet, uint64_t start_index, uint64_t width, uint64_t rank, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {

                assert(inputSet.node_count() > 0);
                this->clear();



                uint64_t end_index = start_index + width - 1;
                assert(end_index < inputSet.node_count());

                uint64_t node_capacity = width * 2;
                uint64_t width_sum = inputSet.get_child_rank(end_index + 1) - inputSet.get_child_rank(start_index);
                uint64_t child_capacity = width_sum * 2;
                this->shrink(node_capacity, child_capacity);

                RINTERVAL intv;
                for (uint64_t i = start_index; i <= end_index; i++)
                {
                    uint64_t width = inputSet.get_width(i);
                    inputSet.get_stnode(i, intv);
                    this->_stnode_count += em.computeNextLCPIntervalSet(intv, inputSet.childVec, rank, width);
                    em.move_st_internal_nodes(this->childVec, w_builder);

                    rank += width;
                }
            }
            void build_bits()
            {
                w_builder.push_back(true);
                //sdsl::bit_vector v;
                leftmost_child_bits.resize(w_builder.size());
                for (uint64_t i = 0; i < w_builder.size(); i++)
                {
                    leftmost_child_bits[i] = w_builder[i];
                }
                //this->leftmost_child_bits.swap(v);

                sdsl::bit_vector::select_1_type b_sel(&leftmost_child_bits);
                leftmost_child_bits_selecter.set_vector(&leftmost_child_bits);
                leftmost_child_bits_selecter.swap(b_sel);
                this->builded = true;
            }

            void print()
            {
                std::cout << "[" << this->node_count() << ", " << this->children_count() << "]" << std::endl;
            }
            void print_info()
            {
                std::cout << "PRINT TREE" << std::endl;
                /*
                for (uint64_t i = 0; i < this->node_count(); i++)
                {
                    this->stnodeVec[i].print2(this->_RLBWTDS->_fposDS);
                }
                */
                std::cout << std::endl;
            }

            void add(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &item)
            {
                //assert(item.stnodeVec.size() == item._widthVec.size());
                //assert(item.stnodeVec.size() == item.maximal_repeat_check_vec.size());

                //assert(this->stnodeVec.size() == this->_widthVec.size());
                //assert(this->stnodeVec.size() == this->maximal_repeat_check_vec.size());
                for (uint64_t i = 0; i < item.maximal_repeat_check_vec.size(); i++)
                {
                    this->maximal_repeat_check_vec.push_back(item.maximal_repeat_check_vec[i]);
                }
                for (uint64_t i = 0; i < item.children_count(); i++)
                {
                    this->childVec.push_back(item.childVec[i]);
                }
                for (uint64_t i = 0; i < item.w_builder.size(); i++)
                {
                    this->w_builder.push_back(item.w_builder[i]);
                }
                this->_stnode_count += item._stnode_count;

                item.clear();
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool