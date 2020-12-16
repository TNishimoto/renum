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
        class SuccinctSortedSTChildren
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

            SuccinctSortedSTChildren()
            {
            }
            SuccinctSortedSTChildren(RLBWTDS *__RLBWTDS) : _RLBWTDS(__RLBWTDS)
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
            void print(){
                std::cout << "@SuccinctSortedSTChildren" << std::endl;
                for(uint64_t i=0;i<this->_stnode_count;i++){
                    uint64_t l = this->get_stnode_left_index(i);
                    uint64_t r = this->get_stnode_right_index(i);
                    uint64_t left = this->get_child_node_left(l);
                    uint64_t right = this->get_child_node_right(r);
                    std::cout << "[" << left << ", " << right << "], c = ";
                    for(uint64_t x = l;x<=r ;x++){
                        left = this->get_child_node_left(x);
                        right = this->get_child_node_right(x);
                        std::cout << "[" << left << ", " << right << "]";
                    }
                    if(i + 1 < this->_stnode_count){
                    std::cout << std::endl;
                    }

                }
                std::cout << "[END]" << std::endl;
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
            uint64_t get_stnode_start_rindex(uint64_t i)
            {
                assert(i < this->_stnode_count);

                uint64_t L = this->leftmost_child_bits_selecter(i + 1);
                uint64_t x = this->get_child_node_left(L);
                return this->_RLBWTDS->get_lindex_containing_the_position(x);
            }
            uint64_t get_stnode_end_rindex(uint64_t i)
            {
                assert(i < this->_stnode_count);

                uint64_t R = this->leftmost_child_bits_selecter(i + 2) - 1;
                uint64_t x = this->get_child_node_right(R);
                return this->_RLBWTDS->get_lindex_containing_the_position(x);
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

    } // namespace lcp_on_rlbwt
} // namespace stool