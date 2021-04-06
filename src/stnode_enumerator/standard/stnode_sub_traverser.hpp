#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <deque>
#include <vector>
#include <array>

#include <type_traits>
#include "../explicit_weiner_link_computer_on_rlbwt.hpp"
#include "../explicit_weiner_link_computer.hpp"

#include "../../../module/stool/src/io.h"
#include "../stnode_vector.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace stnode_on_rlbwt
    {

        template <typename INDEX_SIZE, typename RLBWTDS, typename CHAR = uint8_t>
        class STNodeSubTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;

            uint64_t _first_child_flag_vec_count = 0;
            uint64_t _maximal_repeat_check_vec_count = 0;
            std::vector<INDEX_SIZE> childs_vec;
            std::vector<CHAR> edge_char_vec;
            std::vector<bool> first_child_flag_vec;
            std::vector<bool> maximal_repeat_check_vec;
            bool store_edge_chars = false;
            //RLBWTDS *_RLBWTDS = nullptr;

        public:
            STNodeSubTraverser()
            {
            }
            STNodeSubTraverser(uint64_t size, bool _store_edge_chars)
            {
                this->childs_vec.resize(size * 2);
                this->first_child_flag_vec.resize(size * 2);
                this->maximal_repeat_check_vec.resize(size);
                this->store_edge_chars = _store_edge_chars;
                if (this->store_edge_chars)
                {
                    this->edge_char_vec.resize(size * 2);
                }
            }
            /*
            RLBWTDS *get_rlbwtds() const
            {
                return this->_RLBWTDS;
            }
            */
            uint64_t capacity()
            {
                return this->first_child_flag_vec.size();
            }

        private:
            uint64_t get_first_child_pointer() const
            {
                return 1;
            }

        public:
            inline uint64_t get_child_left_boundary(uint64_t child_end_pointer) const
            {
                assert(child_end_pointer > 0);
                assert(child_end_pointer < this->childvec_size());
                if (this->first_child_flag_vec[child_end_pointer - 1])
                {
                    return this->childs_vec[child_end_pointer - 1];
                }
                else
                {
                    return this->childs_vec[child_end_pointer - 1] + 1;
                }
            }
            inline uint64_t get_child_right_boundary(uint64_t child_end_pointer) const
            {
                assert(child_end_pointer < this->childvec_size());
                return this->childs_vec[child_end_pointer];
            }
            bool has_edge_characters() const
            {
                return this->store_edge_chars;
            }
            CHAR get_edge_character(INDEX_SIZE child_index) const
            {
                return this->edge_char_vec[child_index];
            }
            uint64_t get_integer_array_size() const
            {
                return this->_first_child_flag_vec_count;
            }
            uint64_t increment(uint64_t L, uint64_t &left, uint64_t &right) const
            {
                assert(L > 0);
                assert(this->first_child_flag_vec[L - 1]);
                assert(!this->first_child_flag_vec[L]);

                uint64_t R = L + 1;
                while (R < this->childvec_size() && !this->first_child_flag_vec[R])
                {
                    R++;
                }
                left = this->get_child_left_boundary(L);
                right = this->get_child_right_boundary(R - 1);

                return R + 1;
            }

            void get_lcp_intervals(uint64_t lcp, std::vector<stool::LCPInterval<uint64_t>> &output) const
            {
                uint64_t size = this->node_count();
                uint64_t L = this->get_first_child_pointer();
                uint64_t left;
                uint64_t right;
                for (uint64_t i = 0; i < size; i++)
                {
                    L = this->increment(L, left, right);
                    output.push_back(stool::LCPInterval<uint64_t>(left, right, lcp));
                }
            }
            void load(std::ifstream &file)
            {
                (void)file;
                assert(false);
                throw -1;
                /*
                file.read((char *)&this->_stnode_count, sizeof(uint64_t));

                //stool::IO::load(file, this->childs_vec, false);
                stool::IO::load_deque(file, this->childs_vec);

                stool::IO::load_bits(file, this->first_child_flag_vec);
                stool::IO::load_bits(file, this->maximal_repeat_check_vec);
                */
            }
            void write(std::ofstream &out)
            {
                (void)out;
                assert(false);
                throw -1;
                /*
                uint64_t _node_count = this->_stnode_count;
                out.write(reinterpret_cast<const char *>(&_node_count), sizeof(uint64_t));
                //stool::IO::write(out, this->childs_vec, false);
                stool::IO::write_deque(out, this->childs_vec);

                stool::IO::write_bits(out, this->first_child_flag_vec);
                stool::IO::write_bits(out, this->maximal_repeat_check_vec);
                */
            }

            uint64_t childvec_size() const
            {
                return this->_first_child_flag_vec_count;
            }

            uint64_t children_count() const
            {
                return this->_first_child_flag_vec_count - this->_maximal_repeat_check_vec_count;
            }
            uint64_t node_count() const
            {
                return this->_maximal_repeat_check_vec_count;
            }

            void clear()
            {
                this->_maximal_repeat_check_vec_count = 0;
                this->_first_child_flag_vec_count = 0;
            }
            void swap(STNodeSubTraverser<INDEX_SIZE, RLBWTDS> &item)
            {
                this->childs_vec.swap(item.childs_vec);
                this->first_child_flag_vec.swap(item.first_child_flag_vec);
                this->maximal_repeat_check_vec.swap(item.maximal_repeat_check_vec);
                this->edge_char_vec.swap(item.edge_char_vec);

                std::swap(this->_maximal_repeat_check_vec_count, item._maximal_repeat_check_vec_count);
                std::swap(this->_first_child_flag_vec_count, item._first_child_flag_vec_count);
                std::swap(this->store_edge_chars, item.store_edge_chars);

                //std::swap(this->_RLBWTDS, item._RLBWTDS);
            }

            uint64_t read_node(uint64_t L, std::pair<INDEX_SIZE, INDEX_SIZE> &output_node, std::vector<std::pair<INDEX_SIZE, INDEX_SIZE>> &output_children, std::vector<CHAR> &output_edge_chars)
            {

                uint64_t _left = 0, _right = 0;
                uint64_t nextL = this->increment(L, _left, _right);
                output_node.first = _left;
                output_node.second = _right;
                uint64_t _count = nextL - L - 1;

                uint64_t pointer = L;
                for (uint64_t i = 0; i < _count; i++)
                {
                    assert(pointer < this->childvec_size());

                    uint64_t left = this->get_child_left_boundary(pointer);
                    uint64_t right = this->get_child_right_boundary(pointer);
                    if (this->store_edge_chars)
                    {
                        CHAR c = this->get_edge_character(pointer);
                        output_edge_chars.push_back(c);
                    }
                    assert(left <= right);
                    output_children.push_back(std::pair<INDEX_SIZE, INDEX_SIZE>(left, right));
                    pointer++;
                }
                return nextL;
            }

            void computeNextSTNodes(ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS> &em, STNodeVector<INDEX_SIZE> &tmp)
            {
                tmp.clear();
                //RINTERVAL intv;
                //RINTERVAL child;

                uint64_t size = this->node_count();
                uint64_t L = this->get_first_child_pointer();

                std::pair<INDEX_SIZE, INDEX_SIZE> output_node;
                std::vector<std::pair<INDEX_SIZE, INDEX_SIZE>> output_children;
                //std::vector<uint8_t> output_chars;
                std::vector<CHAR> output_edge_chars;

                for (uint64_t i = 0; i < size; i++)
                {
                    output_children.clear();
                    //output_chars.clear();
                    output_edge_chars.clear();
                    assert(this->first_child_flag_vec[L - 1]);
                    L = this->read_node(L, output_node, output_children, output_edge_chars);
                    em.executeWeinerLinkSearch(output_node, output_children, this->store_edge_chars ? &output_edge_chars : nullptr, tmp);

                    /*
                    if (this->store_edge_chars)
                    {
                    }
                    else
                    {
                        em.executeWeinerLinkSearch(output_node, output_children, nullptr, output_chars);
                    }
                    */
                    //WeinerLinkCommonFunctions::output(this->store_edge_chars, tmp);
                }
            }

            std::pair<uint64_t, uint64_t> countNextLCPIntervalSet(ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS> &em)
            {
                STNodeVector<INDEX_SIZE> tmp;
                this->computeNextSTNodes(em, tmp);

                uint64_t next_st_count = tmp.maximal_repeat_check_vec.size();
                uint64_t next_children_count = tmp.first_child_flag_vec.size() - tmp.maximal_repeat_check_vec.size();

                return std::pair<uint64_t, uint64_t>(next_st_count, next_children_count);
            }

            void print(uint64_t spaceSize = 0)
            {
                std::string space = std::string(spaceSize, ' ');

                std::cout << space << "(STANDARD_SUB)[STNODE_COUNT, CHILDREN_COUNT] = [" << this->node_count() << ", " << this->children_count() << "]" << std::endl;
                STNodeVector<INDEX_SIZE> item;
                this->to_stnode_vector(item);

                stool::Printer::print("child_vec", item.childs_vec);
                stool::Printer::print_bits("first_child_flag_vec", item.first_child_flag_vec);
                stool::Printer::print_bits("maximal_repeat_check_vec", item.maximal_repeat_check_vec);
            }
            void print_info()
            {
                uint64_t L = this->get_first_child_pointer();
                uint64_t left = 0, right = 0;
                for (uint64_t i = 0; i < this->node_count(); i++)
                {
                    L = increment(L, left, right);
                    std::cout << "[" << left << ", " << right << "]";
                }

                std::cout << std::endl;
            }
            void to_stnode_vector(STNodeVector<INDEX_SIZE> &item)
            {
                for (uint64_t i = 0; i < this->childvec_size(); i++)
                {
                    item.childs_vec.push_back(this->childs_vec[i]);
                    item.first_child_flag_vec.push_back(this->first_child_flag_vec[i]);
                }
                if (this->store_edge_chars)
                {
                    for (uint64_t i = 0; i < this->edge_char_vec.size(); i++)
                    {
                        item.edge_char_vec.push_back(this->edge_char_vec[i]);
                    }
                }
                for (uint64_t i = 0; i < this->node_count(); i++)
                {
                    item.maximal_repeat_check_vec.push_back(this->maximal_repeat_check_vec[i]);
                }
            }
            uint64_t get_using_memory()
            {
                uint64_t x1 = this->childs_vec.size() * sizeof(INDEX_SIZE);
                uint64_t x2 = (this->first_child_flag_vec.size() * 1) / 8;
                //uint64_t x3 = (this->leftmost_child_bits.size() * 1) / 8;
                uint64_t x4 = (maximal_repeat_check_vec.size() * 1) / 8;
                return x1 + x2 + x4;
            }

            void import(STNodeVector<INDEX_SIZE> &item)
            {
                this->clear();
                for (uint64_t i = 0; i < item.childs_vec.size(); i++)
                {
                    this->childs_vec[i] = item.childs_vec[i];
                }
                for (uint64_t i = 0; i < item.first_child_flag_vec.size(); i++)
                {
                    this->first_child_flag_vec[i] = item.first_child_flag_vec[i];
                }
                if (this->store_edge_chars)
                {
                    for (uint64_t i = 0; i < item.edge_char_vec.size(); i++)
                    {
                        this->edge_char_vec[i] = item.edge_char_vec[i];
                    }
                }
                for (uint64_t i = 0; i < item.maximal_repeat_check_vec.size(); i++)
                {
                    this->maximal_repeat_check_vec[i] = item.maximal_repeat_check_vec[i];
                }

                this->_maximal_repeat_check_vec_count = item.maximal_repeat_check_vec.size();
                this->_first_child_flag_vec_count = item.first_child_flag_vec.size();
            }
            void move_push(STNodeVector<INDEX_SIZE> &item)
            {
                //assert(this->store_edge_chars == item.store_edge_chars);
                auto pair = item.compute_import_positions(this->capacity() - this->childvec_size());
                assert(pair.second > 0);

                uint64_t k = pair.first + pair.second;

                uint64_t p = item.first_child_flag_vec.size() - k;
                for (uint64_t i = 0; i < k; i++)
                {
                    this->childs_vec[_first_child_flag_vec_count] = item.childs_vec[p + i];
                    this->first_child_flag_vec[_first_child_flag_vec_count] = item.first_child_flag_vec[p + i];
                    if (this->store_edge_chars)
                    {
                        this->edge_char_vec[_first_child_flag_vec_count] = item.edge_char_vec[p + i];
                    }
                    this->_first_child_flag_vec_count++;
                }

                for (uint64_t i = 0; i < k; i++)
                {
                    item.childs_vec.pop_back();
                    item.first_child_flag_vec.pop_back();
                    if (this->store_edge_chars)
                    {
                        item.edge_char_vec.pop_back();
                    }
                }
                uint64_t q = item.maximal_repeat_check_vec.size() - pair.first;
                for (uint64_t i = 0; i < pair.first; i++)
                {
                    this->maximal_repeat_check_vec[_maximal_repeat_check_vec_count] = item.maximal_repeat_check_vec[q + i];
                    this->_maximal_repeat_check_vec_count++;
                }
                for (uint64_t i = 0; i < pair.first; i++)
                {
                    item.maximal_repeat_check_vec.pop_back();
                }
            }

            bool check_maximal_repeat(INDEX_SIZE node_index) const
            {
                return this->maximal_repeat_check_vec[node_index];
            }
        };

    } // namespace stnode_on_rlbwt
} // namespace stool