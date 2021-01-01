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
#include "../weiner_link_emulator.hpp"
#include "../../../module/stool/src/io.h"
#include "../stnode_vector.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {

        template <typename INDEX_SIZE, typename RLBWTDS>
        class STNodeSubTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;

            uint64_t _first_child_flag_vec_count = 0;
            uint64_t _maximal_repeat_check_vec_count = 0;
            std::vector<INDEX_SIZE> childs_vec;
            std::vector<bool> first_child_flag_vec;
            std::vector<bool> maximal_repeat_check_vec;
            RLBWTDS *_RLBWTDS = nullptr;

        public:
            STNodeSubTraverser()
            {
            }
            STNodeSubTraverser(uint64_t size, RLBWTDS *__RLBWTDS) : _RLBWTDS(__RLBWTDS)
            {
                this->childs_vec.resize(size * 2);
                this->first_child_flag_vec.resize(size);
                this->maximal_repeat_check_vec.resize(size);
            }
            RLBWTDS *get_rlbwtds() const
            {
                return this->_RLBWTDS;
            }
            uint64_t capacity()
            {
                return this->first_child_flag_vec.size();
            }
            void pop(std::deque<INDEX_SIZE> &item1, std::deque<bool> &item2, std::deque<bool> &item3)
            {
                assert(false);
                throw -1;
                /*
                uint64_t size = this->first_child_flag_vec.size();
                for (uint64_t i = 0; i < size; i++)
                {
                    item1.push_back(this->childs_vec[0]);
                    this->childs_vec.pop_front();

                    item1.push_back(this->childs_vec[0]);
                    this->childs_vec.pop_front();

                    item2.push_back(this->first_child_flag_vec[0]);
                    this->first_child_flag_vec.pop_front();
                }
                for (uint64_t i = 0; i < this->_stnode_count; i++)
                {
                    item3.push_back(this->maximal_repeat_check_vec[0]);
                    this->maximal_repeat_check_vec.pop_front();
                }
                this->_stnode_count = 0;
                */
            }

        private:
            inline uint64_t get_child_start_position(uint64_t child_end_pointer) const
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
            inline uint64_t get_child_end_position(uint64_t child_end_pointer) const
            {
                assert(child_end_pointer < this->childvec_size());
                return this->childs_vec[child_end_pointer];
            }
            uint64_t get_first_child_pointer() const
            {
                return 1;
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
                left = this->get_child_start_position(L);
                right = this->get_child_end_position(R - 1);

                return R + 1;
            }

        public:
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

                std::swap(this->_maximal_repeat_check_vec_count, item._maximal_repeat_check_vec_count);
                std::swap(this->_first_child_flag_vec_count, item._first_child_flag_vec_count);
                std::swap(this->_RLBWTDS, item._RLBWTDS);
            }
            uint64_t write_maximal_repeats(uint64_t lcp, std::ofstream &out) const
            {
                std::vector<stool::LCPInterval<INDEX_SIZE>> buffer;
                uint64_t size = this->node_count();
                uint64_t L = this->get_first_child_pointer();
                uint64_t left;
                uint64_t right;
                uint64_t count = 0;
                for (uint64_t i = 0; i < size; i++)
                {
                    L = this->increment(L, left, right);
                    if (this->maximal_repeat_check_vec[i])
                    {
                        buffer.push_back(stool::LCPInterval<INDEX_SIZE>(left, right, lcp));
                        count++;
                        if (buffer.size() >= 8000)
                        {
                            out.write(reinterpret_cast<const char *>(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                            buffer.clear();
                        }
                    }
                }

                if (buffer.size() >= 1)
                {
                    out.write(reinterpret_cast<const char *>(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                    buffer.clear();
                }
                return count;
            }

            void computeNextSTNodes(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em, STNodeVector<INDEX_SIZE> &tmp)
            {
                tmp.clear();
                RINTERVAL intv;
                RINTERVAL child;

                uint64_t size = this->node_count();
                uint64_t L = this->get_first_child_pointer();
                uint64_t _left = 0, _right = 0;

                for (uint64_t i = 0; i < size; i++)
                {
                    assert(this->first_child_flag_vec[L - 1]);

                    uint64_t nextL = this->increment(L, _left, _right);
                    this->_RLBWTDS->to_rinterval(_left, _right, intv);
                    uint64_t _count = nextL - L - 1;

                    em.clear();
#if DEBUG
                    if (this->_RLBWTDS->str_size() < 100)
                    {
                        std::cout << "Search input = ";
                        //intv.print2(this->_RLBWTDS->_fposDS);
                        std::cout << std::endl;
                    }
#endif
                    em.computeSTNodeCandidates(intv);

                    uint64_t pointer = L;

                    for (uint64_t i = 0; i < _count; i++)
                    {
                        assert(pointer < this->childvec_size());

                        uint64_t left = this->get_child_start_position(pointer);
                        uint64_t right = this->get_child_end_position(pointer);
                        assert(left <= right);

                        this->_RLBWTDS->to_rinterval(left, right, child);
                        em.computeSTChildren(child);
                        pointer++;
                    }

                    em.fit(false);

#if DEBUG
                    if (this->_RLBWTDS->stnc != nullptr)
                    {
                        em.verify_next_lcp_interval(_left, _right);
                    }
#endif
                    tmp.import(em, this->_RLBWTDS);
                    L = nextL;
                }
            }

            std::pair<uint64_t, uint64_t> countNextLCPIntervalSet(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                STNodeVector<INDEX_SIZE> tmp;
                this->computeNextSTNodes(em, tmp);

                uint64_t next_st_count = tmp.maximal_repeat_check_vec.size();
                uint64_t next_children_count = tmp.first_child_flag_vec.size() - tmp.maximal_repeat_check_vec.size();

                return std::pair<uint64_t, uint64_t>(next_st_count, next_children_count);
            }

            void print()
            {
                std::cout << "[STNODE_COUNT, CHILDREN_COUNT] = [" << this->node_count() << ", " << this->children_count() << "]" << std::endl;
                STNodeVector<INDEX_SIZE> item;
                this->convert(item);

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
            void convert(STNodeVector<INDEX_SIZE> &item)
            {
                for (uint64_t i = 0; i < this->childs_vec.size(); i++)
                {
                    item.childs_vec.push_back(this->childs_vec[i]);
                }
                for (uint64_t i = 0; i < this->first_child_flag_vec.size(); i++)
                {
                    item.first_child_flag_vec.push_back(this->first_child_flag_vec[i]);
                }
                for (uint64_t i = 0; i < this->maximal_repeat_check_vec.size(); i++)
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
                for (uint64_t i = 0; i < item.maximal_repeat_check_vec.size(); i++)
                {
                    this->maximal_repeat_check_vec[i] = item.maximal_repeat_check_vec[i];
                }

                this->_maximal_repeat_check_vec_count = item.maximal_repeat_check_vec.size();
                this->_first_child_flag_vec_count = item.first_child_flag_vec.size();
            }
            void move_push(STNodeVector<INDEX_SIZE> &item)
            {
                auto pair = item.compute_import_positions(this->capacity() - this->childvec_size());
                assert(pair.second > 0);

                uint64_t k = pair.first + pair.second;

                uint64_t p = item.first_child_flag_vec.size() - k;
                for (uint64_t i = 0; i < k; i++)
                {
                    this->childs_vec[_first_child_flag_vec_count] = item.childs_vec[p + i];
                    this->first_child_flag_vec[_first_child_flag_vec_count] = item.first_child_flag_vec[p + i];
                    this->_first_child_flag_vec_count++;
                }

                for (uint64_t i = 0; i < k; i++)
                {
                    item.childs_vec.pop_back();
                    item.first_child_flag_vec.pop_back();
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
        };

    } // namespace lcp_on_rlbwt
} // namespace stool