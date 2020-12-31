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
#include "../../../module/stool/src/io.h"
#include "../stnode_vector.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        template <typename INDEX_SIZE, typename RLBWTDS>
        class SingleSTNodeTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;

            uint64_t _stnode_count = 0;
            std::deque<INDEX_SIZE> childs_vec;
            std::deque<bool> first_child_flag_vec;
            std::deque<bool> maximal_repeat_check_vec;
            int64_t lcp = -1;
            RLBWTDS *_RLBWTDS = nullptr;
            ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> em;

        public:
            SingleSTNodeTraverser()
            {
            }
            SingleSTNodeTraverser(RLBWTDS *__RLBWTDS) : _RLBWTDS(__RLBWTDS)
            {
            }
            void initialize(RLBWTDS *__RLBWTDS)
            {
                this->_RLBWTDS = __RLBWTDS;
                this->lcp = -1;
                this->em.initialize(__RLBWTDS);
            }

        private:
            inline uint64_t get_child_start_position(uint64_t i) const
            {
                return this->childs_vec[(i * 2)];
            }
            inline uint64_t get_child_end_position(uint64_t i) const
            {
                return this->childs_vec[(i * 2) + 1];
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

                    if (left < st_left)
                    {
                        st_left = left;
                    }
                    if (right > st_right)
                    {
                        st_right = right;
                    }

                    this->childs_vec.push_back(left);
                    this->childs_vec.push_back(right);
                    this->first_child_flag_vec.push_back(j == 0);
                }
                uint64_t x = this->_RLBWTDS->get_lindex_containing_the_position(st_left);
                uint64_t d = this->_RLBWTDS->get_run(x);
                bool isMaximalRepeat = (this->_RLBWTDS->get_lpos(x) + d) <= st_right;
                this->maximal_repeat_check_vec.push_back(isMaximalRepeat);
                this->_stnode_count++;
            }
            uint64_t increment(uint64_t L, uint64_t &left, uint64_t &right) const
            {
                assert(this->first_child_flag_vec[L]);
                assert(!this->first_child_flag_vec[L + 1]);

                uint64_t R = L + 1;
                while (R < this->first_child_flag_vec.size() && !this->first_child_flag_vec[R])
                {
                    R++;
                }

                R--;

                left = this->get_child_start_position(L);
                right = this->get_child_end_position(R);

                return R + 1;
            }

        public:
            void get_lcp_intervals(std::vector<stool::LCPInterval<uint64_t>> &output)
            {
                uint64_t size = this->maximal_repeat_check_vec.size();
                uint64_t L = 0;
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
                file.read((char *)&this->_stnode_count, sizeof(uint64_t));

                stool::IO::load_deque(file, this->childs_vec);
                stool::IO::load_bits(file, this->first_child_flag_vec);
                stool::IO::load_bits(file, this->maximal_repeat_check_vec);
            }
            void write(std::ofstream &out)
            {
                uint64_t _node_count = this->_stnode_count;
                out.write(reinterpret_cast<const char *>(&_node_count), sizeof(uint64_t));
                stool::IO::write_deque(out, this->childs_vec);
                stool::IO::write_bits(out, this->first_child_flag_vec);
                stool::IO::write_bits(out, this->maximal_repeat_check_vec);
            }

            int64_t get_current_lcp() const
            {
                return lcp;
            }
            uint64_t child_count() const
            {
                return this->first_child_flag_vec.size();
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
            void process()
            {
                if (this->lcp == -1)
                {
                    this->first_compute();
                }
                else
                {
                    this->computeNextLCPIntervalSet();
                }
            }
            void first_compute()
            {
                this->first_child_flag_vec.clear();

                //em.computeFirstLCPIntervalSet();
                auto r = em.getFirstChildren();
                for (uint64_t i = 0; i < r.size(); i++)
                {
                    auto &it = r[i];
                    this->childs_vec.push_back(it.first);
                    this->childs_vec.push_back(it.second);
                    this->first_child_flag_vec.push_back(i == 0);
                }
                this->maximal_repeat_check_vec.push_back(true);
                this->_stnode_count++;
                this->lcp = 0;
            }
            bool computeNextLCPIntervalSet()
            {

                bool isSplit = false;
                RINTERVAL intv;

                uint64_t size = this->_stnode_count;
                uint64_t L = 0;
                uint64_t _left = 0, _right = 0;

                for (uint64_t i = 0; i < size; i++)
                {
                    assert(this->first_child_flag_vec[L]);

                    uint64_t nextL = this->increment(L, _left, _right);
                    assert(_left <= _right);
                    this->_RLBWTDS->to_rinterval(_left, _right, intv);

                    //this->get_stnode(0, x - 1, intv);
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

                    for (uint64_t i = 0; i < nextL; i++)
                    {
                        RINTERVAL child;
                        uint64_t left = this->get_child_start_position(0);
                        uint64_t right = this->get_child_end_position(0);
                        assert(left <= right);

                        this->_RLBWTDS->to_rinterval(left, right, child);
                        em.computeSTChildren(child);
                        this->childs_vec.pop_front();
                        this->childs_vec.pop_front();

                        this->first_child_flag_vec.pop_front();
                    }

                    em.fit(false);
#if DEBUG
                    if (this->_RLBWTDS->stnc != nullptr)
                    {
                        em.verify_next_lcp_interval(_left, _right);
                    }
#endif
                    for (uint64_t i = 0; i < em.indexCount; i++)
                    {
                        auto c = em.indexVec[i];
                        auto &currentVec = em.childrenVec[c];
                        uint64_t count = currentVec.size();
                        this->add(c, count, em);
                    }
                    this->maximal_repeat_check_vec.pop_front();
                    this->_stnode_count--;
                }

                this->lcp++;

                return isSplit;
            }
            uint64_t write_maximal_repeats(std::ofstream &out)
            {
                std::vector<stool::LCPInterval<INDEX_SIZE>> buffer;
                uint64_t size = this->maximal_repeat_check_vec.size();
                uint64_t L = 0;
                uint64_t left;
                uint64_t right;
                uint64_t count = 0;
                for (uint64_t i = 0; i < size; i++)
                {
                    L = this->increment(L, left, right);
                    if (this->maximal_repeat_check_vec[i])
                    {
                        buffer.push_back(stool::LCPInterval<INDEX_SIZE>(left, right, this->lcp));
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

            void print()
            {
                std::cout << "[STNODE_COUNT, CHILDREN_COUNT] = [" << this->_stnode_count << ", " << this->first_child_flag_vec.size() << "]" << std::endl;
                stool::Printer::print("child_vec", this->childs_vec);
                stool::Printer::print_bits("first_child_flag_vec", this->first_child_flag_vec);
                stool::Printer::print_bits("maximal_repeat_check_vec", this->maximal_repeat_check_vec);
            }
            void print_info()
            {
                uint64_t L = 0;
                uint64_t left = 0, right = 0;
                for (uint64_t i = 0; i < this->node_count(); i++)
                {

                    //uint64_t left = this->get_stnode2_start_position(L);
                    //uint64_t right = this->get_stnode2_end_position(L);
                    L = increment(L, left, right);
                    std::cout << "[" << left << ", " << right << "]";
                }

                std::cout << std::endl;
                /*
                for(auto &it : this->childs_vec){
                    std::cout << it << ", " << std::flush;
                }
                */

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
            void convert_to_vector(STNodeVector<INDEX_SIZE> &output)
            {
                output.childs_vec.reserve(this->childs_vec.size());
                while (this->childs_vec.size() > 0)
                {
                    output.childs_vec.push_back(this->childs_vec[0]);
                    this->childs_vec.pop_front();
                }

                output.first_child_flag_vec.reserve(this->first_child_flag_vec.size());
                while (this->first_child_flag_vec.size() > 0)
                {
                    output.first_child_flag_vec.push_back(this->first_child_flag_vec[0]);
                    this->first_child_flag_vec.pop_front();
                }
                output.maximal_repeat_check_vec.reserve(this->maximal_repeat_check_vec.size());

                while (this->maximal_repeat_check_vec.size() > 0)
                {
                    output.maximal_repeat_check_vec.push_back(this->maximal_repeat_check_vec[0]);
                    this->maximal_repeat_check_vec.pop_front();
                }
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool