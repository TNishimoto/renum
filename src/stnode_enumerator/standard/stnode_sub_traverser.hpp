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

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        /*
        int64_t debug_sum_counter = 0;
        int64_t debug_peak_counter = 0;
        */
        //std::mutex mtx2;

        template <typename INDEX_SIZE, typename RLBWTDS>
        class STNodeSubTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;

            uint64_t _stnode_count = 0;
            std::deque<INDEX_SIZE> childs_vec;
            std::deque<bool> first_child_flag_vec;
            std::deque<bool> maximal_repeat_check_vec;
            RLBWTDS *_RLBWTDS = nullptr;

        public:
            STNodeSubTraverser()
            {
                maximal_repeat_check_vec.resize(0);
            }
            STNodeSubTraverser(RLBWTDS *__RLBWTDS) : _RLBWTDS(__RLBWTDS)
            {
                maximal_repeat_check_vec.resize(0);
            }
            /*
            ~STNodeSubTraverser()
            {
            }
            */
            void pop(std::deque<INDEX_SIZE> &item1, std::deque<bool> &item2, std::deque<bool> &item3)
            {
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

            void add(uint8_t c, uint64_t count, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em,
                     std::vector<INDEX_SIZE> &_child_vec, std::vector<bool> &_first_child_flag_vec, std::vector<bool> &_maximal_repeat_check_vec)
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

                    _child_vec.push_back(left);
                    _child_vec.push_back(right);
                    _first_child_flag_vec.push_back(j == 0);
                }
                uint64_t x = this->_RLBWTDS->get_lindex_containing_the_position(st_left);
                uint64_t d = this->_RLBWTDS->get_run(x);
                bool isMaximalRepeat = (this->_RLBWTDS->get_lpos(x) + d) <= st_right;
                _maximal_repeat_check_vec.push_back(isMaximalRepeat);
                this->_stnode_count++;
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
            void get_lcp_intervals(uint64_t lcp, std::vector<stool::LCPInterval<uint64_t>> &output)
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

            bool check_maximal_repeat(uint64_t st_index) const
            {
                return this->maximal_repeat_check_vec[st_index];
            }
            /*
            uint64_t get_stnode(uint64_t L, stool::LCPInterval<uint64_t> &output, uint64_t lcp)
            {
                if (!this->first_child_flag_vec[L])
                {
                    std::cout << L << std::endl;
                    this->print();
                }
                assert(this->first_child_flag_vec[L]);
                assert(!this->first_child_flag_vec[L + 1]);

                RINTERVAL tmp;
                uint64_t left = 0, right = 0;
                uint64_t newL = increment(L, left, right);

                output.i = left;
                output.j = right;
                output.lcp = lcp;

                return newL;
            }
            */

            uint64_t children_count() const
            {
                return this->childs_vec.size() / 2;
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
                //this->first_child_flag_vec.clear();

                //this->leftmost_child_bits.clear();
            }
            void swap(STNodeSubTraverser<INDEX_SIZE, RLBWTDS> &item)
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

                //em.move_st_internal_nodes(this->childs_vec, this->first_child_flag_vec);
                /*
                for (uint64_t i = 0; i < em.indexCount; i++)
                {
                    auto c = em.indexVec[i];
                    em.move_st_internal_nodes(this->childs_vec, this->first_child_flag_vec, this->maximal_repeat_check_vec, c);
                    this->_stnode_count++;
                }
                */
            }
            bool computeNextLCPIntervalSet(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {

                bool isSplit = false;
                RINTERVAL intv;

                assert(this->children_count() == this->first_child_flag_vec.size());
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

                assert(this->children_count() == this->first_child_flag_vec.size());

                return isSplit;
            }
            uint64_t write_maximal_repeats(uint64_t lcp, std::ofstream &out)
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
            bool computeNextLCPIntervalSetForParallelProcessing(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                std::vector<INDEX_SIZE> _child_vec;
                std::vector<bool> _first_child_flag_vec;
                std::vector<bool> _maximal_repeat_check_vec;
                RINTERVAL intv;

                assert(this->children_count() == this->first_child_flag_vec.size());
                uint64_t size = this->_stnode_count;
                uint64_t L = 0;
                uint64_t _left = 0, _right = 0;

                for (uint64_t i = 0; i < size; i++)
                {
                    assert(this->first_child_flag_vec[L]);

                    uint64_t nextL = this->increment(L, _left, _right);
                    this->_RLBWTDS->to_rinterval(_left, _right, intv);

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

                    for (uint64_t i = L; i < nextL; i++)
                    {
                        RINTERVAL child;
                        uint64_t left = this->get_child_start_position(i);
                        uint64_t right = this->get_child_end_position(i);
                        assert(left <= right);

                        this->_RLBWTDS->to_rinterval(left, right, child);
                        em.computeSTChildren(child);
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
                        this->add(c, count, em, _child_vec, _first_child_flag_vec, _maximal_repeat_check_vec);
                    }
                    L = nextL;
                }

                uint64_t minSize = std::min(this->first_child_flag_vec.size(), _first_child_flag_vec.size());
                for (uint64_t i = 0; i < minSize; i++)
                {
                    this->childs_vec[(i * 2)] = _child_vec[(i * 2)];
                    this->childs_vec[(i * 2) + 1] = _child_vec[(i * 2) + 1];
                    this->first_child_flag_vec[i] = _first_child_flag_vec[i];
                }
                for (uint64_t i = minSize; i < _first_child_flag_vec.size(); i++)
                {
                    this->childs_vec.push_back(_child_vec[(i * 2)]);
                    this->childs_vec.push_back(_child_vec[(i * 2) + 1]);
                    this->first_child_flag_vec.push_back(_first_child_flag_vec[i]);
                }
                while (this->first_child_flag_vec.size() > _first_child_flag_vec.size())
                {
                    this->childs_vec.pop_back();
                    this->childs_vec.pop_back();
                    this->first_child_flag_vec.pop_back();
                }

                uint64_t minSTSize = std::min(this->maximal_repeat_check_vec.size(), _maximal_repeat_check_vec.size());
                for (uint64_t i = 0; i < minSTSize; i++)
                {
                    this->maximal_repeat_check_vec[i] = _maximal_repeat_check_vec[i];
                }
                for (uint64_t i = minSTSize; i < _maximal_repeat_check_vec.size(); i++)
                {
                    this->maximal_repeat_check_vec.push_back(_maximal_repeat_check_vec[i]);
                }
                while (this->maximal_repeat_check_vec.size() > _maximal_repeat_check_vec.size())
                {
                    this->maximal_repeat_check_vec.pop_back();
                }
                this->_stnode_count = this->maximal_repeat_check_vec.size();

                assert(this->children_count() == this->first_child_flag_vec.size());

                return true;
            }
            /*
            std::pair<uint64_t, uint64_t> test(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em){
                auto st = this->countNextLCPIntervalSet(em);
                //this->clear();
                std::deque<INDEX_SIZE> tmp1;
                std::deque<bool> tmp2;
                while(tmp2.size() < st.second){
                    tmp1.push_back(0);
                    tmp1.push_back(0);
                    tmp2.push_back(true);
                }
                uint64_t minSize = std::min(this->first_child_flag_vec.size(), tmp2.size());
                for(uint64_t i=0;i<minSize;i++){
                    this->childs_vec[(i*2)] = tmp1[(i*2)];
                    this->childs_vec[(i*2)+1] = tmp1[(i*2)+1];
                    this->first_child_flag_vec[i] = tmp2[i];
                }
                for(uint64_t i=minSize;i<tmp2.size();i++){
                    this->childs_vec.push_back(tmp1[(i*2)]);
                    this->childs_vec.push_back(tmp1[(i*2)+1]);
                    this->first_child_flag_vec.push_back(tmp2[i]);
                }
                while(this->first_child_flag_vec.size() > st.second){
                    this->childs_vec.pop_back();
                    this->childs_vec.pop_back();
                    this->first_child_flag_vec.pop_back();
                }
                return st;

            }
            */

            std::pair<uint64_t, uint64_t> countNextLCPIntervalSet(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {

                RINTERVAL intv;

                uint64_t next_st_count = 0;
                uint64_t next_children_count = 0;

                assert(this->children_count() == this->first_child_flag_vec.size());
                uint64_t size = this->_stnode_count;
                uint64_t L = 0;
                uint64_t _left = 0, _right = 0;

                for (uint64_t i = 0; i < size; i++)
                {
                    assert(this->first_child_flag_vec[L]);

                    uint64_t nextL = this->increment(L, _left, _right);
                    this->_RLBWTDS->to_rinterval(_left, _right, intv);

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

                    for (uint64_t i = L; i < nextL; i++)
                    {
                        RINTERVAL child;
                        uint64_t left = this->get_child_start_position(i);
                        uint64_t right = this->get_child_end_position(i);
                        assert(left <= right);

                        this->_RLBWTDS->to_rinterval(left, right, child);
                        em.computeSTChildren(child);
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
                        next_st_count++;
                        next_children_count += count;
                    }
                    L = nextL;
                }

                return std::pair<uint64_t, uint64_t>(next_st_count, next_children_count);
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
            void merge(STNodeSubTraverser<INDEX_SIZE, RLBWTDS> &item)
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
            void split(STNodeSubTraverser<INDEX_SIZE, RLBWTDS> &item)
            {
                uint64_t k = this->_stnode_count / 2;
                uint64_t x = 0;
                while (this->first_child_flag_vec.size() > 0)
                {
                    if (this->first_child_flag_vec[0])
                    {
                        x++;
                        if (x > k)
                        {
                            break;
                        }
                    }
                    item.childs_vec.push_back(this->childs_vec[0]);
                    this->childs_vec.pop_front();

                    item.childs_vec.push_back(this->childs_vec[0]);
                    this->childs_vec.pop_front();

                    item.first_child_flag_vec.push_back(this->first_child_flag_vec[0]);
                    this->first_child_flag_vec.pop_front();
                }

                assert(x - 1 == k);
                for (uint64_t i = 0; i < k; i++)
                {
                    item.maximal_repeat_check_vec.push_back(this->maximal_repeat_check_vec[0]);
                    this->maximal_repeat_check_vec.pop_front();
                }
                item._stnode_count += k;
                this->_stnode_count -= k;
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool