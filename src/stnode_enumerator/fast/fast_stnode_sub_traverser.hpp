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
#include "../explicit_weiner_link_computer_on_rlbwt.hpp"

namespace stool
{
    namespace stnode_on_rlbwt
    {

        template <typename INDEX_SIZE, typename RLBWTDS, typename CHAR = uint8_t>
        class FastSTNodeSubTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;

            //uint64_t _stnode_count = 0;
            std::deque<INDEX_SIZE> childs_vec;
            std::deque<CHAR> edge_char_vec;

            std::deque<bool> first_child_flag_vec;
            std::deque<bool> maximal_repeat_check_vec;
            std::deque<uint64_t> child_width_vec;
            std::deque<uint64_t> stnode_count_vec;
            std::deque<bool> processed_flag_vec;
            uint64_t current_lcp = 0;
            uint64_t parent_current_lcp = 0;

            RLBWTDS *_RLBWTDS = nullptr;
            stool::stnode_on_rlbwt::RLE<CHAR> *rlbwt;

            bool store_edge_chars = false;

        public:
            FastSTNodeSubTraverser()
            {
            }
            FastSTNodeSubTraverser(RLBWTDS *__RLBWTDS, bool _store_edge_chars) : _RLBWTDS(__RLBWTDS), store_edge_chars(_store_edge_chars)
            {
                this->rlbwt = this->_RLBWTDS->get_rlbwt();
            }
            bool has_edge_characters() const {
                return this->store_edge_chars;
            }
            void set_parent_current_lcp(uint64_t lcp)
            {
                this->parent_current_lcp = lcp;
            }
            void import(STNodeVector<INDEX_SIZE> &item, uint64_t lcp, uint64_t num)
            {
                std::vector<stool::CharInterval<INDEX_SIZE, uint8_t>> tmp;
                for (uint64_t i = 0; i < num; i++)
                {
                    tmp.clear();
                    item.get_last(tmp);
                    for (uint64_t j = 0; j < tmp.size(); j++)
                    {
                        this->first_child_flag_vec.push_back(j == 0);
                        uint64_t lb = this->rlbwt->get_lindex_containing_the_position(tmp[j].i);
                        uint64_t ld = tmp[j].i - this->rlbwt->get_lpos(lb);
                        uint64_t rb = this->rlbwt->get_lindex_containing_the_position(tmp[j].j);
                        uint64_t rd = tmp[j].j - this->rlbwt->get_lpos(rb);
                        this->childs_vec.push_back(lb);
                        this->childs_vec.push_back(ld);
                        this->childs_vec.push_back(rb);
                        this->childs_vec.push_back(rd);
                        if(this->store_edge_chars){
                            this->edge_char_vec.push_back(tmp[j].c);
                        }
                    }
                    this->maximal_repeat_check_vec.push_back(item.maximal_repeat_check_vec[item.maximal_repeat_check_vec.size() - 1]);
                    item.pop();
                }

                assert(this->maximal_repeat_check_vec.size() == num);
                this->current_lcp = lcp;
                this->child_width_vec.push_back(this->first_child_flag_vec.size());
                this->stnode_count_vec.push_back(num);
                this->processed_flag_vec.push_back(false);
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
            inline CHAR read_child(uint64_t i, RINTERVAL &output) const
            {
                output.beginIndex = this->childs_vec[(i * 4)];
                output.beginDiff = this->childs_vec[(i * 4) + 1];
                output.endIndex = this->childs_vec[(i * 4) + 2];
                output.endDiff = this->childs_vec[(i * 4) + 3];
                return this->store_edge_chars ? this->edge_char_vec[i] : 0;
            }

            inline void add_child(RINTERVAL &output, CHAR c, bool isFirst)
            {
                this->childs_vec.push_back(output.beginIndex);
                this->childs_vec.push_back(output.beginDiff);
                this->childs_vec.push_back(output.endIndex);
                this->childs_vec.push_back(output.endDiff);
                if(this->store_edge_chars){
                    this->edge_char_vec.push_back(c);
                }
                this->first_child_flag_vec.push_back(isFirst);
            }
            inline void call_st_node(uint64_t count)
            {
                uint64_t endIndex = this->children_count() - 1;
                uint64_t startIndex = this->children_count() - 1 - (count - 1);

                uint64_t fst = this->childs_vec[(startIndex * 4)];
                uint64_t end = this->childs_vec[(endIndex * 4) + 2];
                bool isMaximalRepeat = fst != end;

                this->maximal_repeat_check_vec.push_back(isMaximalRepeat);
            }
            inline void call_level_nodes(uint64_t stnode_count, uint64_t childnode_count)
            {
                this->processed_flag_vec.push_back(stnode_count == 0);
                this->stnode_count_vec.push_back(stnode_count);
                this->child_width_vec.push_back(childnode_count);
            }

            void add(uint8_t c, uint64_t count, ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS> &em)
            {
                RINTERVAL copy;

                for (uint64_t j = 0; j < count; j++)
                {
                    CHAR edge_char = em.get_child(c, j, copy);

                    this->add_child(copy, edge_char, j == 0);
                }
                this->call_st_node(count);
            }

            uint64_t children_count() const
            {
                return this->first_child_flag_vec.size();
            }
            uint64_t node_count() const
            {
                return this->maximal_repeat_check_vec.size();
            }

        public:
            uint64_t get_current_lcp() const
            {
                return this->current_lcp;
            }
            uint64_t get_last_lcp() const
            {
                return current_lcp + this->child_width_vec.size() - 1;
            }

            void pop_level(std::deque<INDEX_SIZE> &item1, std::deque<bool> &item2, std::deque<bool> &item3, std::deque<CHAR> &output_edge_char_vec)
            {

                assert(this->child_width_vec.size() >= 2);
                for (uint64_t i = 0; i < this->child_width_vec[0]; i++)
                {
                    item1.push_back(this->childs_vec[0]);
                    this->childs_vec.pop_front();

                    item1.push_back(this->childs_vec[0]);
                    this->childs_vec.pop_front();

                    item1.push_back(this->childs_vec[0]);
                    this->childs_vec.pop_front();

                    item1.push_back(this->childs_vec[0]);
                    this->childs_vec.pop_front();

                    item2.push_back(this->first_child_flag_vec[0]);
                    this->first_child_flag_vec.pop_front();

                    if(this->store_edge_chars){                        
                        output_edge_char_vec.push_back(this->edge_char_vec[0]);
                        this->edge_char_vec.pop_front();
                    }
                }
                for (uint64_t i = 0; i < this->stnode_count_vec[0]; i++)
                {
                    item3.push_back(this->maximal_repeat_check_vec[0]);
                    this->maximal_repeat_check_vec.pop_front();
                }
                child_width_vec.pop_front();
                stnode_count_vec.pop_front();
                processed_flag_vec.pop_front();
                this->current_lcp++;
                assert(this->child_width_vec.size() >= 1);
            }

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
            uint64_t current_children_count() const
            {
                if (this->child_width_vec.size() == 0)
                {
                    return 0;
                }
                else
                {
                    return this->child_width_vec[0];
                }
            }
            uint64_t current_node_count() const
            {
                if (this->stnode_count_vec.size() == 0)
                {
                    return 0;
                }
                else
                {
                    return this->stnode_count_vec[0];
                }
            }
            uint64_t last_node_count() const
            {
                if (this->stnode_count_vec.size() == 0)
                {
                    return 0;
                }
                else
                {
                    return this->stnode_count_vec[this->stnode_count_vec.size() - 1];
                }
            }
            uint64_t finished_level_count()
            {
                if (this->child_width_vec.size() == 0)
                {
                    return 0;
                }
                else
                {
                    return this->child_width_vec.size() - 1;
                }
            }

            void clear()
            {
                //this->_stnode_count = 0;
                this->childs_vec.clear();
                this->maximal_repeat_check_vec.clear();
                this->first_child_flag_vec.clear();
            }
            void first_compute(ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS> &em)
            {
                this->first_child_flag_vec.clear();

                //em.computeFirstLCPIntervalSet();
                auto r = em.getFirstChildren();
                RINTERVAL intv;
                for (uint64_t i = 0; i < r.size(); i++)
                {
                    auto &it = r[i];

                    this->_RLBWTDS->to_rinterval(it.i, it.j, intv);
                    this->add_child(intv, it.c, i == 0);
                }
                this->call_st_node(r.size());
                this->call_level_nodes(1, r.size());
            }
            void computeNextLCPIntervalSet(ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS> &em, uint64_t limit)
            {
                if (this->child_width_vec.size() == 0)
                {
                    return;
                }
                uint64_t level = this->child_width_vec.size() - 1;
                uint64_t index = this->children_count() - this->child_width_vec[this->child_width_vec.size() - 1];

                uint64_t processed_level_counter = 0;

                while (level < this->child_width_vec.size())
                {

                    if (!this->processed_flag_vec[level])
                    {

                        this->computeNextLCPIntervalSet(em, index, level);
                        processed_level_counter++;
                    }
                    index += this->child_width_vec[level];

                    level++;
                    if (processed_level_counter >= 1 && this->children_count() > limit)
                    {
                        break;
                    }
                }
            }

            uint64_t read_st_node(uint64_t L, RINTERVAL &node, std::vector<RINTERVAL> &children, std::vector<CHAR> &output_chars)
            {
                uint64_t nextL = this->read_st_node(L, node);
                RINTERVAL child;
                while (L < nextL)
                {
                    CHAR c = this->read_child(L, child);
                    L++;
                    children.push_back(child);
                    if(this->store_edge_chars){
                    output_chars.push_back(c);
                    }
                }
                return L;
            }

            uint64_t computeNextLCPIntervalSet(ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS> &em, uint64_t index, uint64_t level_index)
            {
                assert(!this->processed_flag_vec[level_index]);
                RINTERVAL intv;
                std::vector<RINTERVAL> children;
                std::vector<uint8_t> output_chars;
                std::vector<CHAR> output_edge_chars;

                assert(this->children_count() == this->first_child_flag_vec.size());
                uint64_t c_count = this->child_width_vec[level_index];
                uint64_t L = index;
                assert(this->first_child_flag_vec[L]);

                uint64_t addedChildrenNodeCounter = 0;
                uint64_t addedSTNodeCounter = 0;

                uint64_t currentProcessCount = 0;

                while (currentProcessCount < c_count)
                {
                    children.clear();
                    output_chars.clear();
                    output_edge_chars.clear();
                    assert(this->first_child_flag_vec[L]);
                    L = this->read_st_node(L, intv, children, output_edge_chars);
                    currentProcessCount += children.size();
                    #if DEBUG
                        em.checker_on = false;
                    #endif

                    em.executeWeinerLinkSearch(intv, children, this->store_edge_chars ? &output_edge_chars : nullptr, output_chars);

                    for (auto &c : output_chars)
                    {
                        auto &currentVec = em.childrenVec[c];
                        uint64_t count = currentVec.size();
                        this->add(c, count, em);
                        addedChildrenNodeCounter += count;
                        addedSTNodeCounter++;
                    }
                }
                this->processed_flag_vec[level_index] = true;
                this->call_level_nodes(addedSTNodeCounter, addedChildrenNodeCounter);
                //assert(this->first_child_flag_vec[L]);

                return L;
            }

            void print()
            {
                std::cout << "[STNODE_COUNT, CHILDREN_COUNT] = [" << this->node_count() << ", " << this->children_count() << "]" << std::endl;
                stool::Printer::print_bits("processed_flag_vec", this->processed_flag_vec);
                stool::Printer::print_bits("first_child_flag_vec", this->first_child_flag_vec);
                stool::Printer::print_bits("maximal_repeat_check_vec", this->maximal_repeat_check_vec);
                stool::Printer::print("child_vec", this->childs_vec);

                stool::Printer::print("child_width_vec", this->child_width_vec);
                stool::Printer::print("stnode_count_vec", this->stnode_count_vec);

                RINTERVAL output;
                output.beginIndex = std::numeric_limits<INDEX_SIZE>::max();
                output.endIndex = std::numeric_limits<INDEX_SIZE>::max();
                output.beginDiff = std::numeric_limits<INDEX_SIZE>::max();
                output.endDiff = std::numeric_limits<INDEX_SIZE>::max();

                uint64_t L = 0;
                for (uint64_t i = 0; i < this->stnode_count_vec.size(); i++)
                {
                    for (uint64_t j = 0; j < this->stnode_count_vec[i]; j++)
                    {
                        L = this->read_st_node(L, output);

                        auto intv = output.get_lcp_interval(this->current_lcp + i, *this->_RLBWTDS->get_lpos_vec_pointer());
                        std::cout << intv.to_string();
                    }
                    std::cout << std::endl;
                }
            }
            void print_info()
            {
                uint64_t L = 0;
                RINTERVAL intv;
                for (uint64_t i = 0; i < this->current_node_count(); i++)
                {
                    L = read_st_node(L, intv);
                    intv.print2(*this->_RLBWTDS->get_lpos_vec_pointer());
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
                assert(this->node_count() == k);
            }
#endif
            void swap(FastSTNodeSubTraverser<INDEX_SIZE, RLBWTDS> &item)
            {
                this->childs_vec.swap(item.childs_vec);
                this->first_child_flag_vec.swap(item.first_child_flag_vec);

                this->maximal_repeat_check_vec.swap(item.maximal_repeat_check_vec);

                std::swap(this->_RLBWTDS, item._RLBWTDS);
                std::swap(this->rlbwt, item.rlbwt);

            }
        public:
            bool is_empty()
            {
                return this->first_child_flag_vec.size() == 0;
            }
            void split(FastSTNodeSubTraverser<INDEX_SIZE, RLBWTDS> &item)
            {
                //this->print();

                uint64_t size = this->stnode_count_vec[this->stnode_count_vec.size() - 1];
                uint64_t split_size = size / 2;
                uint64_t x = this->first_child_flag_vec.size();
                uint64_t p = 0;
                uint64_t c_count = 0;
                uint64_t st_count = this->node_count();
                while (p < split_size)
                {
                    x--;
                    c_count++;
                    if (this->first_child_flag_vec[x])
                    {
                        p++;
                    }
                }

                uint64_t start = st_count - p;
                for (uint64_t i = start; i < st_count; i++)
                {
                    item.maximal_repeat_check_vec.push_back(this->maximal_repeat_check_vec[i]);
                }

                for (uint64_t i = x; i < this->first_child_flag_vec.size(); i++)
                {
                    item.childs_vec.push_back(this->childs_vec[(i * 4)]);
                    item.childs_vec.push_back(this->childs_vec[(i * 4) + 1]);
                    item.childs_vec.push_back(this->childs_vec[(i * 4) + 2]);
                    item.childs_vec.push_back(this->childs_vec[(i * 4) + 3]);
                    item.first_child_flag_vec.push_back(this->first_child_flag_vec[i]);
                }
                item.child_width_vec.push_back(this->first_child_flag_vec.size() - x);
                item.stnode_count_vec.push_back(split_size);
                item.processed_flag_vec.push_back(false);
                item.current_lcp = this->current_lcp + this->stnode_count_vec.size() - 1;

                for (uint64_t i = start; i < st_count; i++)
                {
                    this->maximal_repeat_check_vec.pop_back();
                }
                uint64_t pop_count = this->first_child_flag_vec.size() - x;
                for (uint64_t i = 0; i < pop_count; i++)
                {
                    this->childs_vec.pop_back();
                    this->childs_vec.pop_back();
                    this->childs_vec.pop_back();
                    this->childs_vec.pop_back();
                    this->first_child_flag_vec.pop_back();
                }
                this->child_width_vec[this->child_width_vec.size() - 1] -= c_count;
                this->stnode_count_vec[this->stnode_count_vec.size() - 1] -= split_size;

                /*
                std::cout << "THSI" << std::endl;

                this->print();
                this->bit_check();
                std::cout << "ITEM" << std::endl;
                item.print();
                */

                assert(!this->first_child_flag_vec[this->first_child_flag_vec.size() - 1]);
                assert(this->processed_flag_vec.size() == this->child_width_vec.size());
            }
        };

    } // namespace stnode_on_rlbwt
} // namespace stool