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
#include "stool/include/io.hpp"
#include "../stnode_vector.hpp"
#include "../node_iterator.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace renum
    {
        template <typename INDEX_SIZE, typename INTERVAL_SEARCH, typename CHAR = uint8_t>
        class DFSTraverser
        {
            /*
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using ITERATOR = STNodeIterator<SingleSTNodeTraverser>;
            using DEPTH_ITERATOR = STDepthIterator<SingleSTNodeTraverser>;
            */

            uint64_t current_stnode_id = 0;
            stool::renum::STNodeVector<INDEX_SIZE, CHAR> stack;
            bool store_edge_chars = false;
            INTERVAL_SEARCH *em;

        public:
            class iterator
            {
                uint64_t _current_stnode_id = 0;

                DFSTraverser *traverser;

            public:
                iterator() = default;

                iterator(DFSTraverser *_traverser, bool is_start) : traverser(_traverser)
                {

                    if (is_start)
                    {
                        this->_current_stnode_id = _traverser->get_current_stnode_id();
                    }
                    else
                    {
                        this->_current_stnode_id = UINT64_MAX;
                    }
                }
                iterator &operator++()
                {
                    traverser->compute();
                    _current_stnode_id = traverser->get_current_stnode_id();
                    return *this;
                }

                iterator &operator++(int)
                {
                    ++(*this);
                    return *this;
                }

                bool operator!=(const iterator &rhs) const
                {
                    return this->_current_stnode_id != rhs._current_stnode_id;
                }
                void print() const
                {
                    uint64_t left = this->get_left();
                    uint64_t right = this->get_right();
                    uint64_t lcp = this->get_lcp();

                    std::cout << "[" << left << ", " << right << ", " << lcp << "]" << std::endl;
                }

                INDEX_SIZE get_children_count() const
                {
                    return traverser->get_children_count(*this);
                }
                INDEX_SIZE get_left() const
                {
                    return traverser->stack.get_last_left_boundary();
                }
                INDEX_SIZE get_right() const
                {
                    return traverser->stack.get_last_right_boundary();
                }
                INDEX_SIZE get_lcp() const
                {
                    return traverser->stack.get_last_depth();
                }
                bool is_maximal_repeat()
                {
                    return traverser->stack.maximal_repeat_check_vec[traverser->stack.maximal_repeat_check_vec.size() - 1];
                }
                uint64_t child_count() const
                {
                    return traverser->stack.get_last_width() - 1;
                }

                /*
                INDEX_SIZE get_child_left_boundary(uint64_t ith_child) const
                {
                    return traverser->get_child_left_boundary(*this, ith_child);
                }
                INDEX_SIZE get_child_right_boundary(uint64_t ith_child) const
                {
                    return traverser->get_child_right_boundary(*this, ith_child);
                }

                INDEX_SIZE get_edge_character(uint64_t ith_child) const
                {
                    return traverser->get_edge_character(*this, ith_child);
                }
                bool is_maximal_repeat() const
                {
                    return traverser->check_maximal_repeat(*this);
                }

                std::pair<INDEX_SIZE, INDEX_SIZE> operator*() const
                {
                    uint64_t left = traverser->get_left(*this);
                    uint64_t right = traverser->get_right(*this);
                    return std::pair<INDEX_SIZE, INDEX_SIZE>(left, right);
                }

                
                bool has_edge_characters() const
                {
                    return this->traverser->has_edge_characters();
                }
                */
            };
            using index_type = INDEX_SIZE;

            INTERVAL_SEARCH *get_interval_search_deta_structure() const
            {
                return this->em;
            }

            bool has_edge_characters() const
            {
                return this->store_edge_chars;
            }

            DFSTraverser()
            {
            }
            void initialize(INTERVAL_SEARCH *_interval_search, bool _store_edge_chars)
            {
                this->stack.clear();
                this->em = _interval_search;
                this->store_edge_chars = _store_edge_chars;
                this->first_compute();
                /*
                this->em = new ExplicitWeinerLinkComputerOnRLBWT<RLBWTDS>();
                this->em->initialize(__RLBWTDS);
                */
            }

            void first_compute()
            {
                std::vector<CharInterval<INDEX_SIZE, uint8_t>> r = em->getFirstChildren();
                uint64_t k = 0;
                for (uint64_t i = 0; i < r.size(); i++)
                {
                    auto &it = r[i];
                    if (i == 0)
                    {
                        stack.childs_vec.push_back(it.i);
                        stack.first_child_flag_vec.push_back(true);
                        if (store_edge_chars)
                        {
                            stack.edge_char_vec.push_back(0);
                        }
                    }
                    stack.childs_vec.push_back(it.j);
                    stack.first_child_flag_vec.push_back(false);
                    k++;
                    if (store_edge_chars)
                    {
                        stack.edge_char_vec.push_back(it.c);
                    }
                }
                stack.maximal_repeat_check_vec.push_back(true);
                stack.depth_vec.push_back(0);
                stack.width_vec.push_back(k);
                assert(stack.first_child_flag_vec[0]);
                this->current_stnode_id = 0;
            }
            uint64_t get_current_stnode_id() const
            {
                return this->current_stnode_id;
            }
            bool compute()
            {
                em->executeWeinerLinkSearch(stack);
                bool b = stack.maximal_repeat_check_vec.size() > 0;
                if (b)
                {
                    this->current_stnode_id++;
                }
                else
                {
                    this->current_stnode_id = UINT64_MAX;
                }
                return b;
            }
            uint64_t get_stack_size() const
            {
                return this->stack.maximal_repeat_check_vec.size();
            }
            DFSTraverser::iterator begin()
            {
                this->initialize(this->em, this->store_edge_chars);
                return DFSTraverser::iterator(this, true);
            }
            DFSTraverser::iterator end()
            {
                return DFSTraverser::iterator(this, false);
            }

            uint64_t get_input_text_length()
            {
                return this->em->get_input_text_length();
            }
            /*
            void clear(){
                this->stack.clear();
                this->current_stnode_id = 0;
                this->first_compute();
            }
            */

            /*
            bool is_finished() const
            {
                return this->get_current_lcp() >= 0 && this->node_count() == 0;
            }

            int64_t get_current_lcp() const
            {
                return lcp;
            }
            uint64_t get_integer_array_size() const
            {
                return this->first_child_flag_vec.size();
            }
            uint64_t child_count() const
            {
                return this->first_child_flag_vec.size() - this->node_count();
            }
            uint64_t node_count() const
            {
                return this->_stnode_count;
            }
            DEPTH_ITERATOR begin()
            {
                this->clear();
                this->succ();
                return DEPTH_ITERATOR(this, true);
            }
            DEPTH_ITERATOR end()
            {
                return DEPTH_ITERATOR(this, false);
            }
            ITERATOR node_end_iterator()
            {
                return ITERATOR(this, false);
            }

            void clear()
            {
                this->_stnode_count = 0;
                this->childs_vec.clear();
                this->maximal_repeat_check_vec.clear();
                this->first_child_flag_vec.clear();
            }
            bool succ()
            {
                if (this->is_finished())
                {
                    return false;
                }
                else
                {

                    if (this->lcp == -1)
                    {
                        this->first_compute();
                    }
                    else
                    {
                        this->computeNextLCPIntervalSet();
                    }
                    return true;
                }
            }
            uint64_t get_first_child_pointer() const
            {
                return 1;
            }

            uint64_t get_input_text_length()
            {
                return this->em->get_input_text_length();
            }
            void print() const
            {
                std::cout << "[STNODE_COUNT, CHILDREN_COUNT] = [" << this->node_count() << ", " << this->child_count() << "]" << std::endl;
                stool::Printer::print("child_vec", this->childs_vec);
                stool::Printer::print_bits("first_child_flag_vec", this->first_child_flag_vec);
                stool::Printer::print_bits("maximal_repeat_check_vec", this->maximal_repeat_check_vec);
            }
            void print_info() const
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
            uint64_t get_using_memory() const
            {
                uint64_t x1 = this->childs_vec.size() * sizeof(INDEX_SIZE);
                uint64_t x2 = (this->first_child_flag_vec.size() * 1);
                //uint64_t x3 = (this->leftmost_child_bits.size() * 1) / 8;
                uint64_t x4 = (maximal_repeat_check_vec.size() * 1);
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

                if (this->store_edge_chars)
                {
                    //output.store_edge_chars = true;
                    output.edge_char_vec.reserve(this->edge_char_vec.size());
                    while (this->edge_char_vec.size() > 0)
                    {
                        output.edge_char_vec.push_back(this->edge_char_vec[0]);
                        this->edge_char_vec.pop_front();
                    }
                }
            }

        private:
            

            void pop_node(std::pair<INDEX_SIZE, INDEX_SIZE> &output_node, std::vector<std::pair<INDEX_SIZE, INDEX_SIZE>> &output_children, std::vector<CHAR> &output_edge_chars)
            {
                assert(this->first_child_flag_vec[0]);

                uint64_t L = 1;
                uint64_t _left = 0, _right = 0;
                uint64_t nextL = this->increment(L, _left, _right);
                output_node.first = _left;
                output_node.second = _right;
                uint64_t _count = nextL - L - 1;
                for (uint64_t i = 0; i < _count; i++)
                {
                    uint64_t left = this->get_child_left_boundary(L);
                    uint64_t right = this->get_child_right_boundary(L);
                    if (this->store_edge_chars)
                    {
                        CHAR c = this->get_edge_character(L);
                        output_edge_chars.push_back(c);
                        this->edge_char_vec.pop_front();
                    }
                    output_children.push_back(std::pair<INDEX_SIZE, INDEX_SIZE>(left, right));

                    assert(left <= right);
                    this->childs_vec.pop_front();
                    this->first_child_flag_vec.pop_front();
                }
                this->childs_vec.pop_front();
                this->first_child_flag_vec.pop_front();
                if (store_edge_chars)
                {
                    this->edge_char_vec.pop_front();
                }

                this->maximal_repeat_check_vec.pop_front();
                this->_stnode_count--;
            }

            void import(stool::renum::STNodeVector<INDEX_SIZE, CHAR> &item)
            {
                for (auto &it : item.childs_vec)
                {
                    this->childs_vec.push_back(it);
                }
                for (auto it : item.first_child_flag_vec)
                {
                    this->first_child_flag_vec.push_back(it);
                }
                for (auto it : item.maximal_repeat_check_vec)
                {
                    this->maximal_repeat_check_vec.push_back(it);
                    this->_stnode_count++;
                }
                for (auto &it : item.edge_char_vec)
                {
                    this->edge_char_vec.push_back(it);
                }
            }
            void computeNextLCPIntervalSet()
            {
                //RINTERVAL intv;
                //RINTERVAL child;

                uint64_t size = this->node_count();
                std::pair<INDEX_SIZE, INDEX_SIZE> output_node;
                std::vector<std::pair<INDEX_SIZE, INDEX_SIZE>> output_children;
                //std::vector<uint8_t> output_chars;
                std::vector<CHAR> output_edge_chars;
                stool::renum::STNodeVector<INDEX_SIZE, CHAR> output_vec;

                for (uint64_t i = 0; i < size; i++)
                {

                    output_children.clear();
                    output_edge_chars.clear();
                    output_vec.clear();

                    this->pop_node(output_node, output_children, output_edge_chars);
                    em->executeWeinerLinkSearch(output_node, output_children, this->store_edge_chars ? &output_edge_chars : nullptr, output_vec);
                    this->import(output_vec);

                }

                this->lcp++;
            }

        public:
            uint64_t increment(uint64_t L, uint64_t &left, uint64_t &right) const
            {
                assert(L > 0);
                assert(this->first_child_flag_vec[L - 1]);
                assert(!this->first_child_flag_vec[L]);

                uint64_t R = L + 1;
                while (R < this->first_child_flag_vec.size() && !this->first_child_flag_vec[R])
                {
                    R++;
                }
                left = this->get_child_left_boundary(L);
                right = this->get_child_right_boundary(R - 1);

                return R + 1;
            }
            uint64_t increment(uint64_t L) const
            {
                assert(L > 0);
                assert(this->first_child_flag_vec[L - 1]);
                assert(!this->first_child_flag_vec[L]);

                uint64_t R = L + 1;
                while (R < this->first_child_flag_vec.size() && !this->first_child_flag_vec[R])
                {
                    R++;
                }
                return R + 1;
            }
            void increment(ITERATOR &iter) const
            {
                this->increment2(iter.child_index, iter.node_index, iter.array_index);
            }
            void increment2(INDEX_SIZE &child_index, INDEX_SIZE &node_index, INDEX_SIZE &array_index) const
            {
                child_index = this->increment(child_index);
                if (child_index >= this->get_integer_array_size())
                {
                    child_index = std::numeric_limits<INDEX_SIZE>::max();
                    node_index = std::numeric_limits<INDEX_SIZE>::max();
                    array_index = std::numeric_limits<INDEX_SIZE>::max();
                }
                else
                {
                    node_index++;
                }
            }
            INDEX_SIZE get_children_count(const ITERATOR &iter) const
            {
                return this->get_children_count(iter.child_index);
            }
            INDEX_SIZE get_children_count(INDEX_SIZE child_index) const
            {
                uint64_t R = this->increment(child_index);
                return R - child_index - 1;
            }
            inline uint64_t get_child_left_boundary(uint64_t child_end_pointer) const
            {
                assert(child_end_pointer > 0);
                assert(child_end_pointer < this->childs_vec.size());
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
                assert(child_end_pointer < this->childs_vec.size());
                return this->childs_vec[child_end_pointer];
            }
            CHAR get_edge_character(uint64_t child_pointer) const
            {
                assert(this->store_edge_chars);
                return this->edge_char_vec[child_pointer];
            }
            INDEX_SIZE get_edge_character(const ITERATOR &iter, uint64_t ith_child) const
            {
                return this->get_edge_character(iter.child_index + ith_child);
            }

            INDEX_SIZE get_child_left_boundary(const ITERATOR &iter, uint64_t ith_child) const
            {
                return this->get_child_left_boundary(iter.child_index + ith_child);
            }
            INDEX_SIZE get_child_right_boundary(const ITERATOR &iter, uint64_t ith_child) const
            {
                return this->get_child_right_boundary(iter.child_index + ith_child);
            }

            INDEX_SIZE get_left(const ITERATOR &iter) const
            {
                return this->get_child_left_boundary(iter.child_index);
            }
            INDEX_SIZE get_left(INDEX_SIZE child_index) const
            {
                return this->get_child_left_boundary(child_index);
            }

            INDEX_SIZE get_right(const ITERATOR &iter) const
            {
                uint64_t left = 0, right = 0;
                this->increment(iter.child_index, left, right);
                return right;
            }
            INDEX_SIZE get_right(INDEX_SIZE child_index) const
            {
                uint64_t left = 0, right = 0;
                this->increment(child_index, left, right);
                return right;
            }

            void set_current_first_iterator(ITERATOR &it) const
            {
                this->set_current_first_iterator(it.child_index, it.node_index, it.array_index);
            }
            void set_current_first_iterator(INDEX_SIZE &child_index, INDEX_SIZE &node_index, INDEX_SIZE &array_index) const
            {
                if (this->is_finished())
                {
                    child_index = std::numeric_limits<INDEX_SIZE>::max();
                    node_index = std::numeric_limits<INDEX_SIZE>::max();
                    array_index = std::numeric_limits<INDEX_SIZE>::max();
                }
                else
                {
                    child_index = 1;
                    node_index = 0;
                    array_index = std::numeric_limits<INDEX_SIZE>::max();
                }
            }

            bool check_maximal_repeat(INDEX_SIZE node_index) const
            {
                return this->maximal_repeat_check_vec[node_index];
            }
            bool check_maximal_repeat(ITERATOR &iter) const
            {
                return this-
                >check_maximal_repeat(iter.node_index);
            }
            */
        };

    } // namespace renum
} // namespace stool