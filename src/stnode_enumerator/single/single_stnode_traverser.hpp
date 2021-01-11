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
#include "bit_deque.hpp"
#include "../node_iterator.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        template <typename INDEX_SIZE, typename INTERVAL_SEARCH, typename CHAR = uint8_t>
        class SingleSTNodeTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using ITERATOR = STNodeIterator<SingleSTNodeTraverser>;
            using DEPTH_ITERATOR = STDepthIterator<SingleSTNodeTraverser>;

            uint64_t _stnode_count = 0;
            std::deque<INDEX_SIZE> childs_vec;
            std::deque<CHAR> edge_char_vec;

            BitDeque first_child_flag_vec;
            BitDeque maximal_repeat_check_vec;
            bool store_edge_chars = false;
            //std::deque<bool> first_child_flag_vec;
            //std::deque<bool> maximal_repeat_check_vec;
            int64_t lcp = -1;
            INTERVAL_SEARCH *em;

        public:
            using index_type = INDEX_SIZE;
            
            INTERVAL_SEARCH* get_interval_search_deta_structure() const {
                return this->em;
            }

            bool has_edge_characters() const
            {
                return this->store_edge_chars;
            }
            SingleSTNodeTraverser()
            {
            }
            void initialize(INTERVAL_SEARCH *_interval_search, bool _store_edge_chars)
            {
                this->em = _interval_search;
                this->store_edge_chars = _store_edge_chars;
                this->lcp = -1;
                /*
                this->em = new ExplicitWeinerLinkEmulator<RLBWTDS>();
                this->em->initialize(__RLBWTDS);
                */
            }

            bool is_finished() const
            {
                return this->get_current_lcp() >= 0 && this->node_count() == 0;
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
            void first_compute()
            {
                this->first_child_flag_vec.clear();
                this->edge_char_vec.clear();

                //em.computeFirstLCPIntervalSet();
                std::vector<CharInterval<INDEX_SIZE>> r = em->getFirstChildren();
                for (uint64_t i = 0; i < r.size(); i++)
                {
                    auto &it = r[i];
                    if (i == 0)
                    {
                        this->childs_vec.push_back(it.i);
                        this->first_child_flag_vec.push_back(true);
                        if (store_edge_chars)
                        {
                            this->edge_char_vec.push_back(0);
                        }
                    }
                    this->childs_vec.push_back(it.j);
                    this->first_child_flag_vec.push_back(false);
                    if (store_edge_chars)
                    {
                        this->edge_char_vec.push_back(it.c);
                    }
                }
                this->maximal_repeat_check_vec.push_back(true);
                assert(this->first_child_flag_vec[0]);
                this->_stnode_count++;
                this->lcp = 0;
            }

            void pop_node()
            {
                assert(this->first_child_flag_vec[0]);

                uint64_t L = 1;
                uint64_t _left = 0, _right = 0;
                uint64_t nextL = this->increment(L, _left, _right);
                //output_node.first = _left;
                //output_node.second = _right;
                uint64_t _count = nextL - L - 1;
                for (uint64_t i = 0; i < _count; i++)
                {
                    //uint64_t left = this->get_child_left_boundary(L);
                    //uint64_t right = this->get_child_right_boundary(L);
                    if (this->store_edge_chars)
                    {
                        //CHAR c = this->get_edge_character(L);
                        //output_edge_chars.push_back(c);
                        this->edge_char_vec.pop_front();
                    }
                    //output_children.push_back(std::pair<INDEX_SIZE, INDEX_SIZE>(left, right));

                    //assert(left <= right);
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
            void import(stool::lcp_on_rlbwt::STNodeVector<INDEX_SIZE, CHAR> &item){
                for(auto &it : item.childs_vec){
                    this->childs_vec.push_back(it);
                }
                for(auto it : item.first_child_flag_vec){
                    this->first_child_flag_vec.push_back(it);
                }
                for(auto it : item.maximal_repeat_check_vec){
                    this->maximal_repeat_check_vec.push_back(it);
                                    this->_stnode_count++;

                }
                for(auto &it : item.edge_char_vec){
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
                std::vector<uint8_t> output_chars;
                std::vector<CHAR> output_edge_chars;

                stool::lcp_on_rlbwt::STNodeVector<INDEX_SIZE, CHAR> tmp;
                auto beg = ITERATOR(this, true);
                //auto beg = this->begin().begin();
                for (uint64_t i = 0; i < size; i++)
                {
                    tmp.clear();
                    WeinerLinkCommonFunctions::compute_weiner_links(*this->em, beg, tmp);

                    this->import(tmp);
                    this->pop_node();

                    /*
                    output_children.clear();
                    output_chars.clear();
                    output_edge_chars.clear();

                    if (this->store_edge_chars)
                    {
                        em->executeWeinerLinkSearch(output_node, output_children, &output_edge_chars, output_chars);
                    }
                    else
                    {
                        em->executeWeinerLinkSearch(output_node, output_children, nullptr, output_chars);
                    }

                    for (auto &c : output_chars)
                    {
                        this->add(c, *em);
                    }
                    */
                }

                this->lcp++;
            }

            void add(uint8_t c, INTERVAL_SEARCH &em)
            {
                RINTERVAL copy;
                uint64_t st_left = UINT64_MAX;
                uint64_t st_right = 0;
                std::pair<INDEX_SIZE, INDEX_SIZE> child;

                uint64_t count = em.get_width(c);

                for (uint64_t j = 0; j < count; j++)
                {
                    CHAR edge_char = em.get_child(c, j, child);
                    uint64_t left = child.first;
                    uint64_t right = child.second;

                    if (left < st_left)
                    {
                        st_left = left;
                    }
                    if (right > st_right)
                    {
                        st_right = right;
                    }
                    if (j == 0)
                    {
                        this->childs_vec.push_back(left);
                        this->first_child_flag_vec.push_back(true);
                        if (store_edge_chars)
                        {
                            this->edge_char_vec.push_back(0);
                        }
                    }
                    this->childs_vec.push_back(right);
                    this->first_child_flag_vec.push_back(false);
                    if (store_edge_chars)
                    {
                        this->edge_char_vec.push_back(edge_char);
                    }
                }
                assert(this->first_child_flag_vec[0]);
                bool isMaximalRepeat = this->em->checkMaximalRepeat(st_left, st_right);
                this->maximal_repeat_check_vec.push_back(isMaximalRepeat);
                this->_stnode_count++;
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
                return this->check_maximal_repeat(iter.node_index);
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool