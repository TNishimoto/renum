#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <deque>
#include <vector>
#include <cassert>

#include <type_traits>

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace stnode_on_rlbwt
    {

        template <typename TRAVERSER, typename CHAR = uint8_t>
        class STNodeIterator
        {
            const TRAVERSER *traverser;

        public:
            using INDEX_SIZE = typename TRAVERSER::index_type;
            INDEX_SIZE node_index;
            INDEX_SIZE child_index;
            INDEX_SIZE array_index;

            STNodeIterator() = default;

            STNodeIterator(const TRAVERSER *_traverser, bool is_start) : traverser(_traverser)
            {

                if (is_start)
                {
                    traverser->set_current_first_iterator(*this);
                }
                else
                {
                    this->node_index = std::numeric_limits<INDEX_SIZE>::max();
                    this->child_index = std::numeric_limits<INDEX_SIZE>::max();
                    this->array_index = std::numeric_limits<INDEX_SIZE>::max();

                    //this->lcp = std::numeric_limits<INDEX_SIZE>::max();
                }
            }

            INDEX_SIZE get_children_count() const
            {
                return traverser->get_children_count(*this);
            }
            INDEX_SIZE get_left() const
            {
                return traverser->get_left(*this);
            }
            INDEX_SIZE get_right() const
            {
                return traverser->get_right(*this);
            }
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

            STNodeIterator &operator++()
            {
                traverser->increment(*this);
                return *this;
            }

            STNodeIterator &operator++(int)
            {
                ++(*this);
                return *this;
            }

            bool operator!=(const STNodeIterator &rhs) const
            {
                return this->node_index != rhs.node_index || this->child_index != rhs.child_index;
            }
            void print() const
            {
                std::cout << "[" << this->get_left() << ", " << this->get_right() << "]" << std::flush;
                uint64_t w = this->get_children_count();
                for (uint64_t i = 0; i < w; i++)
                {
                    uint64_t left = this->get_child_left_boundary(i);
                    uint64_t right = this->get_child_right_boundary(i);
                    if (traverser->has_edge_characters())
                    {
                        uint8_t c = this->get_edge_character(i);
                        std::string charStr = c == 0 ? "$(0)" : std::string(1, c);
                        std::cout << "(" << left << ", " << right << ", " << charStr << ")" << std::flush;
                    }
                    else
                    {
                        std::cout << "(" << left << ", " << right << ")" << std::flush;
                    }
                }
            }
            bool has_edge_characters() const
            {
                return this->traverser->has_edge_characters();
            }
        };
        template <typename TRAVERSER>
        class STDepthIterator
        {
        private:
            using INDEX_SIZE = typename TRAVERSER::index_type;

            TRAVERSER *traverser;
            INDEX_SIZE depth;

        public:
            STDepthIterator() = default;

            STDepthIterator(TRAVERSER *_traverser, bool is_start) : traverser(_traverser)
            {

                if (is_start)
                {
                    this->depth = 0;
                }
                else
                {
                    this->depth = std::numeric_limits<INDEX_SIZE>::max();
                }
            }
            uint64_t get_depth() const
            {
                return this->depth;
            }
            uint64_t child_count() const
            {
                return this->traverser->child_count();
            }
            uint64_t node_count() const
            {
                return this->traverser->node_count();
            }

            STNodeIterator<TRAVERSER> operator*() const
            {
                return STNodeIterator<TRAVERSER>(traverser, true);
            }

            STDepthIterator &operator++()
            {

                if (traverser->succ())
                {

                    this->depth++;
                    return *this;
                }
                else
                {
                    this->depth = std::numeric_limits<INDEX_SIZE>::max();
                    return *this;
                }
            }
            STDepthIterator &operator++(int)
            {
                ++(*this);
                return *this;
            }

            STNodeIterator<TRAVERSER> begin() const
            {
                return STNodeIterator<TRAVERSER>(traverser, true);
            }
            STNodeIterator<TRAVERSER> end() const
            {
                return STNodeIterator<TRAVERSER>(traverser, false);
            }
            bool operator!=(const STDepthIterator &rhs) const
            {
                return this->get_depth() != rhs.get_depth();
            }
            void print(){
                std::cout << "\033[31m [BEGIN]Nodes: " << std::endl;
                for(auto node_it = this->begin();node_it != this->end();++node_it){
                    std::cout << "   ";
                    node_it.print();
                    std::cout << std::endl;
                }
                std::cout << "[END] \033[39m" << std::endl;

            }
        };

        class STDepthIteratorErrorChecker
        {
            public:
            template <typename TRAVERSER>
            static void error_check(TRAVERSER &it)
            {
                it.print();
                /*
                std::cout << "Error Check" << std::endl;
                it.print();
                uint64_t w = it.get_children_count();

                uint64_t prevR = 0;

                for (uint64_t i = 0; i < w; i++)
                {
                    uint64_t L = it.get_child_left_boundary(i);
                    uint64_t R = it.get_child_right_boundary(i);

                    std::cout << "[R/L/prevR]= " << R << "/" << L << "/" << prevR << std::endl;
                    if (i != 0)
                    {
                        assert(L - 1 == prevR);
                    }
                    prevR = R;
                }
                */
            }
        };

    } // namespace stnode_on_rlbwt
} // namespace stool