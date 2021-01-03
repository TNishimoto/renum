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
    namespace lcp_on_rlbwt
    {

        template <typename TRAVERSER>
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
            /*
                INDEX_SIZE get_lcp() const {
                    return traverser->get_current_lcp();
                }
                INDEX_SIZE get_left() const {
                    if(this->divided_array){
                        return traverser->get_child_start_position(this->child_index);
                    }else{
                        return traverser->get_child_start_position(this->child_index);
                    }
                }
                INDEX_SIZE get_right() const {
                    uint64_t left = 0, right = 0;
                    traverser->increment(this->child_index, left, right);
                    return right;
                }
                */
            /*
            bool isEnd() const
            {
                return this->node_index == std::numeric_limits<INDEX_SIZE>::max();
            }
            */

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
        };
        template <typename TRAVERSER>
        class STDepthIterator
        {
            TRAVERSER *traverser;

        public:
            using INDEX_SIZE = typename TRAVERSER::index_type;
            INDEX_SIZE depth;

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

            STNodeIterator<TRAVERSER> operator*() const
            {
                return STNodeIterator<TRAVERSER>(traverser, true);
            }

            STDepthIterator &operator++()
            {
                if (traverser->succ())
                {
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
                return this->depth != rhs.depth;
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool