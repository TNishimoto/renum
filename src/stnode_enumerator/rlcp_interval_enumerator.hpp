#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>
#include "application.hpp"
/*
namespace stool
{
    namespace lcp_on_rlbwt
    {
        template <typename RLBWT_STR, typename INDEX_SIZE, typename RLBWTDS, typename FPOSDS>
        class RLCPIntervalEnumerator
        {
            using CHAR = typename RLBWT_STR::char_type;
            using CHARVEC = typename RLBWT_STR::char_vec_type;
            using UCHAR = typename std::make_unsigned<CHAR>::type;
            using RINTERVAL = RInterval<INDEX_SIZE>;

            NextRIntervalStorageConstructor<RLBWT_STR, INDEX_SIZE, RLBWTDS> wds;
            //RLBWTDataStructures<RLBWT_STR, INDEX_SIZE, FPOSDS> *_RLBWTDS;
            RIntervalStorage<INDEX_SIZE> rintervalStorage;
            RIntervalStorage<INDEX_SIZE> rintervalTmpStorage;
            uint64_t current_lcp = 0;

        public:
            void create_next_lcp_intervals()
            {
                if (current_lcp == 0)
                {
                    this->wds.computeFirstLCPIntervalSet(this->rintervalTmpStorage);
                    this->rintervalStorage.swap(this->rintervalTmpStorage);
                }
                else
                {
                    this->wds.computeNextLCPIntervalSet(this->rintervalStorage, this->rintervalTmpStorage);
                    this->rintervalStorage.swap(this->rintervalTmpStorage);
                }
                current_lcp++;
            }
            class iterator
            {
            private:
                RLCPIntervalEnumerator *parent;
                uint64_t _index = 0;
                uint64_t _lcp = 0;

            public:
                iterator() = default;
                iterator(RLCPIntervalEnumerator *__parent, bool isEnd) : parent(__parent)
                {
                    if (isEnd)
                    {
                        this->_index = UINT64_MAX;
                        this->_lcp = UINT64_MAX;
                    }
                }

            public:
                iterator &operator++()
                {
                    if (this->_index < this->parent->rintervalStorage.lcpIntvCount)
                    {
                        this->_index++;
                    }
                    else
                    {
                        this->parent->create_next_lcp_intervals();
                        if (this->parent->rintervalStorage.lcpIntvCount == 0)
                        {
                            this->_lcp = UINT64_MAX;
                            this->_index = UINT64_MAX;
                        }
                        else
                        {
                            this->_lcp++;
                            this->_index = 0;
                        }
                    }

                    return *this;
                }
                std::pair<RINTERVAL, bool> operator*()
                {
                    bool b = this->_index == 0;
                    return std::pair<RINTERVAL, bool>(this->parent->rintervalStorage.lcpIntvVec[this->_index], b);
                }
                bool operator!=(const iterator &rhs)
                {
                    return (_lcp != rhs._lcp) || (_index != rhs._index);
                }
            };
            iterator begin()
            {
                this->current_lcp = 0;
                this->wds.clear();
                auto it = iterator(this, false);
                return it;
            }
            iterator end()
            {
                auto it = iterator(this, true);
                return it;
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool
*/
