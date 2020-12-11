#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>

#include "stool/src/debug.hpp"
#include "stool/src/elias_fano_vector.hpp"
#include "../../beller/char_interval.hpp"
#include <sdsl/rmq_support.hpp> // include header for range minimum queries

namespace stool
{
    namespace lcp_on_rlbwt
    {


        template <typename CHAR_VEC, typename INDEX_SIZE>
        class RangeDistinctDataStructure
        {
        private:
            using CHAR = typename CHAR_VEC::value_type;

            const CHAR_VEC *_char_vec;
            const stool::WT *wt;
            std::vector<stool::EliasFanoVector> positionVec;
            //std::vector<INDEX_SIZE> rankVec;
            INDEX_SIZE size;
            sdsl::rmq_succinct_sada<> RMQ;
            sdsl::rmq_succinct_sada<> RmQ;
            std::vector<INDEX_SIZE> tmpRangeDistinctResult;
            std::stack<INDEX_SIZE> tmpSearchStack;

            //std::unordered_map<CHAR, uint64_t> tmpRangeDistinctResult;

            INDEX_SIZE get_next(INDEX_SIZE i)
            {
                CHAR c = (*_char_vec)[i];
                //INDEX_SIZE rank = rankVec[i];
                //INDEX_SIZE rank = wt->rank(i, c);
                INDEX_SIZE rank = wt->rank(i + 1, c);

                if (positionVec[(uint8_t)c].size() == rank + 1)
                {
                    return std::numeric_limits<INDEX_SIZE>::max();
                }
                else
                {
                    return positionVec[(uint8_t)c][rank + 1];
                }
            }
            int64_t get_prev(INDEX_SIZE i)
            {
                CHAR c = (*_char_vec)[i];
                //INDEX_SIZE rank = rankVec[i];
                //INDEX_SIZE rank = wt->rank(i, c);
                INDEX_SIZE rank = wt->rank(i + 1, c);

                //assert(rank == result);

                if (rank == 0)
                {
                    return -1;
                }
                else
                {
                    return positionVec[(uint8_t)c][rank - 1];
                }
            }
            std::vector<INDEX_SIZE> construct_next_vector()
            {
                std::vector<INDEX_SIZE> r;
                r.resize(size, 0);
                for (INDEX_SIZE i = 0; i < size; i++)
                {
                    r[i] = this->get_next(i);
                }
                return r;
            }
            std::vector<INDEX_SIZE> construct_rev_next_vector()
            {
                std::vector<INDEX_SIZE> r;
                r.resize(size, 0);
                for (INDEX_SIZE i = 0; i < size; i++)
                {
                    INDEX_SIZE p = this->get_next(i);
                    if (p == std::numeric_limits<INDEX_SIZE>::max())
                    {
                        r[i] = 0;
                    }
                    else
                    {
                        r[i] = size - p;
                    }
                }
                return r;
            }

            std::vector<int64_t> construct_prev_vector()
            {
                std::vector<int64_t> r;
                r.resize(size, 0);
                for (INDEX_SIZE i = 0; i < size; i++)
                {
                    r[i] = this->get_prev(i);
                }
                return r;
            }

            void search_less(INDEX_SIZE x, INDEX_SIZE i, INDEX_SIZE j, std::stack<INDEX_SIZE> &output)
            {
                std::vector<INDEX_SIZE> r;
                INDEX_SIZE p = RmQ(i, j);
                int64_t value = this->get_prev(p);
                //std::cout << "RMQ[" << i << "," << j << "]=" << p << std::endl;
                if (value >= (int64_t)x)
                {
                    return;
                }
                else
                {
                    output.push(p);
                    if (p > i)
                    {
                        search_less(x, i, p - 1, output);
                    }

                    if (p < j)
                    {
                        search_less(x, p + 1, j, output);
                    }
                }
            }
            void search_than(INDEX_SIZE x, INDEX_SIZE i, INDEX_SIZE j, std::stack<INDEX_SIZE> &output)
            {
                std::vector<INDEX_SIZE> r;
                INDEX_SIZE p = RMQ(i, j);
                INDEX_SIZE value = this->get_next(p);
                //std::cout << "RMQ[" << i << "," << j << "]=" << p << std::endl;
                if (value <= x)
                {
                    return;
                }
                else
                {
                    output.push(p);
                    if (p > i)
                    {
                        search_than(x, i, p - 1, output);
                    }

                    if (p < j)
                    {
                        search_than(x, p + 1, j, output);
                    }
                }
            }

        public:
            RangeDistinctDataStructure()
            {
            }
            void preprocess(const CHAR_VEC *__char_vec, const stool::WT *_wt)
            {
                wt = _wt;
                int32_t charMaxSize = ((int32_t)UINT8_MAX) + 1;
                this->_char_vec = __char_vec;
                this->positionVec.resize(charMaxSize);
                size = _char_vec->size();

                std::vector<std::vector<INDEX_SIZE>> positionSeqVec;
                positionSeqVec.resize(charMaxSize, std::vector<INDEX_SIZE>());
                tmpRangeDistinctResult.resize(charMaxSize, 0);

                //this->rankVec.resize(size, 0);

                for (INDEX_SIZE i = 0; i < size; i++)
                {
                    uint8_t c = (uint8_t)(*_char_vec)[i];
                    //this->rankVec[i] = positionSeqVec[c].size();

                    positionSeqVec[c].push_back(i);
                }

                for (INDEX_SIZE i = 0; i < positionSeqVec.size(); i++)
                {
                    if (positionSeqVec[i].size() > 0)
                    {
                        positionVec[i].construct(&positionSeqVec[i]);
                    }
                }

                auto next_vec = this->construct_rev_next_vector();
                sdsl::rmq_succinct_sada<> next_rmq(&next_vec);
                this->RMQ.swap(next_rmq);

                next_vec.resize(0);
                next_vec.shrink_to_fit();

                auto prev_vec = this->construct_prev_vector();
                sdsl::rmq_succinct_sada<> prev_rmq(&prev_vec);
                this->RmQ.swap(prev_rmq);
            }

            uint64_t range_distinct(INDEX_SIZE i, INDEX_SIZE j, std::vector<CharInterval<INDEX_SIZE>> &output)
            {
                uint64_t count = 0;

                search_less(i, i, j, tmpSearchStack);
                while (tmpSearchStack.size() > 0)
                {
                    INDEX_SIZE p = tmpSearchStack.top();
                    uint8_t c = (uint8_t)(*_char_vec)[p];
                    tmpRangeDistinctResult[c] = p;
                    tmpSearchStack.pop();
                }
                search_than(j, i, j, tmpSearchStack);

                while (tmpSearchStack.size() > 0)
                {
                    INDEX_SIZE p = tmpSearchStack.top();
                    uint8_t c = (uint8_t)(*_char_vec)[p];
                    auto pair = stool::CharInterval<INDEX_SIZE>(tmpRangeDistinctResult[c], p, c);
                    output[count++] = pair;

                    tmpSearchStack.pop();
                }
                return count;
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool