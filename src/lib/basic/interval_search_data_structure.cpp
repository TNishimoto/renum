#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>
#include "../../include/basic/interval_search_data_structure.hpp"

template <typename CHAR>
void stool::IntervalSearchDataStructure<CHAR>::initialize(stool::WT *_wt, std::vector<uint64_t> *_C, CHAR _lastChar)
{
    this->wt = _wt;
    this->C = _C;
    this->lastChar = _lastChar;
    /*
            stool::FMIndex::constructC(bwt, C);
            construct_im(wt, bwt);
            */

    cs.resize(256, 0);
    cs1.resize(256, 0);
    cs2.resize(256, 0);
}
template <typename CHAR>
template <typename INDEX_SIZE>
uint64_t stool::IntervalSearchDataStructure<CHAR>::getIntervals(INDEX_SIZE i, INDEX_SIZE j, std::vector<CharInterval<INDEX_SIZE, CHAR>> &output)
{
    using CHARINTV = CharInterval<INDEX_SIZE, CHAR>;
    std::vector<CHARINTV> r;
    uint64_t k;
    uint64_t newJ = j + 1 == wt->size() ? wt->size() : j + 2;
    uint64_t p = 0;

    //std::cout << "@[" << i << "/" << j << ", " << wt->size() << "]" << std::endl;

    sdsl::interval_symbols(*wt, i + 1, newJ, k, cs, cs1, cs2);

    bool b = j + 1 < wt->size();
    for (INDEX_SIZE x = 0; x < k; x++)
    {
        INDEX_SIZE left = (*C)[cs[x]] + cs1[x];
        INDEX_SIZE right = left + (cs2[x] - cs1[x] - 1);

        //uint64_t right = C[cs[x]] + cs2[x]+1;

        if (j + 1 == wt->size() && cs[x] == lastChar)
        {
            right++;
            b = true;
        }

        //std::cout << ((int)cs[x]) << "/" << left << "/" << right << std::endl;
        output[p++] = CHARINTV(left, right, cs[x]);
        //r.push_back(CHARINTV(left, right, cs[x]));
    }
    if (!b)
    {
        INDEX_SIZE num = wt->rank(wt->size(), lastChar) + 1;
        INDEX_SIZE left = (*C)[lastChar] + num - 1;
        INDEX_SIZE right = left;
        output[p++] = CHARINTV(left, right, lastChar);
        //r.push_back(CHARINTV(left, right, lastChar));
    }

    return p;
}

template class stool::IntervalSearchDataStructure<uint8_t>;

template uint64_t stool::IntervalSearchDataStructure<uint8_t>::getIntervals(uint32_t, uint32_t, std::vector<CharInterval<uint32_t, uint8_t>> &);
template uint64_t stool::IntervalSearchDataStructure<uint8_t>::getIntervals(uint64_t, uint64_t, std::vector<CharInterval<uint64_t, uint8_t>> &);
