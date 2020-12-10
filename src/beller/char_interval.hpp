#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wt_gmr.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>

#include <queue>
#include "../module/stool/src/sa_bwt_lcp.hpp"
#include "fmindex.hpp"

//#include "sa_lcp.hpp"
using namespace std;
using namespace sdsl;

namespace stool
{
    template <typename index_type>
    class CharInterval
    {
    public:
        index_type i;
        index_type j;
        uint8_t c;
        CharInterval()
        {
        }
        CharInterval(index_type _i, index_type _j, uint8_t _c) : i(_i), j(_j), c(_c)
        {
        }

        std::string to_string() const
        {
            std::string s = "a";
            s[0] = c;
            return "[" + std::to_string(i) + ", " + std::to_string(j) + ", " + s + "]";
        }
        /*
        static uint64_t LF(uint64_t i, int_vector<> &bwt, std::vector<uint64_t> &C, wt_gmr<> &wt)
        {
            uint8_t c = bwt[i];
            uint64_t cNum = wt.rank(i, c);
            return C[c] + cNum;
        }

        static std::vector<CharInterval> getIntervals_naive(uint64_t i, uint64_t j, int_vector<> &bwt, std::vector<uint64_t> &C, wt_gmr<> &wt)
        {

            std::vector<CharInterval> r;
            std::map<uint8_t, std::pair<uint64_t, uint64_t>> mapper;
            for (uint64_t x = i; x <= j; x++)
            {
                uint8_t c = bwt[x];
                auto f = mapper.find(c);
                if (f == mapper.end())
                {
                    mapper[c] = std::pair<uint64_t, uint64_t>(x, x);
                }
                else
                {
                    mapper[c] = std::pair<uint64_t, uint64_t>((*f).second.first, x);
                }
            }
            for (auto it : mapper)
            {
                uint8_t c = it.first;
                uint64_t l_intv = it.second.first;
                uint64_t r_intv = it.second.second;
                //std::cout << "LF: " << l_intv << ", " << r_intv << std::endl;
                uint64_t new_l = LF(l_intv, bwt, C, wt);
                uint64_t new_r = LF(r_intv, bwt, C, wt);
                CharInterval ci;
                ci.i = new_l;
                ci.j = new_r;
                ci.c = c;
                r.push_back(ci);
            }
            return r;
        }
        */
    };
    class IntervalSearchDataStructure
    {
        public:
        std::vector<uint64_t> C;
        sdsl::wt_huff<> wt;
        uint8_t lastChar;
        std::vector<uint8_t> cs;
        std::vector<uint64_t> cs1;
        std::vector<uint64_t> cs2;

        void initialize(int_vector<> &bwt)
        {
            this->lastChar = bwt[bwt.size() - 1];
            stool::FMIndex::constructC(bwt, C);
            construct_im(wt, bwt);

            cs.resize(256, 0);
            cs1.resize(256, 0);
            cs2.resize(256, 0);
        }
        template <typename INDEX_SIZE>
        uint64_t getIntervals(INDEX_SIZE i, INDEX_SIZE j, std::vector<CharInterval<INDEX_SIZE>> &output)
        {
            using CHARINTV = CharInterval<INDEX_SIZE>;
            std::vector<CHARINTV> r;
            uint64_t k;
            uint64_t newJ = j + 1 == wt.size() ? wt.size() : j + 2;
            uint64_t p =0;


            //std::cout << "@[" << i << "/" << j << ", " << wt.size() << "]" << std::endl;

            sdsl::interval_symbols(wt, i + 1, newJ, k, cs, cs1, cs2);

            bool b = j + 1 < wt.size();
            for (INDEX_SIZE x = 0; x < k; x++)
            {
                INDEX_SIZE left = C[cs[x]] + cs1[x];
                INDEX_SIZE right = left + (cs2[x] - cs1[x] - 1);

                //uint64_t right = C[cs[x]] + cs2[x]+1;

                if (j + 1 == wt.size() && cs[x] == lastChar)
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
                INDEX_SIZE num = wt.rank(wt.size(), lastChar) + 1;
                INDEX_SIZE left = C[lastChar] + num - 1;
                INDEX_SIZE right = left;
                output[p++] = CHARINTV(left, right, lastChar);
                //r.push_back(CHARINTV(left, right, lastChar));
            }

            return p;
        }
    };

} // namespace stool