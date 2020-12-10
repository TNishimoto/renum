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

//#include "sa_lcp.hpp"
using namespace std;
using namespace sdsl;

namespace stool
{

    class FMIndex
    {
    public:
        static uint64_t LF(uint64_t i, int_vector<> &bwt, std::vector<uint64_t> &C, wt_gmr<> &wt)
        {
            uint8_t c = bwt[i];
            uint64_t cNum = wt.rank(i, c);
            return C[c] + cNum;
        }

        template <typename TEXT, typename OUTPUT>
        static void constructC(TEXT &text, OUTPUT &output)
        {
            uint64_t CHARMAX = UINT8_MAX+1;
        
            vector<uint64_t> tmp;

            tmp.resize(CHARMAX, 0);
            output.resize(CHARMAX, 0);

            for (uint64_t i = 0; i < text.size(); i++)
            {
                assert(text[i] >= 0 && text[i] < tmp.size());
                tmp[text[i]]++;
            }

            for (uint64_t i = 1; i < output.size(); i++)
            {
                output[i] = output[i - 1] + tmp[i - 1];
            }
        }
        static void constructSelect(int_vector<> &bwt, wt_gmr<> &wt)
        {
            construct_im(wt, bwt);
        }
        /*
        static void constructBWT(string &bwt, int_vector<> &outputBWT)
        {
            outputBWT.width(8);
            outputBWT.resize(bwt.size());
            for (uint64_t i = 0; i < bwt.size(); i++)
            {
                outputBWT[i] = bwt[i];
            }
        }
        static void constructBWT(std::vector<char> &bwt, int_vector<> &outputBWT)
        {
            outputBWT.width(8);
            outputBWT.resize(bwt.size());
            for (uint64_t i = 0; i < bwt.size(); i++)
            {
                outputBWT[i] = bwt[i];
            }
        }
        */
        template <typename INDEX = uint64_t>
        static void constructBWT(const std::vector<uint8_t> &text, const std::vector<INDEX> &sa, sdsl::int_vector<> &outputBWT)
        {
            outputBWT.width(8);
            outputBWT.resize(text.size());
            
            INDEX n = text.size();
            for (INDEX i = 0; i < n; i++)
            {
                if (sa[i] == 0)
                {
                    outputBWT[i] = (uint8_t)text[n - 1];
                }
                else
                {
                    outputBWT[i] = (uint8_t)text[sa[i] - 1];
                }
            }
        }

        static uint64_t FL(uint64_t i, std::vector<uint64_t> &C, wt_gmr<> &wt)
        {
            uint64_t x = 0;
            for (x = 0; x < C.size(); x++)
            {
                //std::cout << "xxx:" << x << "/" << C[x] << std::endl;
                if (C[x] <= i && i < C[x + 1])
                {
                    break;
                }
            }
            uint8_t c = x;
            uint64_t nth = i - C[x];

            //std::cout << "xxx:" << c << "/" << nth << std::endl;
            return wt.select((nth + 1), c);
        }
        static uint64_t get_start_pos(int_vector<> &bwt)
        {
            for (uint64_t i = 0; i < bwt.size(); i++)
            {
                if (bwt[i] == 0)
                {
                    return i;
                }
            }
            return UINT64_MAX;
        }
        static uint8_t getUpperChar(uint64_t i, std::vector<uint64_t> &C)
        {
            std::vector<uint8_t> chars;
            for (uint64_t x = 0; x < C.size(); x++)
            {
                for (uint64_t p = C[x]; p < C[x + 1]; p++)
                {
                    chars.push_back(x);
                }
            }
            if (i == 0)
            {
                return chars[chars.size() - 1];
            }
            else
            {
                return chars[i - 1];
            }
        }
        static uint8_t getFChar(uint64_t i, std::vector<uint64_t> &C)
        {
            std::vector<uint8_t> chars;
            for (uint64_t x = 0; x < C.size(); x++)
            {
                for (uint64_t p = C[x]; p < C[x + 1]; p++)
                {
                    chars.push_back(x);
                }
            }
            return chars[i];
        }
    };

} // namespace stool