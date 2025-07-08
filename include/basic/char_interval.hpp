#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>
#include <queue>

// #include "sa_lcp.hpp"
// using namespace std;
// using namespace sdsl;

namespace stool
{
    namespace renum
    {

        // using WT = sdsl::wt_huff<>;
        template <typename INDEX, typename CHAR>
        class CharInterval
        {
        public:
            INDEX i;
            INDEX j;
            CHAR c;
            CharInterval(){

            }
            CharInterval(INDEX _i, INDEX _j, CHAR _c): i(_i), j(_j), c(_c){

            }

            std::string to_string() const
        {
            std::string s = "a";
            s[0] = c;
            return "[" + std::to_string(i) + ", " + std::to_string(j) + ", " + s + "]";
        }
        };
    }
} // namespace stool
