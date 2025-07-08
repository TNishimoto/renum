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
            CharInterval();
            CharInterval(INDEX _i, INDEX _j, CHAR _c);

            std::string to_string() const;
        };
    }
} // namespace stool
