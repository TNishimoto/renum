#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>
#include <queue>
#include "char_interval.hpp"

//#include "sa_lcp.hpp"
using namespace std;

namespace stool
{
    template <typename INDEX, typename CHAR>
    CharInterval<INDEX,CHAR>::CharInterval()
    {
    }
    template <typename INDEX, typename CHAR>
    CharInterval<INDEX,CHAR>::CharInterval(INDEX _i, INDEX _j, CHAR _c) : i(_i), j(_j), c(_c)
    {
    }
    template <typename INDEX, typename CHAR>
    std::string CharInterval<INDEX,CHAR>::to_string() const
    {
        std::string s = "a";
        s[0] = c;
        return "[" + std::to_string(i) + ", " + std::to_string(j) + ", " + s + "]";
    }

} // namespace stool
template class stool::CharInterval<uint32_t, uint8_t>;
template class stool::CharInterval<uint64_t, uint8_t>;
