
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>
#include <queue>
#include <iostream>
#include <fstream>

#include "../../include/basic/rinterval.hpp"

//#include "sa_lcp.hpp"
using namespace std;

namespace stool
{
    namespace renum
    {

        template <typename INDEX>
        void RInterval<INDEX>::print() const
        {
            std::cout << "[(" << this->beginIndex << ", " << this->beginDiff << "), (" << this->endIndex << ", " << this->endDiff << ")]" << std::endl;
        }
        template <typename INDEX>
        bool RInterval<INDEX>::is_special() const
        {
            return this->beginIndex == std::numeric_limits<INDEX>::max();
        }
        template <typename INDEX>
        RInterval<INDEX> RInterval<INDEX>::get_special()
        {
            RInterval r;
            r.beginIndex = std::numeric_limits<INDEX>::max();
            return r;
        }
    } // namespace renum

} // namespace stool
template class stool::renum::RInterval<uint32_t>;
template class stool::renum::RInterval<uint64_t>;