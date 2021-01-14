
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>
#include <queue>
#include "rinterval.hpp"

//#include "sa_lcp.hpp"
using namespace std;

namespace stool
{
    namespace stnode_on_rlbwt
    {

        template <typename INDEX>
        void RInterval<INDEX>::print() const
        {
            std::cout << "[(" << this->beginIndex << ", " << this->beginDiff << "), (" << this->endIndex << ", " << this->endDiff << ")]" << std::endl;
        }
        template <typename INDEX>
        bool RInterval<INDEX>::is_special() const
        {
            return this->beginIndex == std::numeric_limits<INDEX_SIZE>::max();
        }
        template <typename INDEX>
        RInterval<INDEX> RInterval<INDEX>::get_special()
        {
            RInterval r;
            r.beginIndex = std::numeric_limits<INDEX_SIZE>::max();
            return r;
        }
    } // namespace stnode_on_rlbwt

} // namespace stool
template class stool::stnode_on_rlbwt::RInterval<uint32_t>;
template class stool::stnode_on_rlbwt::RInterval<uint64_t>;