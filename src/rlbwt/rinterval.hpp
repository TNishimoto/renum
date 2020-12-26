#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>

namespace stool
{
    namespace lcp_on_rlbwt
    {
        template <typename INDEX_SIZE>
        struct RInterval
        {
            INDEX_SIZE beginIndex;
            INDEX_SIZE beginDiff;
            INDEX_SIZE endIndex;
            INDEX_SIZE endDiff;

            void print() const
            {
                std::cout << "[(" << this->beginIndex << ", " << this->beginDiff << "), (" << this->endIndex << ", " << this->endDiff << ")]" << std::endl;
            }
            template <typename FPOSDS>
            void print2(const FPOSDS &fposArray) const
            {
                if (this->is_special())
                {
                    std::cout << "[BOTTOM]" << std::endl;
                }
                else
                {
                    assert(this->beginIndex < fposArray.size());
                    assert(this->endIndex < fposArray.size());
                    //stool::Printer::print(fposArray);

                    INDEX_SIZE begin_pos = fposArray[this->beginIndex] + this->beginDiff;
                    INDEX_SIZE end_pos = fposArray[this->endIndex] + this->endDiff;

                    std::cout << "[(" << this->beginIndex << ", " << this->beginDiff << "), (" << this->endIndex << ", " << this->endDiff << ")]" << std::flush;
                    std::cout << "[" << begin_pos << ", " << end_pos << "]" << std::flush;
                    assert(begin_pos <= end_pos);
                }
            }
            template <typename FPOSDS>
            stool::LCPInterval<uint64_t> get_lcp_interval(uint64_t lcp,const FPOSDS &fposArray) const
            {
                    assert(this->beginIndex < fposArray.size());
                    assert(this->endIndex < fposArray.size());

                    INDEX_SIZE begin_pos = fposArray[this->beginIndex] + this->beginDiff;
                    INDEX_SIZE end_pos = fposArray[this->endIndex] + this->endDiff;
                    return stool::LCPInterval<uint64_t>(begin_pos, end_pos,lcp);
            }
            /*
            template <typename FPOSDS>
            void to_lcp_interval(uint64_t lcp,const FPOSDS &fposArray, stool::LCPInterval<uint64_t> &output) const
            {
                    assert(this->beginIndex < fposArray.size());
                    assert(this->endIndex < fposArray.size());

                    output.i = fposArray[this->beginIndex] + this->beginDiff;
                    output.j = fposArray[this->endIndex] + this->endDiff;
                    output.lcp = lcp;
            }
            */

            bool is_special() const
            {
                return this->beginIndex == std::numeric_limits<INDEX_SIZE>::max();
            }

            static RInterval get_special()
            {
                RInterval r;
                r.beginIndex = std::numeric_limits<INDEX_SIZE>::max();
                return r;
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool