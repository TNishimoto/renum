#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>

#include "../basic/interval_search_data_structure.hpp"
#include "../rlbwt/range_distinct/succinct_range_distinct.hpp"
#include "beller_component.hpp"

//#include "sa_lcp.hpp"
//using namespace std;
//using namespace sdsl;

namespace stool
{
    namespace beller
    {

        template <typename INDEX>
        std::vector<stool::LCPInterval<INDEX>> computeLCPIntervals(IntervalSearchDataStructure<uint8_t> &range)
        {
            stool::beller::BellerComponent<INDEX> comp;
            comp.initialize(range);
            bool isEnd = false;

            stool::beller::OutputStructure<INDEX> os;
            os.is_output_maximal_substrings = false;
            //std::vector<stool::LCPInterval<INDEX>> r;

            while (!isEnd)
            {
                comp.computeLCPIntervals(range, isEnd, os);
                /*
                for (auto it : r2)
                {
                    r.push_back(it);
                }
                */
            }
            return os.outputIntervals;
        }
        template <typename INDEX>
        uint64_t outputMaximalSubstrings(IntervalSearchDataStructure<uint8_t> &range, stool::beller::OutputStructure<INDEX> &os)
        {
            stool::beller::BellerComponent<INDEX> comp;
            comp.initialize(range);

            uint64_t print_ratio = 100;
            uint64_t print_counter = 0;
            uint64_t strSize = range.wt->size();
            uint64_t pvalue = strSize / print_ratio;

            bool isEnd = false;
            while (!isEnd)
            {
                uint64_t pp_num = pvalue * print_counter;

                if (os.lcp_interval_count >= pp_num)
                {
                    std::cout << "[" << (print_counter) << "/" << print_ratio << "] ";
                    std::cout << "LCP = " << comp.lcp;
                    std::cout << ", Peak count = " << os.peak;

                    std::cout << ", Process = [" << os.lcp_interval_count << "/" << strSize << "]";

                    std::cout << std::endl;

                    print_counter++;
                }

                comp.computeLCPIntervals(range, isEnd, os);
            }

            std::cout << "Range Distinct Count = " << comp.debugCounter << "/" << range.wt->size() << std::endl;
            return os.count;
        }

    } // namespace beller
} // namespace stool