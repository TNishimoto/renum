#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>

#include "char_interval.hpp"

//#include "sa_lcp.hpp"
using namespace std;
using namespace sdsl;

namespace stool
{
    namespace beller
    {

        template <typename INDEX_SIZE>
        bool check(std::vector<CharInterval<INDEX_SIZE>> &vec1, std::vector<CharInterval<INDEX_SIZE>> &vec2)
        {
            std::sort(vec1.begin(), vec1.end(), [](const CharInterval<INDEX_SIZE> &lhs, const CharInterval<INDEX_SIZE> &rhs) {
                return lhs.c < rhs.c;
            });
            std::sort(vec2.begin(), vec2.end(), [](const CharInterval<INDEX_SIZE> &lhs, const CharInterval<INDEX_SIZE> &rhs) {
                return lhs.c < rhs.c;
            });

            if (vec1.size() != vec2.size())
            {
                /*
                std::cout << "//" << tmp.size() << ", " << tmpx.size() << std::endl;
                for (auto intv : tmp)
                {
                    std::cout << intv.to_string();
                }
                */
                        std::cout << "Distinct Sizes" << std::endl;

                throw -1;
            }
            else
            {
                for (uint64_t i = 0; i < vec1.size(); i++)
                {
                    if (vec1[i].i != vec2[i].i || vec1[i].j != vec2[i].j || vec1[i].c != vec2[i].c)
                    {
                        std::cout << "Distinct Values" << std::endl;
                        std::cout << vec1[i].to_string() << "/" << vec2[i].to_string() << std::endl;
                        throw -1;
                    }
                }
                return true;
            }
        }

        template <typename INDEX = uint64_t>
        std::vector<LCPInterval<INDEX>> naive_compute_lcp_intervals(const std::vector<INDEX> &sa, const std::vector<INDEX> &lcpArray)
        {
            std::vector<LCPInterval<INDEX>> r;
            for (uint64_t i = 0; i < sa.size(); i++)
            {
                uint64_t limit_lcp = i == 0 ? 0 : lcpArray[i];
                uint64_t current_lcp = sa.size() - sa[i];
                for (uint64_t x = i + 1; x <= sa.size(); x++)
                {
                    uint64_t lcp = x == sa.size() ? 0 : lcpArray[x];
                    //std::cout << i << "/" << x << "/" << current_lcp << "/"<< lcp << "/" << limit_lcp<< std::endl;

                    if (current_lcp > lcp)
                    {
                        r.push_back(LCPInterval<INDEX>(i, x - 1, current_lcp));
                        current_lcp = lcp;
                    }

                    if (current_lcp <= limit_lcp)
                    {
                        break;
                    }
                }
            }
            r.push_back(LCPInterval<INDEX>(0, sa.size() - 1, 0));
            std::sort(
                r.begin(),
                r.end(),
                LCPIntervalPreorderComp<INDEX>());

            return r;
        }

        template <typename INDEX = uint64_t>
        std::vector<LCPInterval<INDEX>> naive_compute_complete_lcp_intervals(const std::vector<INDEX> &sa, const std::vector<INDEX> &lcpArray)
        {
            std::vector<LCPInterval<INDEX>> r = naive_compute_lcp_intervals(sa, lcpArray);

            std::vector<LCPInterval<INDEX>> correct_lcp_intervals;
            for (auto it : r)
            {
                if (it.j - it.i != 0)
                {
                    correct_lcp_intervals.push_back(it);
                }
                else
                {
                    it.lcp = std::numeric_limits<INDEX>::max() - 1;
                    //correct_lcp_intervals.push_back(it);
                }
            }

            return correct_lcp_intervals;
        }

        template <typename INDEX = uint64_t>
        bool equal_check_lcp_intervals(std::vector<LCPInterval<INDEX>> &item1, std::vector<LCPInterval<INDEX>> &item2)
        {
            using LCPINTV = LCPInterval<INDEX>;

            sort(item1.begin(), item1.end(), [](const LCPINTV &lhs, const LCPINTV &rhs) {
                if (lhs.lcp != rhs.lcp)
                {
                    return lhs.lcp < rhs.lcp;
                }
                else
                {
                    return lhs.i < rhs.i;
                }
            });

            sort(item2.begin(), item2.end(), [](const LCPINTV &lhs, const LCPINTV &rhs) {
                if (lhs.lcp != rhs.lcp)
                {
                    return lhs.lcp < rhs.lcp;
                }
                else
                {
                    return lhs.i < rhs.i;
                }
            });

            bool b = true;
            if (item1.size() != item2.size())
            {
                std::cout << "Distinct Size!" << item1.size() << "/" << item2.size() << std::endl;
                b = false;
            }
            else
            {
                for (uint64_t i = 0; i < item1.size(); i++)
                {

                    if (item1[i].i != item2[i].i || item1[i].j != item2[i].j || item1[i].lcp != item2[i].lcp)
                    {
                        std::cout << "Distinct Value!" << std::endl;
                        std::cout << item1[i].to_string() << "/" << item2[i].to_string() << std::endl;

                        b = false;
                        break;
                    }
                }
            }
            if (!b)
            {
                if (item1.size() < 100)
                {
                    std::cout << "Test: " << item1.size() << "/ Collect: " << item2.size() << std::endl;

                    std::cout << "Test Vec:" << std::endl;
                    for (auto it : item1)
                    {
                        std::cout << it.to_string();
                    }
                    std::cout << std::endl;

                    std::cout << "Correct Vec:" << std::endl;
                    for (auto it : item2)
                    {
                        std::cout << it.to_string();
                    }
                    std::cout << std::endl;
                }
                throw -1;
            }
            return true;
        }
    } // namespace beller

} // namespace stool