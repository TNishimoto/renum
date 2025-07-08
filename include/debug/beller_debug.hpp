#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <set>
#include <iostream>
#include "stool/include/stool.hpp"

#include "../basic/char_interval.hpp"

//#include "sa_lcp.hpp"
//using namespace std;
//using namespace sdsl;

namespace stool
{
    namespace renum
    {

        template <typename INDEX_SIZE>
        bool check(std::vector<CharInterval<INDEX_SIZE, uint8_t>> &vec1, std::vector<CharInterval<INDEX_SIZE, uint8_t>> &vec2)
        {
            std::sort(vec1.begin(), vec1.end(), [](const CharInterval<INDEX_SIZE, uint8_t> &lhs, const CharInterval<INDEX_SIZE, uint8_t> &rhs) {
                return lhs.c < rhs.c;
            });
            std::sort(vec2.begin(), vec2.end(), [](const CharInterval<INDEX_SIZE, uint8_t> &lhs, const CharInterval<INDEX_SIZE, uint8_t> &rhs) {
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
        bool equal_check_lcp_intervals(std::vector<LCPInterval<INDEX>> &item1, std::vector<LCPInterval<INDEX>> &item2, std::string name = "")
        {
            using LCPINTV = stool::LCPInterval<INDEX>;

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
                std::cout << name << " Error!" << std::endl;
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
    } // namespace renum

} // namespace stool