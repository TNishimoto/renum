#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <iostream>
#include <fstream>

namespace stool
{
    namespace renum
    {
        template <typename INDEX>
        struct RInterval
        {
            INDEX beginIndex;
            INDEX beginDiff;
            INDEX endIndex;
            INDEX endDiff;

            void print() const
            {
                std::cout << "[(" << this->beginIndex << ", " << this->beginDiff << "), (" << this->endIndex << ", " << this->endDiff << ")]" << std::endl;
            }

            bool is_special() const
            {
                return this->beginIndex == std::numeric_limits<INDEX>::max();
            }

            static RInterval get_special()
            {
                RInterval r;
                r.beginIndex = std::numeric_limits<INDEX>::max();
                return r;
            }
        };

    } // namespace renum
} // namespace stool