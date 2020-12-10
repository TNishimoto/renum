
#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include "stool/src/io.hpp"
#include "stool/src/debug.hpp"
#include "stool/src/sa_bwt_lcp.hpp"

#include <stack>

namespace stool
{
  namespace esaxx
  {

    template <typename CHAR, typename INDEX = uint64_t>
    std::vector<stool::LCPInterval<INDEX>> naive_compute_lcp_intervals(const std::vector<CHAR> &text, const std::vector<INDEX> &sa)
    {
      

      std::vector<stool::LCPInterval<INDEX>>
          r;
      std::vector<INDEX> lcpArray = stool::constructLCP<CHAR, INDEX>(text, sa);
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
            r.push_back(stool::LCPInterval<INDEX>(i, x - 1, current_lcp));
            current_lcp = lcp;
          }

          if (current_lcp <= limit_lcp)
          {
            break;
          }
        }
      }
      r.push_back(stool::LCPInterval<INDEX>(0, sa.size() - 1, 0));

      std::sort(
          r.begin(),
          r.end(),
          LCPIntervalPreorderComp<INDEX>());

      return r;
    }

    template <typename CHAR, typename INDEX = uint64_t>
    std::vector<stool::LCPInterval<INDEX>> naive_compute_complete_lcp_intervals(const std::vector<CHAR> &text, const std::vector<INDEX> &sa)
    {
      std::vector<stool::LCPInterval<INDEX>> r = naive_compute_lcp_intervals(text, sa);

      std::vector<stool::LCPInterval<INDEX>> correct_lcp_intervals;
      for (auto it : r)
      {
        if (it.j - it.i != 0)
        {
          correct_lcp_intervals.push_back(it);
        }
        else
        {
          it.lcp = std::numeric_limits<INDEX>::max() - 1;
          correct_lcp_intervals.push_back(it);
        }
      }

      return correct_lcp_intervals;
    }

    template <typename CHAR, typename INDEX = uint64_t>
    std::vector<stool::LCPInterval<INDEX>> naive_compute_minimal_substrings(const std::vector<CHAR> &text, const std::vector<INDEX> &sa)
    {
      std::vector<stool::LCPInterval<INDEX>> r;

      for (uint64_t b = 0; b < text.size(); b++)
      {
        std::vector<CHAR> pattern;
        pattern.push_back(text[b]);
        stool::LCPInterval<INDEX> charInterval = stool::LCPInterval<INDEX>::computeLCPInterval(text, pattern, sa);
        r.push_back(charInterval);
        uint64_t prevOcc = charInterval.j - charInterval.i + 1;
        for (uint64_t e = b + 1; e < text.size(); e++)
        {
          if (prevOcc == 1)
            break;
          pattern.push_back(text[e]);
          stool::LCPInterval<INDEX> interval = stool::LCPInterval<INDEX>::computeLCPInterval(text, pattern, sa);
          uint64_t occ = interval.j - interval.i + 1;
          if (prevOcc > occ)
          {
            std::vector<CHAR> rightPattern;
            for (uint64_t x = 1; x < pattern.size(); x++)
              rightPattern.push_back(pattern[x]);
            stool::LCPInterval<INDEX> rightInterval = stool::LCPInterval<INDEX>::computeLCPInterval(text, rightPattern, sa);
            uint64_t rightOcc = rightInterval.j - rightInterval.i + 1;
            if (rightOcc > occ)
            {
              r.push_back(interval);
            }
          }
          prevOcc = occ;
        }
      }

      if (text.size() > 0)
        r.push_back(stool::LCPInterval<INDEX>(0, text.size() - 1, 0));


      std::sort(
          r.begin(),
          r.end(),
          LCPIntervalPreorderComp<INDEX>());

      std::vector<stool::LCPInterval<INDEX>> r2;
      for (auto &it : r)
      {
        if (r2.size() == 0)
        {
          r2.push_back(it);
        }
        else if (r2.size() > 0 && r2[r2.size() - 1] != it)
        {
          r2.push_back(it);
        }
      }
      return r2;
    }
    template <typename CHAR, typename INDEX = uint64_t>
    std::vector<stool::LCPInterval<INDEX>> naive_compute_minimal_substrings_with_uint64(const std::vector<CHAR> &text, const std::vector<INDEX> &sa)
    {
      std::vector<uint64_t> text2;
      for (auto &it : text)
        text2.push_back(it);
      return naive_compute_minimal_substrings(text2, sa);
    }

    template <typename CHAR, typename INDEX = uint64_t>
    std::vector<stool::LCPInterval<INDEX>> naive_compute_maximal_substrings(const std::vector<CHAR> &text, const std::vector<INDEX> &sa)
    {
      std::vector<stool::LCPInterval<INDEX>> r;

      std::vector<stool::LCPInterval<INDEX>> lcpIntervals = naive_compute_lcp_intervals(text, sa);

      for (auto interval : lcpIntervals)
      {
        uint64_t fst_index = sa[interval.i];
        CHAR fstBWT = fst_index == 0 ? text[text.size() - 1] : text[fst_index - 1];
        bool b1 = false;
        for (uint64_t x = interval.i + 1; x <= interval.j; x++)
        {
          uint64_t current_index = sa[x];
          CHAR currentBWT = current_index == 0 ? text[text.size() - 1] : text[current_index - 1];
          if (fstBWT != currentBWT)
          {
            b1 = true;
            break;
          }
        }
        if (interval.lcp == text.size())
          b1 = true;
        if (b1)
        {
          r.push_back(interval);
        }
      }

      std::sort(
          r.begin(),
          r.end(),
          LCPIntervalPreorderComp<INDEX>());

      return r;
    }

  } // namespace esaxx
} // namespace stool
