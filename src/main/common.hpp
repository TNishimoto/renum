#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include "stool/src/io.hpp"
#include "stool/src/sa_bwt_lcp.hpp"

namespace stool
{
namespace esaxx
{


template <typename CHAR, typename INDEX, typename SA = std::vector<INDEX>>
std::vector<CHAR> constructBWT(const std::vector<CHAR> &text, const SA &sa)
{
    std::vector<CHAR> bwt;
    bwt.resize(text.size());
    INDEX n = text.size();
    INDEX i = 0;
    for(auto saValue : sa){
        if (saValue == 0)
        {
            bwt[i] = text[n - 1];
        }
        else
        {
            bwt[i] = text[saValue - 1];
        }

        ++i;
    }
    /*
    for (INDEX i = 0; i < text.size(); i++)
    {
        if (sa[i] == 0)
        {
            bwt[i] = text[n - 1];
        }
        else
        {
            bwt[i] = text[sa[i] - 1];
        }
    }
    */
    return bwt;
}
template <typename CHAR, typename INDEX>
void print(std::vector<stool::LCPInterval<INDEX>> &intervals, std::vector<CHAR> &text, std::vector<INDEX> &sa)
{
    std::cout << "id"
              << "\t"
              << "occurrence"
              << "\t"
              << "range(SA)"
              << "\t"
              << "string length"
              << "\t"
              << "string" << std::endl;
    
    for (uint64_t i = 0; i < intervals.size(); i++)
    {
        std::string line = intervals[i].getCSVLine(i, text, sa);
        if (line.size() < 255)
        {
            std::cout << intervals[i].getCSVLine(i, text, sa) << std::endl;
        }
        else
        {
            std::cout << "Omission since the line length larger than 255." << std::endl;
        }
        if (i > 1000)
        {
            std::cout << "Omit remaining lines." << std::endl;
            break;
        }
    }
}
/*
template <typename CHAR, typename INDEX>
std::string toLogLine(uint64_t id, std::vector<CHAR> &text, std::vector<INDEX> &sa, stool::LCPInterval<INDEX> &interval)
{
    //std::cout << interval.getCSVLine(id, text, sa) << std::endl;
    std::string T(text.begin(), text.end());
    std::string log = "";

    std::string mstr = T.substr(sa[interval.i], interval.lcp);
    log.append("\"");
    log.append(mstr);
    log.append("\" ");

    log.append("SA[");
    log.append(std::to_string(interval.i));
    log.append(", ");
    log.append(std::to_string(interval.j));
    log.append("] ");
    //log.append(std::to_string(interval.lcp));

    if (interval.lcp == 0)
        return log;
    log.append("occ: ");

    std::vector<uint64_t> occs;
    for (uint64_t x = interval.i; x <= interval.j; x++)
    {
        uint64_t pos = sa[x];
        occs.push_back(pos);
    }
    std::sort(occs.begin(), occs.end());

    for (uint64_t i = 0; i < occs.size(); i++)
    {
        uint64_t pos = occs[i];
        uint64_t endPos = occs[i] + interval.lcp - 1;
        std::string occ = "[" + std::to_string(pos) + ".." + std::to_string(endPos) + "]";
        log.append(occ);
    }

    return log;
}
*/

template <typename CHAR, typename INDEX>
void writeText(std::string filename, std::vector<stool::LCPInterval<INDEX>> &intervals, std::vector<CHAR> &text, std::vector<INDEX> &sa)
{
    std::string otext = "";

    otext += "id";
    otext += "\t";
    otext += "occurrence";
    otext += "\t";
    otext += "range(SA)";
    otext += "\t";
    otext += "string length";
    otext += "\t";
    otext += "string";
    otext += "\r\n";
    std::cout << "csv"<< intervals.size() << std::endl;
    for (size_t i = 0; i < intervals.size(); i++)
    {
        std::cout << intervals[i].getCSVLine(i, text, sa) << std::endl;
        otext.append(intervals[i].getCSVLine(i, text, sa));

        //otext.append(stool::esaxx::toLogLine<char,INDEX>(i, text, sa, intervals[i]));
        if (i + 1 != intervals.size())
            otext.append("\r\n");
    }
    std::ofstream out(filename, std::ios::out | std::ios::binary);
    out.write((const char *)(&otext[0]), sizeof(char) * otext.size());
}

template <typename CHAR>
void printText(std::vector<CHAR> &text)
{
    std::cout << "\033[42m";
    std::cout << "\033[30m";
    for (uint64_t i = 0; i < text.size(); i++)
    {
        std::cout << (text[i] == 0 ? '$' : (char)text[i]);
    }
    std::cout << "\033[49m";
    std::cout << "\033[39m";    
    std::cout << " Text";
    std::cout << std::endl;
}

template <typename CHAR, typename INDEX>
void printColor(std::vector<stool::LCPInterval<INDEX>> &intervals, std::vector<CHAR> &text, std::vector<INDEX> &sa, bool printFirstOccurrenceFlag = false, int64_t printVerticalLinePosition = -1, std::string intervalName = "")
{
    //uint64_t wholeFstOcc = text.size();
    if(intervalName != ""){
        for(uint64_t i=0;i<=text.size();i++) std::cout << " ";
        std::cout << intervalName << std::endl;
    }
    for (uint64_t i = 0; i < intervals.size(); i++)
    {
        stool::LCPInterval<INDEX> &interval = intervals[i];
        if (interval.lcp == 0)
            continue;

        std::string s;
        std::string ministr;
        s.resize(text.size(), ' ');
        ministr.resize(interval.lcp, ' ');
        uint64_t fstOcc = text.size();
        for (uint64_t x = interval.i; x <= interval.j; x++)
        {
            uint64_t pos = sa[x];
            if (pos < fstOcc)
                fstOcc = pos;
            for (uint64_t l = 0; l < interval.lcp; l++)
            {
                char c = text[pos + l] == 0 ? '$' : text[pos + l];
                s[pos + l] = c;
                ministr[l] = c;
            }
        }
        if (printFirstOccurrenceFlag)
        {
            for (uint64_t x = 0; x < fstOcc; x++)
            {
                if (s[x] == ' ')
                    s[x] = '-';
            }
        }
        if(printVerticalLinePosition != -1 && printVerticalLinePosition < (int64_t)s.size() && s[printVerticalLinePosition] == ' '){
                s[printVerticalLinePosition] = '|';
        }
        s += '(' + ministr + ')';

        std::cout << "\033[36m";
        std::cout << s;
        std::cout << "\033[39m" << std::endl;
        //if(fstOcc < wholeFstOcc) wholeFstOcc = fstOcc;
    }
    //return wholeFstOcc;
}

template <typename CHAR = uint8_t>
bool compare_suffixes_with_uint64(const std::vector<CHAR> &text, const uint64_t x, const uint64_t y)
{
    uint64_t max = x < y ? text.size() - y : text.size() - x;
    for (uint64_t i = 0; i < max; i++)
    {
        uint64_t c1 = text[x + i];
        uint64_t c2 = text[y + i];
        if (c1 != c2)
        {
            return c1 < c2;
        }
    }
    return x > y;
}
template <typename CHAR = uint8_t, typename INDEX = uint64_t>
std::vector<INDEX> construct_naive_SA_with_uint64(const std::vector<CHAR> &text)
{
    std::vector<INDEX> r;
    for (uint64_t i = 0; i < text.size(); i++)
    {
        r.push_back(i);
    }

    std::sort(
        r.begin(),
        r.end(),
        [&](const uint64_t &x, const uint64_t &y) {
            return compare_suffixes_with_uint64(text, x, y);
        });
    return r;
}
bool check_test(std::vector<char> &text)
{
  for (auto &it : text)
  {
    if (it < 0)
    {
      //throw std::logic_error("This text contains minus character!");
      return true;
    }
  }
  return false;
}

} // namespace esaxx
} // namespace stool