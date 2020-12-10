#include <iostream>
#include <string>
#include <memory>
#include "stool/src/print.hpp"
#include "stool/src/cmdline.h"
#include "stool/src/io.hpp"
#include "stool/src/debug.hpp"
#include "../test/naive_algorithms.hpp"
#include "../src/minimal_substrings/minimal_substring_iterator.hpp"

using namespace std;
using INDEX = uint64_t;
template <typename CHAR>
void lcp_interval_test(vector<CHAR> &text)
{
    //stool::Printer::print("text", text);

    std::vector<INDEX> sa = stool::construct_naive_SA<CHAR, INDEX>(text);
    vector<stool::LCPInterval<INDEX>> correct_intervals = stool::esaxx::naive_compute_lcp_intervals<CHAR, INDEX>(text, sa);


    std::vector<INDEX> lcpArray = stool::constructLCP<CHAR, INDEX>(text, sa);
  vector<stool::LCPInterval<INDEX>> test_intervals = stool::esaxx::PostorderSuffixTreeIntervals<INDEX, std::vector<INDEX>, std::vector<INDEX> >::compute_lcp_intervals(sa, lcpArray);
  stool::sort_in_preorder(test_intervals);
    //vector<stool::LCPInterval<INDEX>> test_intervals = stool::esaxx::compute_preorder_lcp_intervals<T, INDEX>(text, sa);

    /*
    for(auto& p : correct_intervals){
        std::cout << p.to_string() << std::endl;
    }
    std::cout << std::endl;


    for(auto& p : test_intervals){
        std::cout << p.to_string() << std::endl;
    }
    std::cout << std::endl;
    */
    stool::equal_check(correct_intervals, test_intervals);
}

template <typename CHAR>
void minimal_substrings_test(vector<CHAR> &text)
{    
    std::vector<INDEX> sa = stool::construct_naive_SA<CHAR, INDEX>(text);
    std::vector<INDEX> lcpArray = stool::constructLCP<CHAR, INDEX>(text, sa);
    vector<stool::LCPInterval<INDEX>> correct_intervals = stool::esaxx::naive_compute_minimal_substrings<CHAR, INDEX>(text, sa);
    vector<stool::LCPInterval<INDEX>> test_intervals = stool::esaxx::compute_minimal_substrings<CHAR, INDEX>(text, sa, lcpArray);
    stool::sort_in_preorder(test_intervals);

    try
    {
        stool::equal_check(correct_intervals, test_intervals);
    }
    catch (const std::logic_error &e)
    {
        std::cout << std::endl;
        std::cout << "Minimal substrings Error!" << std::endl;
        stool::Printer::print("text", text);

        std::cout << "correct_intervals" << std::endl;
        for (auto &p : correct_intervals)
        {
            std::cout << p.to_string() << std::endl;
        }
        std::cout << std::endl;

        std::cout << "test_intervals" << std::endl;
        for (auto &p : test_intervals)
        {
            std::cout << p.to_string() << std::endl;
        }
        std::cout << std::endl;

        throw e;
    }
}
void minimal_substrings_test(uint64_t loop, uint64_t size){
    std::cout << "Minimal substrings tests" << std::endl;

    std::cout << ":Short string" << std::endl;
    for (size_t i = 0; i < 100; i++)
    {
        for (uint64_t alphabet = 1; alphabet < 64; alphabet *= 2)
        {
            if (i % 100 == 0)
                std::cout << "+" << std::flush;
            std::vector<char> text = stool::create_deterministic_integers<char>(i, alphabet, 0, i);
            text.push_back(std::numeric_limits<char>::min());
            minimal_substrings_test(text);
        }
    }
    std::cout << std::endl;

    
    std::cout << ":Char string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<char> text = stool::create_deterministic_integers<char>(size, 64, 0, i);
        text.push_back(std::numeric_limits<char>::min());
        minimal_substrings_test(text);
    }
    std::cout << std::endl;


    std::cout << ":uint8_t string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<uint8_t> text = stool::create_deterministic_integers<uint8_t>(size, 255, 1, i);
        text.push_back(std::numeric_limits<uint8_t>::min());

        minimal_substrings_test(text);
    }
    std::cout << std::endl;

    std::cout << ":int32_t string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<int32_t> text = stool::create_deterministic_integers<int32_t>(size, 255, -255, i);
        text.push_back(std::numeric_limits<int32_t>::min());
        minimal_substrings_test(text);
    }
    std::cout << std::endl;

    std::cout << ":uint32_t string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<uint32_t> text = stool::create_deterministic_integers<uint32_t>(size, 510, 1, i);
        text.push_back(std::numeric_limits<uint32_t>::min());
        minimal_substrings_test(text);
    }
    std::cout << std::endl;

    std::cout << ":int64_t string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<int64_t> text = stool::create_deterministic_integers<int64_t>(size, 1024, -1024, i);
        text.push_back(std::numeric_limits<int64_t>::min());
        minimal_substrings_test(text);
    }
    std::cout << std::endl;

    std::cout << ":uint64_t string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<uint64_t> text = stool::create_deterministic_integers<uint64_t>(size, 2048, 1, i);
        text.push_back(std::numeric_limits<uint64_t>::min());
        minimal_substrings_test(text);
    }
    std::cout << std::endl;
}

int main(int argc, char *argv[])
{

    cmdline::parser p;
    p.add<uint64_t>("size", 'l', "text length", true);
    p.parse_check(argc, argv);
    uint64_t size = p.get<uint64_t>("size");
    uint64_t loop = 10000;


    std::cout << "Suffix tree intervals tests" << std::endl;
    std::cout << ":Short string" << std::endl;
    for (size_t i = 0; i < 100; i++)
    {
        for (uint64_t alphabet = 1; alphabet < 64; alphabet *= 2)
        {
            if (i % 100 == 0)
                std::cout << "+" << std::flush;
            std::vector<char> text = stool::create_deterministic_integers<char>(i, alphabet, 0, i);
            lcp_interval_test(text);
        }
    }
    std::cout << std::endl;

    std::cout << ":char string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<char> text = stool::create_deterministic_integers<char>(size, 64, 0, i);
        lcp_interval_test(text);
    }

    std::cout << std::endl;


    std::cout << ":uint8_t string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<uint8_t> text = stool::create_deterministic_integers<uint8_t>(size, 255, 0, i);
        lcp_interval_test(text);
    }
    std::cout << std::endl;

    std::cout << ":int32_t string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<int32_t> text = stool::create_deterministic_integers<int32_t>(size, 255, -255, i);
        lcp_interval_test(text);
    }
    std::cout << std::endl;

    std::cout << ":uint32_t string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<uint32_t> text = stool::create_deterministic_integers<uint32_t>(size, 510, 0, i);
        lcp_interval_test(text);
    }
    std::cout << std::endl;

    std::cout << ":int64_t string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<int64_t> text = stool::create_deterministic_integers<int64_t>(size, 1024, -1024, i);
        lcp_interval_test(text);
    }
    std::cout << std::endl;

    std::cout << ":uint64_t string" << std::endl;
    for (size_t i = 0; i < loop; i++)
    {
        if (i % 100 == 0)
            std::cout << "+" << std::flush;
        std::vector<uint64_t> text = stool::create_deterministic_integers<uint64_t>(size, 2048, 0, i);
        lcp_interval_test(text);
    }
    std::cout << std::endl;

    minimal_substrings_test(loop, size);
    std::cout << "END" << std::endl;
}