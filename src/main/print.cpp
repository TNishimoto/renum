// License: MIT http://opensource.org/licenses/MIT

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "stool/src/cmdline.h"
#include "divsufsort.h"
#include "divsufsort64.h"
#include "stool/src/io.hpp"
#include "stool/src/sa_bwt_lcp.hpp"
#include "libdivsufsort/sa.hpp"
#include "common.hpp"

//#include "../minimal_substrings/naive_minimal_substrings.hpp"

using namespace std;
using namespace stool;

using INDEXTYPE = int64_t;
using CHAR = char;



int main(int argc, char *argv[])
{

    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file name", true);
    p.add<string>("lcp_interval_file", 'l', "LCP interval file name", true);
    //p.add<string>("tree_file", 't', "file type", false, "NULL");

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    string intervalFile = p.get<string>("lcp_interval_file");
    //string type = p.get<string>("tree_file");

    std::vector<char> T = stool::load_text(inputFile); // input text
    std::vector<LCPInterval<uint64_t>> intervals;
    stool::load_vector<LCPInterval<uint64_t>>(intervalFile, intervals, false, true);

    std::vector<uint64_t> sa = stool::construct_suffix_array(T);
    stool::esaxx::print<char, uint64_t>(intervals, T, sa);
    stool::esaxx::printText<char>(T);
    stool::esaxx::printColor<char, uint64_t>(intervals, T, sa, true);

        /*
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
            intervals[i].print(i, T, sa);
        }
        */
        

    return 0;
}