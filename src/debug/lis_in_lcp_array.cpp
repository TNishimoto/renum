// License: MIT http://opensource.org/licenses/MIT
/*
  This code was copied from https://takeda25.hatenablog.jp/entry/20101202/1291269994 and I modified it.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "stool/src/cmdline.h"
#include "../test/old_postorder_maximal_substrings.hpp"
#include "libdivsufsort/sa.hpp"
#include "../postorder_maximal_substring_intervals.hpp"
#include "../../module/rlbwt_iterator/src/include/rlbwt_iterator.hpp"
#include "../forward_bwt.hpp"
using namespace std;
using CHAR = char;
using INDEX = uint64_t;

vector<uint64_t> longest_increasing_subsequence(const vector<uint64_t> &a, bool strict)
{
    vector<uint64_t> lis;
    for (auto &p : a)
    {
        vector<uint64_t>::iterator it;
        if (strict)
            it = lower_bound(begin(lis), end(lis), p);
        else
            it = upper_bound(begin(lis), end(lis), p);
        if (end(lis) == it)
            lis.emplace_back(p);
        else
            *it = p;
    }
    /*
  for(auto value : lis){
    std::cout << value << ", ";
  }
  std::cout << std::endl;
  */
    return lis;
}

int main(int argc, char *argv[])
{

    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file name", true);
    p.add<string>("output_file", 'o', "output file name", false, "");

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    string outputFile = p.get<string>("output_file");


    if (outputFile.size() == 0)
    {
        outputFile = inputFile + ".lis.txt";
    }


    std::ofstream out(outputFile, std::ios::out | std::ios::binary);
    if (!out)
    {
        throw std::runtime_error("Cannot open the output file!");
    }

    vector<CHAR> T = stool::load_char_vec_from_file(inputFile, true); // input text
    std::cout << "Constructing SA..." << std::endl;
    std::vector<INDEX> sa = stool::construct_suffix_array(T);
    std::cout << "Constructing LCP array..." << std::endl;
    std::vector<INDEX> lcpArray = stool::constructLCP<CHAR, INDEX>(T, sa);

    /** LIS **/
    std::cout << "Computing the longest increasing subsequence in the LCP array..." << std::endl;
    auto lis_vec = longest_increasing_subsequence(lcpArray, true);

    string s = "";
    for(uint64_t i=0;i<lis_vec.size();i++){
        out << std::to_string(lis_vec[i]);
        if(i < lis_vec.size() - 1){
            out << ", ";
        }

    }
    out.close();


    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "File \t\t\t\t\t : " << inputFile << std::endl;
    std::cout << "Output \t\t\t\t\t : " << outputFile << std::endl;
    std::cout << "The length of the input text \t\t : " << T.size() << std::endl;
    std::cout << "LIS: " << lis_vec.size() << std::endl;

    std::cout << "_______________________________________________________" << std::endl;
    std::cout << "\033[39m" << std::endl;
}