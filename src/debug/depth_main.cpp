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
#include "../main/common.hpp"
#include "libdivsufsort/sa.hpp"
#include "../postorder_maximal_substring_intervals.hpp"
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp_dac.hpp>
#include <sdsl/lcp_support_sada.hpp>
#include "../../module/rlbwt_iterator/src/include/rlbwt_iterator.hpp"
#include "../forward_bwt.hpp"
using namespace std;
using CHAR = char;
using INDEX = uint64_t;


uint64_t input_text_size = 0;
std::vector<std::pair<std::string, uint64_t>> execution_time_messages;

uint64_t iterateMS(string filename){
  vector<CHAR> T = stool::load_char_vec_from_file(filename, true); // input text
  input_text_size = T.size();

  std::vector<INDEX> sa = stool::construct_suffix_array(T);
  using BWT = stool::esaxx::ForwardBWT<CHAR, std::vector<CHAR>, std::vector<INDEX>>;
  BWT bwt(&T,&sa);

  std::vector<INDEX> lcpArray = stool::constructLCP<CHAR, INDEX>(T, sa);
  stool::esaxx::PostorderSuffixTreeIntervals<INDEX, std::vector<INDEX>, std::vector<INDEX>> st;
  st.construct(&sa, &lcpArray);

    auto beg = st.begin();
    auto end = st.end();

    uint64_t maxStack = 0;
    while(beg != end ){
        if(maxStack < beg.get_incomplete_stack_size()){
            maxStack = beg.get_incomplete_stack_size();
        }
        ++beg;
    }
    return maxStack;

}

int main(int argc, char *argv[])
{

  cmdline::parser p;
  p.add<string>("input_file", 'i', "input file name", true);

  p.parse_check(argc, argv);
  string inputFile = p.get<string>("input_file");

  auto start = std::chrono::system_clock::now();
  uint64_t ms_count = iterateMS(inputFile);



  std::cout << "\033[31m";
  std::cout << "______________________RESULT______________________" << std::endl;
  std::cout << "File \t\t\t\t\t : " << inputFile << std::endl;
  std::cout << "The length of the input text \t\t : " << input_text_size << std::endl;
  std::cout << "Nest \t : " << ms_count << std::endl; 
  std::cout << "_______________________________________________________" << std::endl;
  std::cout << "\033[39m" << std::endl;
}