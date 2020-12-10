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


int main(int argc, char *argv[])
{

  cmdline::parser p;
  p.add<string>("input_file", 'i', "input file name", true);
  p.add<bool>("mode", 'm', "0 or 1(standard or faster algorithm)", false, true);
  
  p.parse_check(argc, argv);
  string inputFile = p.get<string>("input_file");
  bool mode = p.get<bool>("mode");


  auto start = std::chrono::system_clock::now();

  using LCP = stool::rlbwt::ForwardLCPArray<std::vector<INDEX>>;
  stool::rlbwt::RLBWT<std::vector<CHAR>, std::vector<INDEX> > rlestr = stool::rlbwt::Constructor::load_RLBWT_from_file<CHAR, INDEX>(inputFile);
  LCP lcpArray;
  auto stime = lcpArray.construct_from_rlbwt(&rlestr, mode);
  
  auto end = std::chrono::system_clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


  std::cout << "\033[31m";
  std::cout << "______________________RESULT______________________" << std::endl;
  std::cout << "mode \t\t\t\t\t : " <<( mode ? "Faster algorithm" : "Standard algorithm") << std::endl; 

  std::cout << "RLBWT File \t\t\t\t\t : " << inputFile << std::endl;
  std::cout << "The length of the input text \t\t : " << rlestr.str_size() << std::endl;
  std::cout << "Excecution time \t\t\t : " << elapsed << "[ms]" << std::endl;

  
  
  std::cout << "_______________________________________________________" << std::endl;
  std::cout << "\033[39m" << std::endl;
}