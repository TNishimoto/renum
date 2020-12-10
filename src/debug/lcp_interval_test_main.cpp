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
#include "../../module/rlbwt_iterator/src/include/rlbwt_iterator.hpp"

using namespace std;
using CHAR = uint8_t;
using INDEX = uint64_t;
#include "../test/old_postorder_maximal_substrings.hpp"
#include "../main/common.hpp"
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp_dac.hpp>
#include <sdsl/lcp_support_sada.hpp>
#include "../postorder_suffix_tree_intervals.hpp"
#include "../test/naive_algorithms.hpp"
#include "../postorder_maximal_substring_intervals.hpp"
#include "../src/minimal_substrings/minimal_substring_iterator.hpp"
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp_dac.hpp>
#include <sdsl/lcp_support_sada.hpp>
#include "../forward_bwt.hpp"
#include "../postorder_lcp_intervals.hpp"

template <typename T, typename U>
bool super_equal_check(const T &vec1, const U &vec2)
{
  if (vec1.size() != vec2.size())
  {
    std::string s = std::string("String sizes are different!") + ", collect = " + std::to_string(vec1.size()) + ", test = " + std::to_string(vec2.size());

    throw std::logic_error(s);
  }
  for (uint64_t i = 0; i < vec1.size(); i++)
  {
    if (vec1[i] != vec2[i])
    {
      std::string msg = "collect_vec[" + std::to_string(i) + "] != test_vec[" + std::to_string(i) + "]";

      throw std::logic_error("Values are different! " + msg);
    }
  }
  return true;
}


int main(int argc, char *argv[])
{

  cmdline::parser p;
  p.add<string>("input_file", 'i', "input file name", true);

  p.parse_check(argc, argv);
  string filename = p.get<string>("input_file");

  vector<char> text = stool::load_char_vec_from_file(filename, true); // input text
  bool is_contained_minus_character = stool::esaxx::check_test(text);
  if (is_contained_minus_character)
  {
    std::cout << "This text contains minus character!" << std::endl;
  }

  using SA_PLAIN = std::vector<INDEX>;
  using LCP_PLAIN = std::vector<INDEX>;
  using LCPINTERVALS = vector<stool::LCPInterval<INDEX>>;
  using LCPINTERVAL = stool::LCPInterval<INDEX>;

  SA_PLAIN sa_naive = stool::esaxx::construct_naive_SA_with_uint64<char, INDEX>(text);
  LCP_PLAIN lcp_naive = stool::constructLCP<char, INDEX>(text, sa_naive);


  //std::cout << "LCP interval test: " << test_name << std::flush;
  vector<stool::LCPInterval<INDEX>> correct_intervals = stool::esaxx::naive_compute_complete_lcp_intervals<char, INDEX>(text, sa_naive);

  vector<stool::LCPInterval<INDEX>> test_intervals = stool::esaxx::PostorderLCPIntervals<INDEX, LCP_PLAIN>::compute_lcp_intervals(lcp_naive);

  
  std::vector<LCPINTERVAL> test_lcp_intervals;
  
  for(auto it : test_intervals){
      test_lcp_intervals.push_back(it);
  }
  
  stool::sort_in_preorder(test_lcp_intervals);
  stool::equal_check(correct_intervals, test_lcp_intervals);
  std::cout << "[OK!]" << std::endl;


  using LCP_SDSL = sdsl::lcp_dac<>;
  LCP_SDSL lcp_sdsl;
  construct(lcp_sdsl, filename, 1);
  vector<stool::LCPInterval<INDEX>> test_intervals_sdsl = stool::esaxx::PostorderLCPIntervals<INDEX, LCP_SDSL>::compute_lcp_intervals(lcp_sdsl);
  stool::sort_in_preorder(test_intervals_sdsl);
  stool::equal_check(correct_intervals, test_intervals_sdsl);
  std::cout << "[OK!]" << std::endl;


  using LCP_RLBWT = stool::rlbwt::ForwardLCPArray<std::vector<INDEX>>;
  stool::rlbwt::RLBWT<std::vector<char>, std::vector<INDEX>> rlestr;
  stool::rlbwt::Constructor::construct_from_file<>(rlestr, filename);
  LCP_RLBWT lcp_rlbwt;
  lcp_rlbwt.construct_from_rlbwt(&rlestr, false);
  vector<stool::LCPInterval<INDEX>> test_intervals_rlbwt = stool::esaxx::PostorderLCPIntervals<INDEX, LCP_RLBWT>::compute_lcp_intervals(lcp_rlbwt);
  stool::sort_in_preorder(test_intervals_rlbwt);
  stool::equal_check(correct_intervals, test_intervals_rlbwt);
  std::cout << "[OK!]" << std::endl;


}