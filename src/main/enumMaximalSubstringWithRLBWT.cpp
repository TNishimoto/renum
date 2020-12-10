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
#include "common.hpp"
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

/*
void iterateMSwithSDSL(){

}
*/



uint64_t input_text_size = 0;
/*
double sa_construction_time = 0;
double lcp_array_construction_time = 0;
double bwt_construction_time = 0;
double ms_construction_time = 0;
*/
std::vector<std::pair<std::string, uint64_t>> execution_time_messages;

uint64_t iterateMSWithRLBWT(string filename, std::ofstream &out){

  using LCP = stool::rlbwt::ForwardLCPArray<std::vector<INDEX>>;
  //using SA = stool::rlbwt::ForwardSA<std::vector<INDEX>>;
  auto start_prep = std::chrono::system_clock::now();
  stool::rlbwt::RLBWT<std::vector<CHAR>, std::vector<INDEX> > rlestr = stool::rlbwt::Constructor::load_RLBWT_from_file<CHAR, INDEX>(filename);
  auto end_prep = std::chrono::system_clock::now();
  double prep_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_prep - start_prep).count();
  execution_time_messages.push_back(std::pair<std::string, uint64_t>("RLBWT loading time\t\t", prep_time));

  input_text_size = rlestr.str_size();

  auto start_lcp = std::chrono::system_clock::now();
  LCP lcpArray;
  auto stime = lcpArray.construct_from_rlbwt(&rlestr, true);
  auto end_lcp = std::chrono::system_clock::now();
  double lcp_array_construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_lcp - start_lcp).count();
  execution_time_messages.push_back(std::pair<std::string, uint64_t>("Sampling SA construction time\t", stime.second));
  execution_time_messages.push_back(std::pair<std::string, uint64_t>("Sampling LCP construction time\t", stime.first));

  execution_time_messages.push_back(std::pair<std::string, uint64_t>("Sampling LCP & SA construction time\t", lcp_array_construction_time));
  
  
  //const SA* sa_pointer = lcpArray.get_ForwardSA();

  using BWT_RLBWT = stool::rlbwt::ForwardBWT<std::vector<CHAR>, std::vector<INDEX>>;
  BWT_RLBWT bwt_rlbwt(&rlestr);
  
  auto start_ms = std::chrono::system_clock::now();
  stool::esaxx::PostorderMaximalSubstringIntervals<CHAR, INDEX, LCP, BWT_RLBWT > pmsi;
  pmsi.construct(&lcpArray, &bwt_rlbwt);
  uint64_t count = 0;
  for (auto it : pmsi)
  {
	  out.write(reinterpret_cast<const char *>(&it), sizeof(stool::LCPInterval<INDEX>));
    ++count;
  }
  auto end_ms = std::chrono::system_clock::now();
  double ms_construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_ms - start_ms).count();
  execution_time_messages.push_back(std::pair<std::string, uint64_t>("MS Construction time\t\t\t", ms_construction_time));
  return count;  
}

int main(int argc, char *argv[])
{

  cmdline::parser p;
  p.add<string>("input_file", 'i', "input RLBWT file name", true);
  p.add<string>("output_file", 'o', "output file name", false, "");
  //p.add<string>("format", 'f', "output format (binary or csv)", false, "binary");

  p.parse_check(argc, argv);
  string inputFile = p.get<string>("input_file");
  string outputFile = p.get<string>("output_file");
  //string format = p.get<string>("format");
  string format = "binary";

  if (format != "binary")
  {
    format = "csv";
  }


  if (outputFile.size() == 0)
  {
    if (format == "csv")
    {
      outputFile = inputFile + ".max.csv";
    }
    else
    {
      outputFile = inputFile + ".max";
    }
  }

  auto start = std::chrono::system_clock::now();


	std::ofstream out(outputFile, std::ios::out | std::ios::binary);
	if (!out){
    throw std::runtime_error("Cannot open the output file!");
  }
  uint64_t ms_count = 0;
  std::vector<stool::LCPInterval<INDEX>> intervals;

  ms_count = iterateMSWithRLBWT(inputFile, out);

  auto end = std::chrono::system_clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


  std::cout << "\033[31m";
  std::cout << "______________________RESULT______________________" << std::endl;
  //std::cout << "mode \t\t\t\t\t : " << mode << std::endl; 
  std::cout << "RLBWT File \t\t\t\t\t : " << inputFile << std::endl;
  std::cout << "Output \t\t\t\t\t : " << outputFile << std::endl;
  //std::cout << "Output format \t\t\t\t : " << format << std::endl;  
  std::cout << "The length of the input text \t\t : " << input_text_size << std::endl;
  //std::cout << "The number of maximum substrings: " << maximumSubstringCount << std::endl;
  std::cout << "The number of maximum substrings \t : " << ms_count << std::endl;
  std::cout << "Excecution time \t\t\t : " << elapsed << "[ms]" << std::endl;
  for(auto it : execution_time_messages){
  std::cout << "|\t " << it.first << " : " << it.second << "[ms]" << std::endl;

  }
  
  std::cout << "_______________________________________________________" << std::endl;
  std::cout << "\033[39m" << std::endl;
}