
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "stool/src/cmdline.h"
#include "stool/src/io.hpp"
//#include "../postorder_suffix_tree.hpp"
#include "libdivsufsort/sa.hpp"
#include "../minimal_substrings/postorder_special_suffix_tree.hpp"
#include "../minimal_substrings/minimal_substring_iterator.hpp"
//#include "../minimal_substrings/naive_minimal_substrings.hpp"
//#include "../minimal_substrings/minimal_substring_tree.hpp"
#include "stool/src/sa_bwt_lcp.hpp"
#include "common.hpp"

using namespace std;
using INDEXTYPE = uint64_t;

int main(int argc, char *argv[])
{
  using CHAR = uint8_t;
  using INDEX = uint64_t;
  cmdline::parser p;
  p.add<string>("input_file", 'i', "input file name", true);
  p.add<string>("output_file", 'o', "output file name", false, "");
  p.add<bool>("print", 'p', "print info", false, true);
  p.add<string>("format", 'f', "output format (binary or csv)", false, "binary");

  p.parse_check(argc, argv);
  string inputFile = p.get<string>("input_file");
  string outputFile = p.get<string>("output_file");
  string format = p.get<string>("format");
  bool isPrint = p.get<bool>("print");
  if (outputFile.size() == 0)
  {
    if (format == "csv")
    {
      outputFile = inputFile + ".min.csv";
    }
    else
    {
      outputFile = inputFile + ".min";
    }
  }
  if (format != "csv")
    format = "binary";


  auto start = std::chrono::system_clock::now();
  vector<CHAR> T = stool::load_text_from_file(inputFile, true); // input text
  std::cout << "Constructing Suffix Array" << std::endl;
  std::vector<INDEX> sa = stool::construct_suffix_array(T);
  //std::cout << "Constructing LCP Array" << std::endl;
  std::vector<INDEX> lcpArray = stool::constructLCP<CHAR, INDEX>(T, sa);
  std::cout << "Constructing BWT" << std::endl;
  std::vector<CHAR> bwt = stool::constructBWT<CHAR, INDEX>(T, sa);
  //vector<stool::LCPInterval<INDEX>> minimalSubstrings = stool::esaxx::MinimalSubstringIterator<CHAR, INDEX, vector<INDEX>>::constructSortedMinimalSubstrings(bwt, sa, lcpArray);
  //stool::write_vector(outputFile, minimalSubstrings, false);

  
  /*
  auto end = std::chrono::system_clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  if (format == "csv")
  {
    //string otext = "";

    std::vector<char> plainT;
    for(uint64_t i=0;i<T.size();i++){
      plainT.push_back((char)T[i]);
    }

    if(isPrint){
      std::cout << "Minimal substrings in the file" << std::endl;
      stool::esaxx::print<char, INDEX>(minimalSubstrings,plainT,sa);
    }
    stool::esaxx::writeText<char, INDEX>(outputFile, minimalSubstrings, plainT, sa);
  }
  else
  {
    stool::write_vector(outputFile, minimalSubstrings, false);
  }
  std::cout << "\033[36m";
  std::cout << "___________RESULT___________" << std::endl;
  std::cout << "File: " << inputFile << std::endl;
  std::cout << "Output: " << outputFile << std::endl;
  std::cout << "Output format: " << format << std::endl;
  std::cout << "The length of the input text: " << T.size() << std::endl;
  std::cout << "The number of minimum substrings: " << minimalSubstrings.size() << std::endl;
  std::cout << "Excecution time : " << elapsed << "ms";
  double charperms = (double)T.size() / elapsed;
  std::cout << "[" << charperms << "chars/ms]" << std::endl;
  std::cout << "_________________________________" << std::endl;
  std::cout << "\033[39m" << std::endl;
  */
  
}