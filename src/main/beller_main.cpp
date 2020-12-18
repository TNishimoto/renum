#include <cassert>
#include <chrono>
#include "../module/stool/src/io.hpp"
#include "../module/stool/src/sa_bwt_lcp.hpp"

#include "../module/stool/src/print.hpp"
#include "../module/stool/src/cmdline.h"
#include "../module/stool/src/debug.hpp"
#include "../module/libdivsufsort/sa.hpp"

//#include "hpp/bwt.hpp"
#include "../beller/fmindex.hpp"
#include "../beller/beller_interval.hpp"
#include "../beller/beller_debug.hpp"

#include "../main/common.hpp"
#include "../test/naive_algorithms.hpp"
//#include "../postorder_maximal_substring_intervals.hpp"
//#include "../forward_bwt.hpp"

using namespace std;
//using namespace stool;
//using namespace stool::rlbwt;

using CHAR = char;
using INDEX = uint64_t;
using LCPINTV = stool::LCPInterval<uint64_t>;

std::vector<uint8_t> load_bwt(std::string filename)
{

  std::ifstream stream;
  stream.open(filename, std::ios::binary);

  std::vector<uint8_t> vec;

  if (!stream)
  {
    std::cerr << "error reading file " << std::endl;
    throw -1;
  }
  uint64_t len;
  stream.seekg(0, std::ios::end);
  uint64_t n = (unsigned long)stream.tellg();
  stream.seekg(0, std::ios::beg);
  len = n / sizeof(uint8_t);

  vec.resize(len, 0);
  stream.read((char *)&(vec)[0], len * sizeof(char));

  return vec;
}

void computeLCPIntervals(std::string inputFile, bool correctCheck)
{

  //string text = "";
  std::cout << "Loading : " << inputFile << std::endl;
  std::vector<uint8_t> text = stool::load_text_from_file(inputFile, true);
  vector<INDEX> sa = stool::construct_suffix_array(text);
  sdsl::int_vector<> bwt;
  stool::FMIndex::constructBWT(text, sa, bwt);

  auto test_Intervals = stool::beller::computeLCPIntervals<uint64_t>(bwt);
  test_Intervals.push_back(LCPINTV(0, text.size() - 1, 0));

  if (correctCheck)
  {
    auto correctLCP = stool::constructLCP(text, sa);
    std::cout << "Correct" << std::endl;
    std::vector<LCPINTV> correct_intervals = stool::beller::naive_compute_complete_lcp_intervals<uint64_t>(sa, correctLCP);
    stool::beller::equal_check_lcp_intervals(test_Intervals, correct_intervals);
    std::cout << "OK!" << std::endl;
  }
}

void computeMaximalSubstrings(std::string inputFile, std::string outputFile, bool correctCheck)
{

  //string text = "";
  auto start = std::chrono::system_clock::now();

  std::cout << "Loading : " << inputFile << std::endl;
  std::vector<uint8_t> text = stool::load_text_from_file(inputFile, true);
  vector<INDEX> sa = stool::construct_suffix_array(text);
  sdsl::int_vector<> bwt;

  stool::FMIndex::constructBWT(text, sa, bwt);

  std::vector<uint64_t> C;

  stool::FMIndex::constructC(bwt, C);

  wt_huff<> wt;
  construct_im(wt, bwt);
  uint64_t ms_count = 0;

  std::ofstream out(outputFile, std::ios::out | std::ios::binary);
  if (!out)
  {
    throw std::runtime_error("Cannot open the output file!");
  }
  uint64_t input_text_size = text.size();

  if (bwt.size() - 10 < UINT32_MAX)
  {
    ms_count = stool::beller::outputMaximalSubstrings<uint32_t>(bwt, out);
  }
  else
  {
    ms_count = stool::beller::outputMaximalSubstrings<uint64_t>(bwt, out);
  }
  auto end = std::chrono::system_clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

  /*
  auto test_Intervals = stool::beller::computeMaximalSubstrings(bwt, C, wt);
  test_Intervals.push_back(LCPINTV(0, text.size() - 1, 0));

  if (correctCheck)
  {
    using BWT = stool::esaxx::ForwardBWT<char, std::vector<char>, std::vector<INDEX>>;
    BWT bwt(&text, &sa);
    vector<stool::LCPInterval<INDEX>> correct_intervals = stool::esaxx::naive_compute_maximal_substrings<char, INDEX>(text, sa);
    stool::beller::equal_check_lcp_intervals(test_Intervals, correct_intervals);
    std::cout << "OK!" << std::endl;
  }
  */
  std::cout << "\033[31m";
  std::cout << "______________________RESULT______________________" << std::endl;
  std::cout << "RLBWT File \t\t\t\t\t : " << inputFile << std::endl;
  std::cout << "Output \t\t\t\t\t : " << outputFile << std::endl;
  std::cout << "The length of the input text \t\t : " << input_text_size << std::endl;
  std::cout << "The number of maximum substrings \t : " << ms_count << std::endl;
  std::cout << "Excecution time \t\t\t : " << elapsed << "[ms]" << std::endl;
  std::cout << "_______________________________________________________" << std::endl;
  std::cout << "\033[39m" << std::endl;
}

int main(int argc, char *argv[])
{
  cmdline::parser p;
  p.add<string>("input_file", 'i', "input file name", true);
  p.add<string>("mode", 'm', "mode", false, "xx");
  p.add<string>("output_file", 'o', "output file name", false, "");

  p.parse_check(argc, argv);
  string inputFile = p.get<string>("input_file");
  string mode = p.get<string>("mode");
  string outputFile = p.get<string>("output_file");
  string format = "binary";

  std::ifstream ifs(inputFile);
  bool inputFileExist = ifs.is_open();
  if (!inputFileExist)
  {
    std::cout << inputFile << " cannot open." << std::endl;
    return -1;
  }

  if (mode == "iv")
  {
    std::vector<uint8_t> text = load_bwt(inputFile);
    sdsl::int_vector<> bwt;
    bwt.width(8);
    bwt.resize(text.size());
    for (uint64_t i = 0; i < text.size(); i++)
    {
      bwt[i] = text[i];
    }
    sdsl::store_to_file(bwt, inputFile + ".iv");
    std::cout << "Finished." << std::endl;
    return 0;
  }
  else if (mode == "wt")
  {
    wt_huff<> wt;
    construct(wt, inputFile);
        std::cout << "WT using memory = " << sdsl::size_in_bytes(wt) / 1000 << "[KB]" << std::endl;

    std::cout << "Finished." << std::endl;
    return 0;
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
  computeMaximalSubstrings(inputFile, outputFile, true);
  //computeLCPIntervals(inputFile, true);
}
