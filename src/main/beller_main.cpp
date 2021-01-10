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
#include "../debug/beller_debug.hpp"

#include "../main/common.hpp"
#include "../debug/naive_algorithms.hpp"
#include <sdsl/wt_algorithm.hpp>

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
uint8_t get_last_char(std::string inputFile, std::vector<uint64_t> &C, sdsl::bit_vector &bv)
{
  sdsl::int_vector<> bwt;
  sdsl::load_from_file(bwt, inputFile);
  uint64_t k = 0;
  bv.resize(bwt.size());
  bool b = true;
  //std::cout << bwt.size() << std::endl;
  for (uint64_t i = 0; i < bwt.size(); i++)
  {
    //std::cout << (uint) bwt[i];
    if (bwt[i] == 0)
    {
      k++;
      std::cout << i << std::endl;
    }
    if (i > 0 && bwt[i] != bwt[i - 1])
    {
      b = !b;
    }
    bv[i] = b;
  }
  //std::cout << std::endl;
  if (k == 0)
  {
    std::cout << "Error: This bwt does not contain 0." << std::endl;
    throw -1;
  }
  else if (k >= 2)
  {
    std::cout << "Error2: This bwt contains 0 twice or more." << std::endl;
    throw -1;
  }
  std::cout << "Constructing array C..." << std::endl;

  stool::FMIndex::constructC(bwt, C);

  return bwt[bwt.size() - 1];
}

void computeLCPIntervals(std::string inputFile, bool correctCheck)
{

  //string text = "";
  std::cout << "Loading : " << inputFile << std::endl;
  std::vector<uint8_t> text = stool::load_text_from_file(inputFile, true);
  vector<INDEX> sa = stool::construct_suffix_array(text);
  sdsl::int_vector<> bwt;
  stool::FMIndex::constructBWT(text, sa, bwt);

  std::vector<uint64_t> C;
  stool::FMIndex::constructC(bwt, C);

  wt_huff<> wt;
  construct_im(wt, bwt);

  uint64_t lastChar = bwt[bwt.size() - 1];

  stool::IntervalSearchDataStructure range;
  range.initialize(&wt, &C, lastChar);
  /*
  std::cout << wt.size() << "/" << bwt.size() << std::endl;
  for(uint64_t i=0;i<bwt.size() - 1;i++){
    uint8_t c1 = wt[i+1];
    uint8_t c2 = bwt[i];
    if(c1 != c2){
    std::cout << c1 << "/" << c2 << std::endl;

    }


  }
  */

  auto test_Intervals = stool::beller::computeLCPIntervals<uint64_t>(range);
  test_Intervals.push_back(LCPINTV(0, text.size() - 1, 0));

  if (correctCheck)
  {
    auto correctLCP = stool::constructLCP(text, sa);
    std::cout << "Correct" << std::endl;
    std::vector<LCPINTV> correct_intervals = stool::esaxx::NaiveAlgorithms::naive_compute_lcp_intervals(text, sa);
    stool::beller::equal_check_lcp_intervals(test_Intervals, correct_intervals);
    std::cout << "OK!" << std::endl;
  }
}

void computeMaximalSubstrings(std::string inputFile, std::string outputFile, bool correctCheck)
{

  //string text = "";
  auto start = std::chrono::system_clock::now();
  std::vector<uint64_t> C;
  sdsl::bit_vector bv;

  std::cout << "Loading : " << inputFile << std::endl;
  uint8_t lastChar = get_last_char(inputFile, C, bv);
  sdsl::bit_vector::rank_1_type bwt_bit_rank1(&bv);

  std::cout << "Constructing Wavelet Tree..." << std::endl;
  wt_huff<> wt;
  construct(wt, inputFile);
  std::cout << "WT using memory = " << sdsl::size_in_bytes(wt) / 1000 << "[KB]" << std::endl;

  uint64_t ms_count = 0;

  stool::IntervalSearchDataStructure range;
  range.initialize(&wt, &C, lastChar);

  std::ofstream out(outputFile, std::ios::out | std::ios::binary);
  if (!out)
  {
    throw std::runtime_error("Cannot open the output file!");
  }
  uint64_t input_text_size = wt.size();

  auto mid = std::chrono::system_clock::now();
  double construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start).count();
  std::cout << "Construction time: " << construction_time << "[ms]" << std::endl;

  std::cout << "Enumerating..." << std::endl;
  uint64_t peak_count = 0;

  if (input_text_size - 10 < UINT32_MAX)
  {
    stool::beller::OutputStructure<uint32_t> os;
    os.is_output_maximal_substrings = true;
    os.out = &out;
    os.bwt_bit_rank1 = &bwt_bit_rank1;

    ms_count = stool::beller::outputMaximalSubstrings<uint32_t>(range, os);
    peak_count = os.peak;
  }
  else
  {
    stool::beller::OutputStructure<uint64_t> os;
    os.is_output_maximal_substrings = true;
    os.out = &out;
    os.bwt_bit_rank1 = &bwt_bit_rank1;

    ms_count = stool::beller::outputMaximalSubstrings<uint64_t>(range, os);
    peak_count = os.peak;

  }
  auto end = std::chrono::system_clock::now();
  double enumeration_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - mid).count();

  double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  double bps = ((double)input_text_size / ((double)elapsed / 1000)) / 1000;

  std::cout << "\033[31m";
  std::cout << "______________________RESULT______________________" << std::endl;
  std::cout << "RLBWT File \t\t\t\t\t : " << inputFile << std::endl;
  std::cout << "Output \t\t\t\t\t : " << outputFile << std::endl;
  std::cout << "The length of the input text \t\t : " << input_text_size << std::endl;
  std::cout << "The number of maximum substrings \t : " << ms_count << std::endl;
  std::cout << "Peak count \t : " << peak_count << std::endl;
  std::cout << "The usage of Wavelet tree : " << sdsl::size_in_bytes(wt) / 1000 << "[KB]" << std::endl;

  std::cout << "Excecution time \t\t\t : " << elapsed << "[ms]" << std::endl;
  std::cout << "Character per second \t\t\t : " << bps << "[KB/s]" << std::endl;

  std::cout << "\t Preprocessing time \t\t : " << construction_time << "[ms]" << std::endl;
  std::cout << "\t Enumeration time \t\t : " << enumeration_time << "[ms]" << std::endl;

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

    uint64_t k = 0;
    for (uint64_t i = 0; i < text.size(); i++)
    {
      bwt[i] = text[i];
      if (bwt[i] == 0)
      {
        k++;
      }
    }
    if (k != 1)
    {
      std::cout << "bwt error" << std::endl;
      throw -1;
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
  else if (mode == "test")
  {
    computeLCPIntervals(inputFile, true);

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
}
