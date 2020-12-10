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
//#include "../postorder_maximal_substrings.hpp"
#include "../main/common.hpp"
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp_dac.hpp>
#include <sdsl/lcp_support_sada.hpp>
#include "../postorder_suffix_tree_intervals.hpp"
#include "../test/naive_algorithms.hpp"
#include "../postorder_maximal_substring_intervals.hpp"

using namespace std;
using CHAR = uint8_t;
using INDEXTYPE = int64_t;
/*
std::vector<stool::LCPInterval<uint64_t>> iterateLCPIntervalsWithSDSL(std::string filename){
    sdsl::csa_sada<> csa;
    construct(csa, filename, 1);

    sdsl::lcp_dac<> dac;
    construct(dac, filename, 1);

    std::vector<stool::LCPInterval<uint64_t>> r;
    stool::esaxx::PostorderSuffixTreeIntervals<uint64_t, sdsl::csa_sada<>, sdsl::lcp_dac<>> gst;
    gst.set(std::move(csa), std::move(dac) );
    for(auto it : gst){
      r.push_back(it);
    }

    return r;
}
*/
/*
std::vector<stool::LCPInterval<uint64_t>> iterateMSIntervalsWithSDSL(std::string filename){
    sdsl::csa_sada<> csa;
    construct(csa, filename, 1);

    sdsl::lcp_dac<> dac;
    construct(dac, filename, 1);

    using CHAR = uint8_t;
    vector<CHAR> T = stool::load_text_from_file(filename, true); // input text
    assert(T.size() == csa.size());
        std::cout << "size = " << T.size() << std::endl;


    stool::DistinctCharacterCheckerOnInterval dcc;
    dcc.construct<CHAR, vector<CHAR>, sdsl::csa_sada<>>(T, csa);
  
    using ST = stool::esaxx::PostorderSuffixTreeIntervals<uint64_t, sdsl::csa_sada<>, sdsl::lcp_dac<>>;
    using STIT = stool::esaxx::PostorderSTIntervalIterator<uint64_t, sdsl::csa_sada<>, sdsl::lcp_dac<>>;

    ST gst;
    gst.set(std::move(csa), std::move(dac) );

    //stool::PostorderMaximalSubstringIntervals<uint64_t, ST, stool::DistinctCharacterCheckerOnInterval> msi;
    //msi.set(std::move(gst), std::move(dcc) );

    std::vector<stool::LCPInterval<uint64_t>> r;
 
    //auto ppp = msi.begin();

    STIT beg = gst.begin();
    STIT end = gst.end();

    stool::PostorderMaximalSubstringIntervalIterator<uint64_t, STIT, stool::DistinctCharacterCheckerOnInterval> ms_begin(beg, dcc);
    stool::PostorderMaximalSubstringIntervalIterator<uint64_t, STIT, stool::DistinctCharacterCheckerOnInterval> ms_end(end, dcc);
    
    while(ms_begin != ms_end ){
      auto it = *ms_begin;
      r.push_back(it);
      ++ms_begin;
    }
    return r;
}
*/
/*
std::vector<stool::LCPInterval<uint64_t>> iterateLCPIntervals(std::string filename){
    using INDEX = uint64_t;
    using CHAR = uint8_t;
    vector<CHAR> T = stool::load_text_from_file(filename, true); // input text
    std::cout << "Constructing Suffix Array" << std::endl;
    std::vector<INDEX> sa = stool::construct_suffix_array(T);
    std::cout << "Constructing LCP Array" << std::endl;
    std::vector<INDEX> lcpArray = stool::constructLCP<CHAR, INDEX>(T, sa);
    std::cout << "Constructing PST" << std::endl;

    stool::esaxx::PostorderSuffixTree<> pst(std::move(sa), std::move(lcpArray) );
    std::vector<stool::LCPInterval<uint64_t>> r;

    for(auto p = pst.begin(); p != pst.end();++p){
      r.push_back(*p);      
    }
    
    return r;
}
*/
/*
std::vector<stool::LCPInterval<uint64_t>> iterateMSIntervals(std::string filename){
    using INDEX = uint64_t;
    using CHAR = uint8_t;

    vector<CHAR> T = stool::load_text_from_file(filename, true); // input text
    std::cout << "size = " << T.size() << std::endl;
    std::vector<INDEX> sa = stool::construct_suffix_array(T);
    std::vector<INDEX> lcpArray = stool::constructLCP<CHAR, INDEX>(T, sa);


    stool::DistinctCharacterCheckerOnInterval dcc;
    dcc.construct<CHAR, vector<CHAR>, std::vector<INDEX> >(T, sa);
  
    using ST = stool::esaxx::PostorderSuffixTreeIntervals<INDEX, std::vector<INDEX>, std::vector<INDEX> >;
    using STIT = stool::esaxx::PostorderSTIntervalIterator<INDEX, std::vector<INDEX>, std::vector<INDEX> >;

    ST gst;
    gst.set(std::move(sa), std::move(lcpArray) );

    stool::PostorderMaximalSubstringIntervals<INDEX, ST, stool::DistinctCharacterCheckerOnInterval> msi;
    msi.set(std::move(gst), std::move(dcc) );

    std::vector<stool::LCPInterval<uint64_t>> r;
    
    for(auto it : msi){
      r.push_back(it);
    }
    

    return r;
}
*/

int main(int argc, char *argv[])
{

  cmdline::parser p;
  p.add<string>("input_file", 'i', "input file name", true);
  //p.add<string>("output_file", 'o', "output file name", false, "");
  //p.add<bool>("print", 'p', "print info", false, true);
  //p.add<string>("format", 'f', "output format (binary or csv)", false, "binary");

  p.parse_check(argc, argv);
  string inputFile = p.get<string>("input_file");
  

  
  
  /*
  auto msVec2 = iterateMSIntervals(inputFile);
  std::cout << msVec2.size() << std::endl;
  for(auto it : msVec2){
    std::cout << it.to_string() << std::endl;
  }
  */
  


}