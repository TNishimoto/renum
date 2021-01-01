#include <cassert>
#include <chrono>
#include <string>

#include "stool/src/io.hpp"
#include "stool/src/cmdline.h"
#include "stool/src/debug.hpp"
#include "../stnode_enumerator/single/bit_deque.hpp"



int main(int argc, char *argv[])
{
    /*
  cmdline::parser p;
  p.add<std::string>("input_file", 'i', "input file name", true);
  p.add<std::string>("mode", 'm', "mode", false, "xx");
  p.add<std::string>("output_file", 'o', "output file name", false, "");

  p.parse_check(argc, argv);
  std::string inputFile = p.get<std::string>("input_file");
  std::string mode = p.get<std::string>("mode");
  std::string outputFile = p.get<std::string>("output_file");
  std::string format = "binary";
  */

    stool::lcp_on_rlbwt::BitDeque bdeq;

    std::cout << "Test[push_back]" << std::endl;
    for(uint64_t i=0;i<20;i++){
        bdeq.push_back(i % 2 == 0);
        auto r = bdeq.to_vector();
        stool::Printer::print(r);        
    }
    std::cout << "Test[pop_back]" << std::endl;
    for(uint64_t i=0;i<20;i++){
        bdeq.pop_back();
        auto r = bdeq.to_vector();
        stool::Printer::print(r);        
    }
    std::cout << "Test[push_front]" << std::endl;
    for(uint64_t i=0;i<20;i++){
        bdeq.push_front(i % 2 == 0);
        auto r = bdeq.to_vector();
        stool::Printer::print(r);
    }
    std::cout << "Test[pop_front]" << std::endl;
    for(uint64_t i=0;i<20;i++){
        bdeq.pop_front();
        auto r = bdeq.to_vector();
        stool::Printer::print(r);        
    }

    /*
    std::cout << "Test[pop_back]" << std::endl;
    for(uint64_t i=0;i<20;i++){
        bdeq.pop_back();
        auto r = bdeq.to_vector();
        stool::Printer::print(r);        
    }
    */

}
