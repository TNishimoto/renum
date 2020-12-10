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
#include "../test/esa.hxx"
#include "stool/src/io.hpp"
#include "stool/src/sa_bwt_lcp.hpp"



using namespace std;
using INDEXTYPE = int64_t;

int main(int argc, char *argv[])
{

  cmdline::parser p;
  p.add<string>("input_file", 'i', "input file name", true);
  p.add<string>("output_file", 'o', "output file name", false, "");

  //p.add<bool>("print", 'p', "print info", false, true);

  p.parse_check(argc, argv);
  string inputFile = p.get<string>("input_file");
  string outputFile = p.get<string>("output_file");
  //bool isPrint = p.get<bool>("print");

  if (outputFile.size() == 0)
  {
    outputFile = inputFile + ".interval";
  }

  vector<char> T = stool::load_text(inputFile); // input text
  INDEXTYPE n = T.size();

  vector<INDEXTYPE> SA(n); // suffix array
  vector<INDEXTYPE> L(n);  // left boundaries of internal node
  vector<INDEXTYPE> R(n);  // right boundaries of internal node
  vector<INDEXTYPE> D(n);  // depths of internal node

  INDEXTYPE alphaSize = 0x100; // This can be very large
  INDEXTYPE nodeNum = 0;

  // Computing internal nodes of the suffix tree of the input file.
  if (esaxx(T.begin(), SA.begin(),
            L.begin(), R.begin(), D.begin(),
            n, alphaSize, nodeNum) == -1)
  {
    return -1;
  }

  //INDEXTYPE size = T.size();
  /*
  if (isPrint)
  {
    std::cout << "The internal nodes of the suffix tree of the file" << std::endl;
    std::cout << "id"
              << "\t\t"
              << "occurrence"
              << "\t"
              << "range(SA)"
              << "\t"
              << "string length"
              << "\t"
              << "string" << std::endl;
  }
  */

  vector<stool::LCPInterval<INDEXTYPE>> buffer;
  ofstream os(outputFile, ios::out | ios::binary);
  if (!os)
    return 1;
  INDEXTYPE nodeCount = 0;

  // Writing and Printing internal nodes.
  
  for (INDEXTYPE i = 0; i < nodeNum; ++i)
  {
    stool::LCPInterval<INDEXTYPE> interval(L[i], R[i], D[i]);
    //INDEXTYPE len = D[i];
    nodeCount++;
    buffer.push_back(interval);
    if (buffer.size() > 8192)
    {
      os.write((const char *)(&buffer[0]), sizeof(stool::LCPInterval<INDEXTYPE>) * buffer.size());
      buffer.clear();
    }
    /*
    if (isPrint)
    {
      if (nodeCount < 1000)
      {
       std::cout << interval.getCSVLine(i, T, SA) << std::endl;
      }
      else if (nodeCount == 1000)
      {
        std::cout << "etc.." << std::endl;
      }
    }
    */
  }
  
  os.write((const char *)(&buffer[0]), sizeof(stool::LCPInterval<INDEXTYPE>) * buffer.size());
  buffer.clear();
  os.close();

  std::cout << "\033[36m";
  std::cout << "___________RESULT___________" << std::endl;
  std::cout << "File: " << inputFile << std::endl;
  std::cout << "Output: " << outputFile << std::endl;
  std::cout << "The length of the input text: " << T.size() << std::endl;
  std::cout << "The number of the internal nodes of the suffix tree of the input file: " << nodeCount << std::endl;
  std::cout << "_________________________________" << std::endl;
  std::cout << "\033[39m" << std::endl;
}