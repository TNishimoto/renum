// License: MIT http://opensource.org/licenses/MIT

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "stool/include/third_party/cmdline.h"
#include "stool/include/strings/sa_bwt_lcp.hpp"
#include "stool/include/strings/lcp_interval.hpp"

#include "../include/rlbwt/bwt_analysis_result.hpp"
#include "../include/rlbwt/rle.hpp"
#include "../include/rlbwt/bwt_decompress.hpp"

#include "../include/basic/fmindex.hpp"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wt_gmr.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>

//#include "divsufsort.h"
//#include "divsufsort64.h"
//#include "stool/src/io.hpp"
//#include "libdivsufsort/sa.hpp"
//#include "common.hpp"

//#include "../minimal_substrings/naive_minimal_substrings.hpp"

//using namespace std;
//using namespace stool;

using INDEXTYPE = int64_t;
using CHAR = char;

class BWTDecompressor
{
public:
    sdsl::int_vector<> bwt;
    sdsl::wt_gmr<> wt;
    std::vector<uint64_t> C;
    void initialize(std::string inputFile)
    {
        std::vector<uint8_t> _bwt;
        stool::bwt::load(inputFile, _bwt);
        this->bwt.resize(_bwt.size());
        for (uint64_t i = 0; i < _bwt.size(); i++)
        {
            this->bwt[i] = _bwt[i];
        }
        stool::renum::FMIndex::constructSelect(this->bwt, wt);
        stool::renum::FMIndex::constructC(this->bwt, this->C);
    }
    void extract(uint64_t i, uint64_t l, std::string &output)
    {
        uint64_t p = i;
        output.resize(l);
        for (uint64_t x = 0; x < l; x++)
        {
            p = stool::renum::FMIndex::FL(p, C, wt);
            output[x] = this->bwt[p];
        }
    }
};

struct LCPIntervalFileInfo
{
    uint8_t print_type;
    uint8_t index_bits;
    uint64_t count;

    void print()
    {
        if (print_type == 0)
        {
            std::cout << "File Type: "
                      << "LCPInterval" << std::endl;
        }
        else
        {
            std::cout << "File Type: "
                      << "???" << std::endl;
        }
        std::cout << "Integer Type: "
                  << "uint" << (8 * index_bits) << "t" << std::endl;
        std::cout << "Count: " << count << std::endl;
    }
};

LCPIntervalFileInfo load_lcp_interval_info(std::ifstream &inp)
{
    LCPIntervalFileInfo info;
    uint8_t print_type = 0;
    inp.read((char *)&print_type, sizeof(uint8_t));
    uint8_t index_bits = 0;
    inp.read((char *)&index_bits, sizeof(uint8_t));
    info.print_type = print_type;
    info.index_bits = index_bits;

    uint64_t len;
    inp.seekg(0, std::ios::end);
    uint64_t n = (unsigned long)inp.tellg();
    inp.seekg(2, std::ios::beg);
    len = (n - 2) / (info.index_bits * 3);
    info.count = len;
    return info;
}

template <typename INDEX>
void load_lcp_intervals(std::ifstream &inp, std::vector<stool::LCPInterval<INDEX>> &output)
{

    uint64_t len;
    inp.seekg(0, std::ios::end);
    uint64_t n = (unsigned long)inp.tellg();
    inp.seekg(2, std::ios::beg);
    len = (n - 2) / (sizeof(stool::LCPInterval<INDEX>));
    output.resize(len);
    inp.read((char *)&output[0], sizeof(stool::LCPInterval<INDEX>) * len);
}

int main(int argc, char *argv[])
{

    cmdline::parser p;
    p.add<std::string>("input_file", 'i', "input bwt file path", true);
    p.add<std::string>("lcp_interval_file", 'l', "LCP interval file path", true);
    p.add<std::string>("output_file", 'o', "output file path (default: lcp_interval_file.interval.log)", false, "");
    p.add<bool>("string_flag", 's', "Output the string represented by each interval if this flag is 1", false, true);

    //p.add<string>("tree_file", 't', "file type", false, "NULL");

    p.parse_check(argc, argv);
    std::string inputFile = p.get<std::string>("input_file");
    std::string intervalFile = p.get<std::string>("lcp_interval_file");
    //string type = p.get<string>("tree_file");
    std::string outputFile = p.get<std::string>("output_file");
    bool string_flag = p.get<bool>("string_flag");

    if (outputFile.size() == 0)
    {
        outputFile = intervalFile + ".interval.log";
    }
    std::ofstream out(outputFile, std::ios::out);
    if (!out)
    {
        throw std::runtime_error("Cannot open the output file!");
    }

    std::ifstream inp;
    inp.open(intervalFile, std::ios::binary);
    bool inputFileExist = inp.is_open();
    if (!inputFileExist)
    {
        std::cout << intervalFile << " cannot open." << std::endl;
        throw std::runtime_error("error");
    }

    /*
    stool::rlbwt2::BWTAnalysisResult analysisResult;
    stool::renum::RLE<uint8_t> rlbwt;
    rlbwt.load(inputFile, analysisResult);
    */
    BWTDecompressor dec;
    dec.initialize(inputFile);

    auto info = load_lcp_interval_info(inp);
    info.print();
    if (string_flag)
    {
        out << "(i, j, LCP, substring)" << std::endl;
    }
    else
    {
        out << "(i, j, LCP)" << std::endl;
    }

    std::string tmpStr;
    if (info.index_bits == 4)
    {
        std::vector<stool::LCPInterval<uint32_t>> lcp_intervals;
        load_lcp_intervals(inp, lcp_intervals);
        inp.close();
        if (lcp_intervals.size() < 100)
        {
            if (string_flag)
            {
                std::cout << "(i, j, LCP, substring)" << std::endl;
            }
            else
            {
                std::cout << "(i, j, LCP)" << std::endl;
            }
        }

        for (auto &it : lcp_intervals)
        {
            out << it.i << ", " << it.j << ", " << it.lcp;
            if (lcp_intervals.size() < 100)
            {
                std::cout << it.i << ", " << it.j << ", " << it.lcp;
            }
            if (string_flag)
            {
                dec.extract(it.i, it.lcp, tmpStr);
                out << ", " << tmpStr;
                if (lcp_intervals.size() < 100)
                {
                    std::cout << ", " << tmpStr;
                }
            }
            out << std::endl;
            if (lcp_intervals.size() < 100)
            {
                std::cout << std::endl;
            }
        }
    }
    else if (info.index_bits == 8)
    {
        std::vector<stool::LCPInterval<uint64_t>> lcp_intervals;
        load_lcp_intervals(inp, lcp_intervals);
        inp.close();

        for (auto &it : lcp_intervals)
        {
            out << it.i << ", " << it.j << ", " << it.lcp;
            if (string_flag)
            {
                dec.extract(it.i, it.lcp, tmpStr);
                out << ", " << tmpStr;
            }
            out << std::endl;
        }
    }
    else
    {
        assert(false);
        throw -1;
    }

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "BWT File \t\t\t\t : " << inputFile << std::endl;
    std::cout << "Interval File \t\t\t\t : " << intervalFile << std::endl;
    std::cout << "Output File \t\t\t\t : " << outputFile << std::endl;
    info.print();

    std::cout << "_______________________________________________________" << std::endl;

    std::cout << "\033[39m" << std::endl;

    return 0;
}