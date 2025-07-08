#include <cassert>
#include <chrono>
#include <stdio.h>
#include "stool/include/io.hpp"
#include "stool/include/sa_bwt_lcp.hpp"

#include "stool/include/print.hpp"
#include "stool/include/cmdline.h"
#include "stool/include/debug.hpp"
#include "../module/libdivsufsort/sa.hpp"

//#include "hpp/bwt.hpp"
#include "../include/basic/interval_search_data_structure.hpp"
//#include "../beller/beller_interval.hpp"
#include "../include/debug/beller_debug.hpp"


#include "../include/debug/naive_algorithms.hpp"
#include "../include/stnode_enumerator/single/single_stnode_traverser.hpp"
#include "../include/stnode_enumerator/application.hpp"

#include <sdsl/wt_algorithm.hpp>
#include "../include/stnode_enumerator/depth_first_search/dfs_traverser.hpp"
#include "../include/basic/sdsl_functions.hpp"

//#include "../postorder_maximal_substring_intervals.hpp"
//#include "../forward_bwt.hpp"

using namespace std;
//using namespace stool;
//using namespace stool::rlbwt;

using CHAR = char;
using INDEX = uint64_t;
using LCPINTV = stool::LCPInterval<uint64_t>;

class DFSApplication
{
public:
    template <typename INDEX_SIZE, typename INTERVAL_SEARCH, typename CHAR = uint8_t>
    static uint64_t outputMaximalSubstrings(std::ofstream &out, stool::renum::DFSTraverser<INDEX_SIZE, INTERVAL_SEARCH, CHAR> &stnodeSequencer, stool::renum::STTreeAnalysisResult &analysis)
    {

        //using INDEX_SIZE = typename stool::renum::DFSTraverser<INDEX_SIZE, INTERVAL_SEARCH, CHAR>::index_type;
        uint8_t print_type = 0;
        out.write((char *)(&print_type), sizeof(print_type));
        uint8_t index_bits = sizeof(INDEX_SIZE);
        out.write((char *)(&index_bits), sizeof(index_bits));

        std::vector<stool::LCPInterval<INDEX_SIZE>> buffer;
        uint64_t count = 0;

        analysis.start(stnodeSequencer.get_input_text_length());

        auto it = stnodeSequencer.begin();
        while (it != stnodeSequencer.end())
        {
            analysis.analyze(it.child_count(), stnodeSequencer.get_stack_size());
            bool b = it.is_maximal_repeat();
            if (b)
            {
                uint64_t left = it.get_left();
                uint64_t right = it.get_right();
                uint64_t lcp = it.get_lcp();

                stool::LCPInterval<INDEX_SIZE> intv(left, right, lcp);
                buffer.push_back(intv);
                count++;
                if (buffer.size() >= 8000)
                {
                    out.write((char *)(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                    buffer.clear();
                }
            }

            it++;
        }
        if (buffer.size() >= 1)
        {
            out.write((char *)(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
            buffer.clear();
        }
        std::cout << "Enumerated" << std::endl;

        std::cout
            << "STOP" << std::endl;
        return count;
        /*
        analysis.start();

        using INDEX_SIZE = typename STNODES::index_type;
        uint8_t print_type = 0;
        out.write((char *)(&print_type), sizeof(print_type));
        uint8_t index_bits = sizeof(INDEX_SIZE);
        out.write((char *)(&index_bits), sizeof(index_bits));

        uint64_t count = 0;

        analysis.input_text_length = stnodeSequencer.get_input_text_length();

        std::vector<stool::LCPInterval<INDEX_SIZE>> buffer;

        auto depth_iterator = stnodeSequencer.begin();
        while (depth_iterator != stnodeSequencer.end())
        {
            analysis.analyze(stnodeSequencer);
            std::vector<stool::LCPInterval<uint64_t>> r2;

            for (auto node_it = depth_iterator.begin(); node_it != depth_iterator.end(); node_it++)
            {
                bool b = stnodeSequencer.check_maximal_repeat(node_it);
                if (b)
                {
                    INDEX_SIZE left = stnodeSequencer.get_left(node_it);
                    INDEX_SIZE right = stnodeSequencer.get_right(node_it);
                    stool::LCPInterval<INDEX_SIZE> intv(left, right, stnodeSequencer.get_current_lcp());
                    buffer.push_back(intv);
                    count++;
                    if (buffer.size() >= 8000)
                    {
                        out.write((char *)(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                        buffer.clear();
                    }
                }
            }
            depth_iterator++;
        }

        if (buffer.size() >= 1)
        {
            out.write((char *)(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
            buffer.clear();
        }
        std::cout << "Enumerated" << std::endl;

        return count;
        */
    }
};
/*
int deleteFile(string fileName)
{
    return !(remove(fileName));
}

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
    if(bwt.size() < 200){
        for(uint64_t i=0;i<bwt.size();i++){
        std::cout << bwt[i] << ", " << std::flush;
        }
        std::cout << std::endl;
    }
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
    std::cout << "Constructing array C..." << std::endl;

    stool::renum::FMIndex::constructC(bwt, C);

    return bwt[bwt.size() - 1];
}
*/
void computeMaximalSubstrings(std::string inputFile, std::string outputFile)
{

    //string text = "";
    auto start = std::chrono::system_clock::now();
    //sdsl::bit_vector bv;

    std::cout << "Loading : " << inputFile << std::endl;
    //std::cout << "Constructing dbit array for maximal repeats : " << std::endl;
    //stool::renum::SDSLFunction::constructDBitArray(inputFile, bv);
    //sdsl::bit_vector::rank_1_type bwt_bit_rank1(&bv);



    std::cout << "Constructing Wavelet Tree..." << std::endl;
    sdsl::wt_huff<> wt;
    std::vector<uint64_t> C;
    uint8_t lastChar = stool::renum::SDSLFunction::load_wavelet_tree(inputFile, wt, C);








    uint64_t ms_count = 0;

    stool::renum::IntervalSearchDataStructure<uint8_t> range;
    range.initialize(&wt, &C, lastChar);

    std::ofstream out(outputFile, std::ios::out | std::ios::binary);
    if (!out)
    {
        throw std::runtime_error("Cannot open the output file!");
    }
    uint64_t input_text_size = wt.size();

    auto mid = std::chrono::system_clock::now();

    std::cout << "Enumerating..." << std::endl;
    stool::renum::STTreeAnalysisResult st_result;

    if (input_text_size - 10 < UINT32_MAX)
    {
        using INDEX_TYPE = uint32_t;

        stool::renum::ExplicitWeinerLinkComputer<INDEX_TYPE> wsearch;
        wsearch.initialize(&range, input_text_size);
        stool::renum::DFSTraverser<INDEX_TYPE, stool::renum::ExplicitWeinerLinkComputer<INDEX_TYPE>> traverser;
        traverser.initialize(&wsearch, false);
        ms_count = DFSApplication::outputMaximalSubstrings(out, traverser, st_result);
    }
    else
    {
        using INDEX_TYPE = uint64_t;

        stool::renum::ExplicitWeinerLinkComputer<INDEX_TYPE> wsearch;
        wsearch.initialize(&range, input_text_size);
        stool::renum::DFSTraverser<INDEX_TYPE, stool::renum::ExplicitWeinerLinkComputer<INDEX_TYPE>> traverser;
        traverser.initialize(&wsearch, false);
        ms_count = DFSApplication::outputMaximalSubstrings(out, traverser, st_result);
    }
    auto end = std::chrono::system_clock::now();
    double enumeration_time = std::chrono::duration_cast<std::chrono::seconds>(end - mid).count();
    double construction_time = std::chrono::duration_cast<std::chrono::seconds>(mid - start).count();
    //std::cout << "Construction time: " << construction_time << "[s]" << std::endl;

    double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    double bps = ((double)input_text_size / ((double)elapsed)) / 1000;

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "BWT File \t\t\t\t : " << inputFile << std::endl;
    std::cout << "Output File \t\t\t\t : " << outputFile << std::endl;
    std::cout << "The length of the input text \t\t : " << input_text_size << std::endl;
    std::cout << "The number of maximum substrings \t : " << ms_count << std::endl;
    std::cout << "The usage of Wavelet tree \t\t : " << sdsl::size_in_bytes(wt) / 1000 << "[KB]" << std::endl;

    //std::cout << "The usage of DBit array \t\t : " << sdsl::size_in_bytes(bv) / 1000 << "[KB]" << std::endl;

    std::cout << "Excecution time \t\t\t : " << elapsed << "[s]" << std::endl;
    std::cout << "Character per second \t\t\t : " << bps << "[KB/s]" << std::endl;

    std::cout << "\t Preprocessing time \t\t : " << construction_time << "[s]" << std::endl;
    std::cout << "\t Enumeration time \t\t : " << enumeration_time << "[s]" << std::endl;

    std::cout << "_______________________________________________________" << std::endl;
    std::cout << "\033[39m" << std::endl;

}


int main(int argc, char *argv[])
{
    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file name", true);
    p.add<string>("output_file", 'o', "output file path (default: input_file_path.dfs.max)", false, "");

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    string outputFile = p.get<string>("output_file");

    std::ifstream ifs(inputFile);
    bool inputFileExist = ifs.is_open();
    if (!inputFileExist)
    {
        std::cout << inputFile << " cannot open." << std::endl;
        return -1;
    }

    if (outputFile.size() == 0)
    {
        outputFile = inputFile + ".dfs.max";
    }
    computeMaximalSubstrings(inputFile, outputFile);
}
