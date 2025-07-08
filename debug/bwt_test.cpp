#include <cassert>
#include <chrono>
#include <stdio.h>
#include "stool/include/stool.hpp"
#include "../module/libdivsufsort/sa.hpp"

//#include "hpp/bwt.hpp"
#include "../basic/interval_search_data_structure.hpp"
//#include "../beller/beller_interval.hpp"
#include "../debug/beller_debug.hpp"

#include "../debug/naive_algorithms.hpp"
#include "../stnode_enumerator/single/single_stnode_traverser.hpp"
#include "../stnode_enumerator/application.hpp"

#include <sdsl/wt_algorithm.hpp>
#include "../stnode_enumerator/depth_first_search/dfs_traverser.hpp"
#include "../basic/sdsl_functions.hpp"

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
void testCArray(std::string filename){
    sdsl::int_vector<> iv;
    sdsl::load_from_file(iv, filename);
    std::vector<uint64_t> C;
    stool::renum::FMIndex::constructC(iv, C);

    uint8_t c = stool::renum::SDSLFunction::get_last_char(filename);
    std::cout << (int)c << std::endl;

    std::cout << "Constructing Wavelet Tree..." << std::endl;
    wt_huff<> wt;
    construct(wt, filename);
    std::cout << "WT using memory = " << sdsl::size_in_bytes(wt) / 1000 << "[KB]" << std::endl;

    std::vector<uint64_t> C2;
    //stool::renum::FMIndex::constructC(bwt, C);
    stool::renum::FMIndex::constructCArray(wt, c, C2);

    stool::equal_check(C, C2);
}

void naiveConstructDBitArray(std::string bwt_iv_file, sdsl::bit_vector &bv){
    sdsl::int_vector<> bwt;
    sdsl::load_from_file(bwt, bwt_iv_file);
    uint64_t k = 0;
    bv.resize(bwt.size());
    bool b = true;
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
}

void testDBitArray(std::string bwt_iv_file){
    sdsl::bit_vector bv1, bv2;
    naiveConstructDBitArray(bwt_iv_file, bv1);
    stool::renum::SDSLFunction::constructDBitArray(bwt_iv_file, bv2);
    if(bv1.size() < 100){
        for(uint64_t i=0;i<bv1.size();i++){
            std::cout << bv1[i];
        }
        std::cout << std::endl;
    }
    if(bv2.size() < 100){
        for(uint64_t i=0;i<bv2.size();i++){
            std::cout << bv2[i];
        }
        std::cout << std::endl;
    }

    if(bv1.size() != bv2.size()){
        std::cout << "Different Size " << bv1.size() << "/" << bv2.size() << std::endl;
            throw -1;

    }

    for(uint64_t i=0;i<bv1.size() ; i++){
        if(bv1[i] != bv2[i]){
            std::cout << "Different Value" << std::endl;
            throw -1;

        }
    }
}
void testBWT_WT(std::string bwt_iv_file){
    sdsl::int_vector<> bwt;
    sdsl::load_from_file(bwt, bwt_iv_file);

    uint8_t c = stool::renum::SDSLFunction::get_last_char(bwt_iv_file);


    wt_huff<> wt;
    construct(wt, bwt_iv_file);
    std::cout << "Size/" << bwt.size() << "/" << wt.size() << std::endl;

    for(uint64_t x=0;x<bwt.size();x++){
        uint8_t c1 = bwt[x];
        uint8_t c2 = stool::renum::SDSLFunction::get_Char(wt, x, c);

        if(c1 != c2){
            std::cout << "Error! " << c1 << "/" << c2 << "/" << x << std::endl;
            throw -1;
        }
    }


}

int main(int argc, char *argv[])
{
    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file name", true);

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");

    std::ifstream ifs(inputFile);
    bool inputFileExist = ifs.is_open();
    if (!inputFileExist)
    {
        std::cout << inputFile << " cannot open." << std::endl;
        return -1;
    }
    std::string ivFile = inputFile + ".iv";

    auto vec = load_bwt(inputFile);
    stool::Printer::print(vec);
    throw -1;
    /*
    testDBitArray(ivFile);
    
    testCArray(ivFile);
    */
   testBWT_WT(ivFile);
    std::cout << "OK!" << std::endl;
}
