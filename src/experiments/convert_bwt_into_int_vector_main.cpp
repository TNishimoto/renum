#include <cassert>
#include <chrono>
#include "../module/stool/src/io.hpp"
#include "../module/stool/src/sa_bwt_lcp.hpp"

#include "../module/stool/src/print.hpp"
#include "../module/stool/src/cmdline.h"
#include "../module/stool/src/debug.hpp"
#include "../module/libdivsufsort/sa.hpp"

//#include "hpp/bwt.hpp"
#include "../basic/fmindex.hpp"
#include "../beller/beller_interval.hpp"
#include "../debug/beller_debug.hpp"


#include "../debug/naive_algorithms.hpp"
#include "../stnode_enumerator/single/single_stnode_traverser.hpp"
#include "../stnode_enumerator/application.hpp"

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
int main(int argc, char *argv[])
{
    cmdline::parser p;
    p.add<string>("input_file", 'i', "input bwt file path (text format)", true);
    p.add<string>("output_file", 'o', "output bwt file path (binary file) (default: input_bwt_file.iv)", false, "");

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
        outputFile = inputFile + ".iv";
    }
    std::vector<uint8_t> text = load_bwt(inputFile);

    sdsl::int_vector<> bwt;
    bwt.width(8);
    bwt.resize(text.size());

    for (uint64_t i = 0; i < text.size(); i++)
    {
        bwt[i] = text[i];
    }
    sdsl::store_to_file(bwt, outputFile);
    std::cout << "Finished." << std::endl;

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "Input BWT File \t\t\t\t\t : " << inputFile << std::endl;
    std::cout << "Output BWT File \t\t\t\t : " << outputFile << std::endl;
    std::cout << "_______________________________________________________" << std::endl;
    std::cout << "\033[39m" << std::endl;

    return 0;

}
