#include <cassert>
#include <chrono>
#include "../module/stool/src/io.hpp"
#include "../module/stool/src/sa_bwt_lcp.hpp"

#include "../module/stool/src/print.hpp"
#include "../module/stool/src/cmdline.h"
#include "../module/stool/src/debug.hpp"
#include "../module/libdivsufsort/sa.hpp"
//#include "../module/rlbwt_iterator/src/include/rlbwt_iterator.hpp"
//#include "module/rlbwt_iterator/src/include/bwt.hpp"

#include "stool/src/io.hpp"
#include "stool/src/cmdline.h"
#include "stool/src/debug.hpp"

#include <sdsl/bit_vectors.hpp>

#include "../debug/beller_debug.hpp"
#include "../debug/naive_algorithms.hpp"
#include "../debug/stnode_checker.hpp"
#include "../debug/hyper_debug.hpp"

#include "../stnode_enumerator/application.hpp"
#include "../stnode_enumerator/rlcp_interval_enumerator.hpp"
#include "../rlbwt/io.hpp"
#include "../rlbwt/fpos_data_structure.hpp"
#include "../rlbwt/bwt_decompress.hpp"
#include "../main/common.hpp"
#include "../stnode_enumerator/weiner_link_search.hpp"

using namespace std;
using namespace stool;

using CHAR = char;
using INDEX = uint64_t;
using LCPINTV = stool::LCPInterval<uint64_t>;

class MaximalRepeatTest
{
public:
    template <typename STNODES>
    static uint64_t test(STNODES &stnodeSequencer, std::vector<uint8_t> &bwt)
    {
        using INDEX_SIZE = typename STNODES::index_type;

        uint64_t count = 0;

        std::vector<stool::LCPInterval<INDEX_SIZE>> buffer;

        auto it = stnodeSequencer.begin();
        while (it != stnodeSequencer.end())
        {

            for (auto node_it = it.begin(); node_it != it.end(); node_it++)
            {
                bool b = stnodeSequencer.check_maximal_repeat(node_it);
                if (b)
                {
                    count++;
                    INDEX_SIZE left = stnodeSequencer.get_left(node_it);
                    INDEX_SIZE right = stnodeSequencer.get_right(node_it);

                    bool mb = false;
                    for (uint64_t i = left + 1; i <= right; i++)
                    {
                        uint8_t l = bwt[i - 1];
                        uint8_t r = bwt[i];

                        if (l != r)
                        {
                            mb = true;
                            break;
                        }
                    }
                    if (mb != b)
                    {
                        std::cout << mb << "/" << b << std::endl;
                        std::cout << "[" << left << ", " << right << "]" << std::endl;

                        std::cout << "BWT on Range: ";
                        for (uint64_t i = left; i <= right; i++)
                        {
                            std::cout << bwt[i];
                        }
                        std::cout << std::endl;
                    }
                    assert(mb == b);
                }
            }
            it++;
        }

        return count;
    }
};


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
uint8_t get_last_char(sdsl::int_vector<> &bwt, std::vector<uint64_t> &C, sdsl::bit_vector &bv)
{
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

void computeMaximalSubstrings(std::string inputFile, string mode, int thread_num)
{
    auto start = std::chrono::system_clock::now();


    sdsl::int_vector<> diff_char_vec;
    stool::EliasFanoVectorBuilder run_bits;
    //std::vector<bool> run_bits;
    auto bwtAnalysis = stool::rlbwt2::load_RLBWT_from_file(inputFile, diff_char_vec, run_bits);
    std::cout << "BWT using memory = " << sdsl::size_in_bytes(diff_char_vec) / 1000 << "[KB]" << std::endl;
    std::cout << "Run bits using memory = " << run_bits.get_using_memory() / 1000 << "[KB]" << std::endl;

    uint64_t data_structure_bytes = 0;
    sdsl::store_to_file(diff_char_vec, inputFile + ".tmp");

    stool::WT wt;
    construct(wt, inputFile + ".tmp");

    string bit_size_mode = "UINT64_t";
    data_structure_bytes += sdsl::size_in_bytes(wt);
    std::cout << "WT using memory = " << sdsl::size_in_bytes(wt) / 1000 << "[KB]" << std::endl;

    uint64_t ms_count = 0;
    stool::lcp_on_rlbwt::STTreeAnalysisResult st_result;

    double construction_time = 0;
    std::chrono::system_clock::time_point mid;

    stool::EliasFanoVector lpos_vec;
    lpos_vec.build_from_builder(run_bits);
    std::cout << "LPOS Vec using memory = " << lpos_vec.get_using_memory() / 1000 << "[KB]" << std::endl;
    data_structure_bytes += lpos_vec.get_using_memory();

    //lpos_vec.build_from_bit_vector(run_bits);
    using LPOSDS = stool::EliasFanoVector;
    using FPOSDS = stool::lcp_on_rlbwt::LightFPosDataStructure;
    FPOSDS fposds = stool::lcp_on_rlbwt::LightFPosDataStructure(diff_char_vec, lpos_vec, wt);
    std::cout << "FPOS Vec using memory = " << fposds.get_using_memory() / 1000 << "[KB]" << std::endl;
    data_structure_bytes += fposds.get_using_memory();

    mid = std::chrono::system_clock::now();
    construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start).count();
    std::cout << "Construction time: " << construction_time << "[ms]" << std::endl;
    std::cout << "Data structure Size \t\t\t : " << (data_structure_bytes / 1000) << "[KB]" << std::endl;
    using RDS = stool::lcp_on_rlbwt::RLEWaveletTree<uint32_t, LPOSDS, FPOSDS>;
    RDS ds = RDS(diff_char_vec, wt, lpos_vec, fposds);

    std::cout << "Enumerate Maximal Substrings..." << std::endl;
    stool::lcp_on_rlbwt::SuffixTreeNodes<uint32_t, RDS> stnodeTraverser;
    stnodeTraverser.use_fast_mode = mode == "1";

    stnodeTraverser.initialize(thread_num, ds);

    std::vector<uint8_t> plain_bwt;
    stool::bwt::load(inputFile, plain_bwt);
    ms_count = MaximalRepeatTest::test(stnodeTraverser, plain_bwt);
    bit_size_mode = "UINT32_t";

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double enumeration_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - mid).count();

    double bps = ((double)bwtAnalysis.str_size / ((double)elapsed / 1000)) / 1000;

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "RLBWT File \t\t\t\t : " << inputFile << std::endl;
    std::cout << "Use Fast \t\t\t\t\t : " << (stnodeTraverser.use_fast_mode ? "True" : "False") << std::endl;

    std::cout << "LPOS and FPos Vector type \t\t\t\t : "
              << "EliasFano" << std::endl;
    std::cout << "Peak children count \t\t\t : " << st_result.max_nodes_at_level << std::endl;
    std::cout << "Data structure Size \t\t\t : " << (data_structure_bytes / 1000) << "[KB]" << std::endl;

    std::cout << "Thread number \t\t\t\t : " << thread_num << std::endl;
    std::cout << "Integer Type \t\t\t\t : " << bit_size_mode << std::endl;

    std::cout << "The length of the input text \t\t : " << bwtAnalysis.str_size << std::endl;
    std::cout << "The number of runs on BWT \t\t : " << bwtAnalysis.run_count << std::endl;
    std::cout << "The number of maximum substrings \t : " << ms_count << std::endl;
    std::cout << "Excecution time \t\t\t : " << elapsed << "[ms]" << std::endl;
    std::cout << "Character per second \t\t\t : " << bps << "[KB/s]" << std::endl;
    std::cout << "\t Preprocessing time \t\t : " << construction_time << "[ms]" << std::endl;
    std::cout << "\t Enumeration time \t\t : " << enumeration_time << "[ms]" << std::endl;

    std::cout << "_______________________________________________________" << std::endl;

    std::cout << "\033[39m" << std::endl;
}

void computeMaximalSubstrings_beller(std::string inputFile)
{
    std::cout << "Loading : " << inputFile << std::endl;

    std::vector<uint8_t> plain_bwt = load_bwt(inputFile);

    sdsl::int_vector<> bwt;
    bwt.width(8);
    bwt.resize(plain_bwt.size());
    for(uint64_t i=0;i<plain_bwt.size();i++){
        bwt[i] = plain_bwt[i];
    }

    //string text = "";
    auto start = std::chrono::system_clock::now();
    std::vector<uint64_t> C;
    sdsl::bit_vector bv;

    uint8_t lastChar = get_last_char(bwt, C, bv);
    sdsl::bit_vector::rank_1_type bwt_bit_rank1(&bv);

    std::cout << "Constructing Wavelet Tree..." << std::endl;
    wt_huff<> wt;
    construct_im(wt, bwt);
    std::cout << "WT using memory = " << sdsl::size_in_bytes(wt) / 1000 << "[KB]" << std::endl;

    uint64_t ms_count = 0;

    stool::IntervalSearchDataStructure range;
    range.initialize(&wt, &C, lastChar);

    uint64_t input_text_size = wt.size();

    auto mid = std::chrono::system_clock::now();
    double construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start).count();
    std::cout << "Construction time: " << construction_time << "[ms]" << std::endl;

    std::cout << "Enumerating..." << std::endl;
    uint64_t peak_count = 0;
    stool::lcp_on_rlbwt::STTreeAnalysisResult st_result;

        using INDEX_TYPE = uint32_t;

        stool::lcp_on_rlbwt::ExplicitWeinerLinkSearch<INDEX_TYPE> wsearch;
        wsearch.initialize(&range, &bwt_bit_rank1, input_text_size);
        stool::lcp_on_rlbwt::SingleSTNodeTraverser<INDEX_TYPE, stool::lcp_on_rlbwt::ExplicitWeinerLinkSearch<INDEX_TYPE>> traverser;
        traverser.initialize(&wsearch);
        ms_count = MaximalRepeatTest::test(traverser, plain_bwt);
    auto end = std::chrono::system_clock::now();
    double enumeration_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - mid).count();

    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double bps = ((double)input_text_size / ((double)elapsed / 1000)) / 1000;

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "RLBWT File \t\t\t\t\t : " << inputFile << std::endl;
    std::cout << "The length of the input text \t\t : " << input_text_size << std::endl;
    std::cout << "The number of maximum substrings \t : " << ms_count << std::endl;
    std::cout << "Peak count \t : " << st_result.max_nodes_at_level << std::endl;
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
#if DEBUG
    std::cout << "\033[31m";
    std::cout << "DEBUG MODE" << std::endl;
    std::cout << "\033[39m" << std::endl;
#endif

    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file name", true);
    p.add<string>("mode", 'm', "mode", false, "1");
    p.add<int>("thread_num", 'p', "thread number", false, -1);

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    string mode = p.get<string>("mode");

    int thread_num = p.get<int>("thread_num");
    if (thread_num < 0)
    {
        thread_num = std::thread::hardware_concurrency();
    }

    std::cout << "Thread number = " << thread_num << std::endl;

    std::ifstream ifs(inputFile);
    bool inputFileExist = ifs.is_open();
    if (!inputFileExist)
    {
        std::cout << inputFile << " cannot open." << std::endl;
        return -1;
    }

    if(mode != "0"){
    computeMaximalSubstrings(inputFile, mode, thread_num);

    }else{
        computeMaximalSubstrings_beller(inputFile);
    }

}
