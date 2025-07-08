#include <cassert>
#include <chrono>
#include "stool/include/io/io.hpp"
#include "stool/include/strings/sa_bwt_lcp.hpp"

#include "stool/include/debug/print.hpp"
#include "stool/include/third_party/cmdline.h"
#include "stool/include/debug/debug.hpp"
#include "libdivsufsort/sa.hpp"
//#include "../module/rlbwt_iterator/src/include/rlbwt_iterator.hpp"
//#include "module/rlbwt_iterator/src/include/bwt.hpp"
#include "../include/include.hpp"
#include <sdsl/bit_vectors.hpp>



using CHAR = char;
using INDEX = uint64_t;
using LCPINTV = stool::LCPInterval<uint64_t>;

template <typename INDEX_SIZE>
class MUSDetectionArray
{
public:
    std::vector<INDEX_SIZE> counter_array;
    std::vector<INDEX_SIZE> indexes;
    std::vector<stool::LCPInterval<INDEX_SIZE>> mus_intervals;

    MUSDetectionArray()
    {
        uint64_t x = UINT8_MAX + 1;
        this->counter_array.resize(x, 0);
    }
    template <typename EM, typename NODE_ITERATOR>
    void build(EM &em, const NODE_ITERATOR &it, uint64_t lcp)
    {
        //using CHAR = typename EM::CHAR;
        for (auto it : indexes)
        {
            this->counter_array[it] = 0;
        }
        indexes.clear();
        mus_intervals.clear();

        uint64_t w = it.get_children_count();
        for (uint64_t i = 0; i < w; i++)
        {
            uint8_t c = it.get_edge_character(i);
            uint64_t L = it.get_child_left_boundary(i);
            uint64_t R = it.get_child_right_boundary(i);

            uint64_t count = R - L + 1;

            this->counter_array[c] = count;

            if (lcp == 0)
            {
                if (count == 1)
                {
                    this->mus_intervals.push_back(stool::LCPInterval<INDEX_SIZE>(L, R, 1));
                }
            }
        }

        stool::renum::STNodeVector<typename EM::INDEX, typename EM::CHAR> output_vec;

        stool::renum::WeinerLinkCommonFunctions::compute_weiner_links(em, it, output_vec);
        for (auto wnode_it = output_vec.begin(); wnode_it != output_vec.end(); wnode_it++)
        {
            uint64_t w = wnode_it.get_children_count();
            for (uint64_t i = 0; i < w; i++)
            {

                uint8_t c = wnode_it.get_edge_character(i);

                if (this->counter_array[c] >= 2)
                {
                    uint64_t L = wnode_it.get_child_left_boundary(i);
                    uint64_t R = wnode_it.get_child_right_boundary(i);
                    uint64_t count = R - L + 1;
                    if (count == 1)
                    {
                        this->mus_intervals.push_back(stool::LCPInterval<INDEX_SIZE>(L, R, lcp + 2));
                    }
                }
            }
        }
    }
};

class MUSEnumerator
{
public:
    template <typename STNODES>
    static std::vector<stool::LCPInterval<uint64_t>> enumerate(STNODES &stnodeSequencer)
    {
        if (!stnodeSequencer.has_edge_characters())
        {
            std::cout << "Error: This sequence does not have edge characters." << std::endl;
            throw -1;
        }

        std::vector<stool::LCPInterval<uint64_t>> r;
        MUSDetectionArray<uint64_t> mus_arr;
        auto em = stnodeSequencer.get_interval_search_deta_structure();

        auto it = stnodeSequencer.begin();

        while (it != stnodeSequencer.end())
        {
            std::vector<stool::LCPInterval<uint64_t>> r2;
            uint64_t lcp = it.get_depth();
            stool::renum::STDepthIteratorErrorChecker::error_check(it);

            for (auto node_it = it.begin(); node_it != it.end(); node_it++)
            {

                if (node_it.is_maximal_repeat())
                {

                    mus_arr.build(*em, node_it, lcp);
                    for (auto &mus_it : mus_arr.mus_intervals)
                    {
                        r2.push_back(mus_it);
                    }
                }
            }

            for (auto &it : r2)
            {
                r.push_back(it);
            }
            it++;
        }
        std::cout << "STOP" << std::endl;
        return r;
    }
    template <typename STNODES>
    static uint64_t online_enumerate(std::ofstream &out, STNODES &stnodeSequencer, stool::renum::STTreeAnalysisResult &analysis)
    {
        if (!stnodeSequencer.has_edge_characters())
        {
            std::cout << "Error: This sequence does not have edge characters." << std::endl;
            throw -1;
        }

        analysis.start(stnodeSequencer.get_input_text_length());

        using INDEX_SIZE = typename STNODES::index_type;
        uint8_t print_type = 0;
        out.write((char *)(&print_type), sizeof(print_type));
        uint8_t index_bits = sizeof(INDEX_SIZE);
        out.write((char *)(&index_bits), sizeof(index_bits));
        uint64_t count = 0;

        std::vector<stool::LCPInterval<uint64_t>> r;
        MUSDetectionArray<uint64_t> mus_arr;
        auto em = stnodeSequencer.get_interval_search_deta_structure();

        auto it = stnodeSequencer.begin();

        std::vector<stool::LCPInterval<INDEX_SIZE>> buffer;

        while (it != stnodeSequencer.end())
        {

            analysis.analyze(stnodeSequencer);
            std::vector<stool::LCPInterval<uint64_t>> r2;
            uint64_t lcp = it.get_depth();

            for (auto node_it = it.begin(); node_it != it.end(); node_it++)
            {
                if (node_it.is_maximal_repeat())
                {
                    //node_it.error_check();
                    //stool::renum::STDepthIteratorErrorChecker::error_check(node_it);
                    mus_arr.build(*em, node_it, lcp);

                    for (auto &mus_it : mus_arr.mus_intervals)
                    {
                        INDEX_SIZE left = mus_it.i;
                        INDEX_SIZE right = mus_it.j;

                        stool::LCPInterval<INDEX_SIZE> intv(left, right, mus_it.lcp);

                        buffer.push_back(intv);
                        count++;
                        if (buffer.size() >= 8000)
                        {
                            out.write((char *)(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                            buffer.clear();
                        }
                    }
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
        return count;
    }
};
void debug(std::string inputFile)
{
    std::cout << "Constuct Test MUSs" << std::endl;
    std::vector<char> text = stool::bwt::decompress_bwt(inputFile);
    std::vector<uint64_t> sa = stool::construct_suffix_array(text);
    /*
    for (uint64_t p = 0; p < text.size(); p++)
    {
        std::cout << (text[p] == 0 ? "$" : std::string(1, text[p]));
    }
    std::cout << std::endl;
    */

    auto MUSs = stool::esaxx::NaiveAlgorithms::naive_compute_MUSs(text, sa);
    std::sort(MUSs.begin(), MUSs.end(), stool::LCPIntervalDepthOrderComp<uint64_t>());
    std::cout << "Constuct Test MUSs[END]" << std::endl;

    for (auto &it : MUSs)
    {
        std::string s = "";
        for (uint64_t p = 0; p < it.lcp; p++)
        {
            s += text[sa[it.i] + p];
        }
        std::cout << s << std::endl;
    }
    for (auto &it : MUSs)
    {
        std::cout << it.to_string();
    }
    std::cout << std::endl;

    stool::rlbwt2::BWTAnalysisResult analysisResult;
    stool::renum::RLE<uint8_t> rlbwt;
    rlbwt.load(inputFile, analysisResult);
    using RDS = stool::renum::RLEWaveletTree<uint32_t>;
    RDS ds = RDS(&rlbwt);

    std::cout << "Enumerate Minimal unique substrings..." << std::endl;
    stool::renum::SuffixTreeNodes<uint32_t, RDS> stnodeTraverser;
    stnodeTraverser.initialize(1, ds, true);
    auto test_MUSs = MUSEnumerator::enumerate(stnodeTraverser);

    std::sort(test_MUSs.begin(), test_MUSs.end(), stool::LCPIntervalDepthOrderComp<uint64_t>());
    /*
    for (auto &it : test_MUSs)
    {
        std::cout << it.to_string();
    }
    std::cout << std::endl;
    */
    stool::renum::equal_check_lcp_intervals(MUSs, test_MUSs);
    std::cout << "OK!" << MUSs.size() << std::endl;
}

void computeMUSs(std::string inputFile, std::string outputFile, int thread_num)
{

    auto start = std::chrono::system_clock::now();
    std::string bit_size_mode = "UINT64_t";
    std::chrono::system_clock::time_point mid;

    std::ofstream out(outputFile, std::ios::out | std::ios::binary);
    if (!out)
    {
        throw std::runtime_error("Cannot open the output file!");
    }
    stool::rlbwt2::BWTAnalysisResult analysisResult;
    stool::renum::RLE<uint8_t> rlbwt;
    rlbwt.load(inputFile, analysisResult);

    uint64_t ms_count = 0;
    stool::renum::STTreeAnalysisResult st_result;
    uint64_t ds_memory_usage = 0;

    if (analysisResult.str_size < UINT32_MAX - 10)
    {
        using RDS = stool::renum::RLEWaveletTree<uint32_t>;
        RDS ds = RDS(&rlbwt);
        mid = std::chrono::system_clock::now();
        ds_memory_usage = ds.get_using_memory();

        std::cout << "Enumerate Minimal unique substrings..." << std::endl;
        stool::renum::SuffixTreeNodes<uint32_t, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds, true);
        //ms_count = stool::renum::Application::outputMaximalSubstrings(out, stnodeTraverser, st_result);
        ms_count = MUSEnumerator::online_enumerate(out, stnodeTraverser, st_result);
        bit_size_mode = "UINT32_t";
    }
    else
    {
        using RDS = stool::renum::RLEWaveletTree<uint64_t>;
        RDS ds = RDS(&rlbwt);
        mid = std::chrono::system_clock::now();
        ds_memory_usage = ds.get_using_memory();

        std::cout << "Enumerate Minimal unique substrings..." << std::endl;
        stool::renum::SuffixTreeNodes<uint64_t, RDS> stnodeTraverser;
        stnodeTraverser.initialize(thread_num, ds, true);
        ms_count = MUSEnumerator::online_enumerate(out, stnodeTraverser, st_result);
    }

    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double enumeration_time = std::chrono::duration_cast<std::chrono::milliseconds>(end - mid).count();
    double construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(mid - start).count();

    double bps = ((double)analysisResult.str_size / ((double)elapsed / 1000)) / 1000;

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "RLBWT File \t\t\t\t : " << inputFile << std::endl;
    std::cout << "Output \t\t\t\t\t : " << outputFile << std::endl;
    std::cout << "Peak children count \t\t\t : " << st_result.max_nodes_at_level << std::endl;
    //std::cout << "Data structure Size \t\t\t : " << (data_structure_bytes / 1000) << "[KB]" << std::endl;

    std::cout << "Thread number \t\t\t\t : " << thread_num << std::endl;
    std::cout << "Integer Type \t\t\t\t : " << bit_size_mode << std::endl;

    std::cout << "The length of the input text \t\t : " << analysisResult.str_size << std::endl;
    std::cout << "The number of runs on BWT \t\t : " << analysisResult.run_count << std::endl;
    std::cout << "The number of MUSs \t : " << ms_count << std::endl;

    std::cout << "\033[33m";
    std::cout << "______________________Execution Time______________________" << std::endl;
    std::cout << "Excecution time \t\t\t : " << elapsed << " [s]" << std::endl;
    std::cout << "Character per second \t\t\t : " << bps << " [KB/s]" << std::endl;
    std::cout << "\t Preprocessing time \t\t : " << construction_time << " [s]" << std::endl;
    std::cout << "\t Enumeration time \t\t : " << enumeration_time << " [s]" << std::endl;

    std::cout << "\033[32m";
    std::cout << "______________________Memory Usage______________________" << std::endl;
    std::cout << "RLBWT \t : " << (rlbwt.get_using_memory() / 1000) << " [KB]" << std::endl;
    std::cout << "WT \t : " << (ds_memory_usage / 1000) << " [KB]" << std::endl;
    std::cout << "Queue \t : " << st_result.peak_memory_of_queue / 1000 << " [KB]" << std::endl;
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
    p.add<std::string>("input_file", 'i', "input file name", true);
    p.add<std::string>("output_file", 'o', "output file path (default: input_file_path.mus)", false, "");
    p.add<int>("thread_num", 'p', "thread number", false, -1);

    p.parse_check(argc, argv);
    std::string inputFile = p.get<std::string>("input_file");
    //string mode = p.get<string>("mode");

    std::string outputFile = p.get<std::string>("output_file");
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

    if (outputFile.size() == 0)
    {
        outputFile = inputFile + ".mus";
    }
    //debug(inputFile);
    computeMUSs(inputFile, outputFile, thread_num);
}
