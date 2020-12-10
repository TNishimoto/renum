
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "stool/src/cmdline.h"
#include "../../module/rlbwt_iterator/src/include/rlbwt_iterator.hpp"

using namespace std;
using CHAR = uint8_t;
using INDEX = uint64_t;
#include "../test/old_postorder_maximal_substrings.hpp"
#include "../main/common.hpp"
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp_dac.hpp>
#include <sdsl/lcp_support_sada.hpp>
#include "../postorder_suffix_tree_intervals.hpp"
#include "../test/naive_algorithms.hpp"
#include "../postorder_maximal_substring_intervals.hpp"
#include "../src/minimal_substrings/minimal_substring_iterator.hpp"
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp_dac.hpp>
#include <sdsl/lcp_support_sada.hpp>
#include "../forward_bwt.hpp"




void test_array(vector<uint64_t> &sa, vector<uint64_t> &lcp, vector<uint64_t> &test_sa, vector<uint64_t> &test_lcp, std::string name)
{
    std::cout << "Test:" << name << std::flush;
    stool::equal_check(sa, test_sa);
    std::cout << "[OK!]" << std::endl;
    stool::equal_check(lcp, test_lcp);
    std::cout << "[OK!]" << std::endl;
    //super_equal_check(bwt, test_bwt);
    //std::cout << "[OK!]" << std::endl;
}

int main(int argc, char *argv[])
{

    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file name", true);

    p.parse_check(argc, argv);
    string filename = p.get<string>("input_file");

    vector<char> text = stool::load_char_vec_from_file(filename, true); // input text
    bool is_contained_minus_character = stool::esaxx::check_test(text);
    if (is_contained_minus_character)
    {
        std::cout << "This text contains minus character!" << std::endl;
    }

    using SA_PLAIN = std::vector<INDEX>;
    using SA_SDSL = sdsl::csa_sada<>;
    using LCP_PLAIN = std::vector<INDEX>;
    using LCPINTERVALS = vector<stool::LCPInterval<INDEX>>;
    using LCP_SDSL = sdsl::lcp_dac<>;
    using LCP_RLBWT = stool::rlbwt::ForwardLCPArray<std::vector<INDEX>>;
    using SA_RLBWT = stool::rlbwt::ForwardSA<std::vector<INDEX>>;
    using BWT_RLBWT = stool::rlbwt::ForwardBWT<std::vector<char>, std::vector<INDEX>>;

    stool::rlbwt::RLBWT<std::vector<char>, std::vector<INDEX>> rlestr;
    stool::rlbwt::Constructor::construct_from_file<>(rlestr, filename);
    LCP_RLBWT lcp_rlbwt;
    lcp_rlbwt.construct_from_rlbwt(&rlestr, false);

    SA_PLAIN sa_naive = stool::esaxx::construct_naive_SA_with_uint64<char, INDEX>(text);
    LCP_PLAIN lcp_naive = stool::constructLCP<char, INDEX>(text, sa_naive);

    SA_SDSL sa_sdsl;
    construct(sa_sdsl, filename, 1);
    std::vector<INDEX> sa_sdsl_vec, lcp_sdsl_vec;
    for (auto it : sa_sdsl)
        sa_sdsl_vec.push_back(it);

    LCP_SDSL lcp_sdsl;
    construct(lcp_sdsl, filename, 1);
    for (auto it : lcp_sdsl)
        lcp_sdsl_vec.push_back(it);

    test_array(sa_naive, lcp_naive, sa_sdsl_vec, lcp_sdsl_vec, "SDSL");

    SA_RLBWT *sa_pointer = const_cast<SA_RLBWT *>(lcp_rlbwt.get_ForwardSA());
    std::vector<uint64_t> _sa_rlbwt;
    for (auto it : *sa_pointer)
        _sa_rlbwt.push_back(it);

    std::vector<uint64_t> _lcp_rlbwt;
    for (auto it : lcp_rlbwt)
        _lcp_rlbwt.push_back(it);

    test_array(sa_naive, lcp_naive, _sa_rlbwt, _lcp_rlbwt, "RLBWT");

}