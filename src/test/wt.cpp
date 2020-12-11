#include <cassert>
#include <chrono>

#include "stool/src/io.hpp"
#include "stool/src/cmdline.h"
#include "stool/src/debug.hpp"
#include <sdsl/bit_vectors.hpp>
#include <sdsl/wt_algorithm.hpp>
#include <sdsl/wavelet_trees.hpp>
#include "../beller/char_interval.hpp"

using namespace std;
using namespace stool;



int main(int argc, char *argv[])
{
    cmdline::parser p;
    p.add<uint64_t>("len", 'i', "length", true);
    p.add<uint64_t>("width", 'w', "width", true);

    p.parse_check(argc, argv);
    uint64_t len = p.get<uint64_t>("len");
    uint64_t w = p.get<uint64_t>("width");
    
    sdsl::int_vector<0> iv;
    iv.width(w);
    iv.resize(len);
    sdsl::util::set_random_bits(iv);
    uint64_t x = 0;
    for(uint64_t i=0;i<iv.size();i++){
        x += iv[i];
    }
    std::cout << x << std::endl;

    stool::WT wt;
    construct_im(wt, iv);

}
