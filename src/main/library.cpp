#include <cassert>
#include <chrono>
#include <stdio.h>
#include "../module/stool/src/io.hpp"
#include "../module/stool/src/sa_bwt_lcp.hpp"

#include "../module/stool/src/print.hpp"
#include "../module/stool/src/cmdline.h"
#include "../module/stool/src/debug.hpp"
#include "../module/libdivsufsort/sa.hpp"

//#include "hpp/bwt.hpp"
#include "../basic/interval_search_data_structure.hpp"
//#include "../beller/beller_interval.hpp"
#include "../debug/beller_debug.hpp"

#include "../main/common.hpp"
#include "../debug/naive_algorithms.hpp"
#include "../stnode_enumerator/single/single_stnode_traverser.hpp"
#include "../stnode_enumerator/application.hpp"

#include <sdsl/wt_algorithm.hpp>


template class stool::IntervalSearchDataStructure<uint8_t>;
template class stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<uint32_t>;
template class stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<uint64_t>;
template class stool::stnode_on_rlbwt::SingleSTNodeTraverser<uint32_t, stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<uint32_t>>;
template class stool::stnode_on_rlbwt::SingleSTNodeTraverser<uint64_t, stool::stnode_on_rlbwt::ExplicitWeinerLinkComputer<uint64_t>>;
