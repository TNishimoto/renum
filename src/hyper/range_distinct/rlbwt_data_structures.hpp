#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>
//#include "./rinterval_storage.hpp"
#include "./range_distinct_on_rlbwt.hpp"
#include "../../test/stnode_checker.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {

        template <typename INDEX_SIZE, typename LPOSDS, typename FPOSDS>
        class RLBWTDataStructures
        {

        public:
                       stool::lcp_on_rlbwt::STNodeChecker *stnc;

            /*
            using CHARVEC = typename RLBWT_STR::char_vec_type;
            */

            using INDEX = INDEX_SIZE;
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using CHAR = uint8_t;
            using UCHAR = typename std::make_unsigned<CHAR>::type;
            const sdsl::int_vector<> &bwt;
            using CHAR_VEC = sdsl::int_vector<>;

            stool::WT &wt;
            const LPOSDS &lpos_vec;
            const FPOSDS &_fposDS;
            //RangeDistinctDataStructureOnRLBWT<INDEX_SIZE> rangeOnRLBWT;

            //const RLBWT_STR &_rlbwt;


            std::vector<stool::LCPInterval<uint64_t>> *collectLCPIntervals;

            bool checkLCPInterval(const RINTERVAL &input)
            {
                if(this->stnc == nullptr){
                    std::cout << "stnc is null" << std::endl;
                    throw -1;
                }
                stool::LCPInterval<uint64_t> intv2 = input.get_lcp_interval(this->stnc->current_lcp, this->_fposDS);
                return this->stnc->check(intv2.i, intv2.j);
            }
            bool checkWeinerLink(const RINTERVAL &input, std::vector<RINTERVAL> &stnodeVec, std::vector<uint64_t> &indexVec, uint64_t indexCount)
            {
                if(this->stnc == nullptr){
                    std::cout << "stnc is null" << std::endl;
                    throw -1;
                }

                stool::LCPInterval<uint64_t> intv2 = input.get_lcp_interval(this->stnc->current_lcp-1, this->_fposDS);

                //std::cout << "PREV = [" <<    intv2.i << ", " << intv2.j << "]" << std::endl;
                //std::cout << "NEXT: " ;                     
                std::vector<stool::LCPInterval<uint64_t>> wlinks;
                for(uint64_t i=0;i<indexCount;i++){
                    stool::LCPInterval<uint64_t> intv = stnodeVec[indexVec[i]].get_lcp_interval(this->stnc->current_lcp, this->_fposDS);
                    //std::cout << "[" << intv.i << ", " << intv.j << "]";
                    wlinks.push_back(intv);
                }
                //std::cout << std::endl;

                return this->stnc->check2(intv2.i, intv2.j, wlinks);
            }


            RLBWTDataStructures(const sdsl::int_vector<> &diff_char_vec,
                                stool::WT &_wt, const LPOSDS &_lpos_vec, const FPOSDS &__fposDS) : bwt(diff_char_vec), wt(_wt), lpos_vec(_lpos_vec), _fposDS(__fposDS)
            {

                //rangeOnRLBWT.initialize(&wt, &bwt);
            }

            INDEX_SIZE get_fpos(INDEX_SIZE index, INDEX_SIZE diff) const
            {
                return _fposDS[index] + diff;
            }
            uint64_t get_lindex_containing_the_position(uint64_t lposition) const
            {
                uint64_t x = this->lpos_vec.rank(lposition + 1) - 1;
                return x;
            }
            uint64_t str_size() const
            {
                return lpos_vec[lpos_vec.size() - 1];
            }
            uint8_t get_char_by_run_index(uint64_t _run_index) const
            {
                return bwt[_run_index];
            }
            uint64_t rle_size() const
            {
                return bwt.size();
            }
            uint64_t get_run(uint64_t i) const
            {
                return lpos_vec[(i + 1)] - lpos_vec[i];
            }
            uint64_t get_lpos(uint64_t i) const
            {
                return lpos_vec[i];
            }

            uint64_t get_end_rle_lposition() const
            {
                //uint64_t min_char = UINT64_MAX;
                //uint64_t min_pos = UINT64_MAX;

                for (INDEX i = 0; i < bwt.size(); i++)
                {
                    if (bwt[i] == 0)
                    {
                        //min_pos = i;
                        return i;
                    }
                }
                std::cout << "No doller character!" << std::endl;
                throw -1;
                //return min_pos;
                //return std::numeric_limits<INDEX>::max();
            }

            RINTERVAL getIntervalOnL(const RINTERVAL &interval) const
            {
                RINTERVAL output;
                INDEX_SIZE begin_pos = this->get_fpos(interval.beginIndex, interval.beginDiff);
                output.beginIndex = this->get_lindex_containing_the_position(begin_pos);
                output.beginDiff = begin_pos - lpos_vec[output.beginIndex];

                INDEX_SIZE end_pos = this->get_fpos(interval.endIndex, interval.endDiff);

                output.endIndex = this->get_lindex_containing_the_position(end_pos);
                output.endDiff = end_pos - lpos_vec[output.endIndex];
                return output;
            }
        };
    } // namespace lcp_on_rlbwt
} // namespace stool