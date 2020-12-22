#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include <queue>
#include "../test/naive_algorithms.hpp"
#include "../rlbwt/bwt_decompress.hpp"
#include "../rlbwt/range_distinct/light_range_distinct.hpp"

//#include "../weiner_interval.hpp"
#include "../rlbwt/rinterval.hpp"
#include "../rlbwt/bwt_decompress.hpp"
#include "../beller/fmindex.hpp"
#include "../beller/beller_debug.hpp"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wt_gmr.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>

#include <mutex>
namespace stool
{
    namespace lcp_on_rlbwt
    {
        template <typename CHAR_VEC, typename INDEX_SIZE>
        class LightRangeDistinctDataStructure2
        {
        private:
            std::mutex mtx;

        public:
            const CHAR_VEC *_char_vec;

            void preprocess(CHAR_VEC *__char_vec)
            {
                this->_char_vec = __char_vec;
            }
            uint64_t range_distinct(INDEX_SIZE i, INDEX_SIZE j, std::vector<CharInterval<INDEX_SIZE>> &output)
            {
                int32_t charMaxSize = ((int32_t)UINT8_MAX) + 1;

                std::vector<int64_t> checker;
                std::vector<uint64_t> beginIndexes;
                std::vector<uint64_t> endIndexes;
                std::vector<uint64_t> charIndexes;
                checker.resize(charMaxSize, -1);
                beginIndexes.resize(charMaxSize, 0);
                endIndexes.resize(charMaxSize, 0);
                charIndexes.resize(charMaxSize, 0);

                std::lock_guard<std::mutex> lock(std::mutex);
                uint64_t count = 0;
                assert(i <= j);
                for (uint64_t x = i; x <= j; x++)
                {
                    uint8_t c = (*this->_char_vec)[x];
                    if (checker[c] == -1)
                    {
                        checker[c] = count;
                        beginIndexes[count] = x;
                        endIndexes[count] = x;
                        charIndexes[count] = c;
                        count++;
                    }
                    else
                    {
                        endIndexes[checker[c]] = x;
                    }
                }
                for (uint64_t x = 0; x < count; x++)
                {
                    output[x] = CharInterval<INDEX_SIZE>(beginIndexes[x], endIndexes[x], charIndexes[x]);
                    checker[charIndexes[x]] = -1;
                }

                assert(count > 0);
                return count;
            }
        };

        class STNodeChecker
        {
        private:
            std::mutex mtx;

        public:
            std::vector<stool::LCPInterval<uint64_t>> lcp_intervals;
            std::vector<std::unordered_map<uint64_t, uint64_t>> maps;
            //std::vector<std::unordered_set<uint64_t>> sets;
            stool::lcp_on_rlbwt::LightRangeDistinctDataStructure2<sdsl::int_vector<>, uint64_t> lrdds;
            sdsl::int_vector<> bwt;
            wt_gmr<> wt;
            std::vector<uint64_t> C;

            uint64_t current_lcp = 0;

            void initialize(string inputFile)
            {
                std::vector<uint8_t> _bwt;
                stool::bwt::load(inputFile, _bwt);
                this->bwt.resize(_bwt.size());
                for (uint64_t i = 0; i < _bwt.size(); i++)
                {
                    this->bwt[i] = _bwt[i];
                }
                stool::FMIndex::constructSelect(this->bwt, wt);
                stool::FMIndex::constructC(this->bwt, this->C);
                this->lrdds.preprocess(&this->bwt);

                std::cout << "Construct Checker" << std::endl;

                std::vector<char> text = stool::bwt::decompress_bwt(inputFile);
                vector<uint64_t> sa = stool::construct_suffix_array(text);

                vector<stool::LCPInterval<uint64_t>> correct_lcp_intervals = stool::esaxx::naive_compute_lcp_intervals<char, uint64_t>(text, sa);
                std::sort(
                    correct_lcp_intervals.begin(),
                    correct_lcp_intervals.end(),
                    stool::LCPIntervalPreorderComp<uint64_t>());

                for (auto &it : correct_lcp_intervals)
                {
                    if (it.lcp < bwt.size())
                    {
                        this->lcp_intervals.push_back(it);
                        if (text.size() < 100)
                        {
                            std::cout << it.to_string();
                        }
                    }
                    else if (bwt[it.i] == 0 || bwt[it.i] == '$')
                    {
                        std::cout << "Doller " << bwt[it.i] << std::endl;
                        it.lcp = bwt.size();
                        this->lcp_intervals.push_back(it);
                    }
                }
                std::cout << std::endl;

                //this->lcp_intervals.swap(correct_lcp_intervals);

                this->maps.resize(text.size());

                for (auto &it : this->lcp_intervals)
                {
                    if (it.lcp < text.size())
                    {
                        this->maps[it.lcp][it.i] = it.j;
                    }
                }
                std::cout << "Construct Checker[Finished]" << std::endl;
            }
            void increment(uint64_t k)
            {

                if (this->maps[this->current_lcp].size() != k)
                {

                    std::cout << "INCREMENT ERROR!" << " LCP = " << this->current_lcp << ", Collect : " << this->maps[this->current_lcp].size() << "/Test: " << k << std::endl;
                    assert(this->maps[this->current_lcp].size() == k);

                    throw -1;
                }

                this->current_lcp++;
            }
            bool check2(uint64_t i, uint64_t j, std::vector<stool::LCPInterval<uint64_t>> &wlinks)
            {
/*
#if DEBUG
                std::cout << "Weiner Link Check [" << i << ", " << j << "]" << std::endl;
#endif
*/

                std::lock_guard<std::mutex> lock(std::mutex);
                std::vector<stool::LCPInterval<uint64_t>> correctWlinks;

                if (this->current_lcp == 0)
                {
                    correctWlinks.push_back(stool::LCPInterval<uint64_t>(0, bwt.size() - 1, this->current_lcp));
                }
                else
                {
                    std::vector<CharInterval<uint64_t>> output;
                    output.resize(256);

                    uint64_t count = this->lrdds.range_distinct(i, j, output);

                    for (uint64_t x = 0; x < count; x++)
                    {
                        uint64_t l = stool::FMIndex::LF(output[x].i, bwt, C, wt);
                        uint64_t r = stool::FMIndex::LF(output[x].j, bwt, C, wt);

                        auto it = this->maps[this->current_lcp].find(l);
                        if (it != this->maps[this->current_lcp].end())
                        {
/*
#if DEBUG
                        std::cout << "Next [" << l << ", " << r << "]" << std::endl;
#endif
*/

                            correctWlinks.push_back(stool::LCPInterval<uint64_t>(l, r, this->current_lcp));
                        }
                    }
                }

                std::sort(
                    wlinks.begin(),
                    wlinks.end(),
                    stool::LCPIntervalPreorderComp<uint64_t>());

                std::sort(
                    correctWlinks.begin(),
                    correctWlinks.end(),
                    stool::LCPIntervalPreorderComp<uint64_t>());

                stool::beller::equal_check_lcp_intervals(wlinks, correctWlinks, "Weiner link check");
                return true;
            }

            bool check_lcp_interval(uint64_t i, uint64_t j)
            {
                std::lock_guard<std::mutex> lock(std::mutex);

                auto it = this->maps[this->current_lcp].find(i);
                if (it == this->maps[this->current_lcp].end())
                {
                    std::cout << "LCP Interval CHECK ERROR!(1)" << std::endl;
                    std::cout << "Test LCP Interval = "
                              << "[" << i << ", " << j << "]" << std::endl;

                    throw -1;
                }
                else
                {
                    if (it->second == j)
                    {
                        return true;
                    }
                    else
                    {
                        std::cout << "LCP Interval CHECK ERROR!(2)" << std::endl;
                        std::cout << "Test LCP Interval = "
                                  << "[" << i << ", " << j << "]" << std::endl;
                        std::cout << "Collect LCP Interval = "
                                  << "[" << i << ", " << it->second << "]" << std::endl;

                        throw -1;
                    }
                }
            }
        };
    } // namespace lcp_on_rlbwt
} // namespace stool