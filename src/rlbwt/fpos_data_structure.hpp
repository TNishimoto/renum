#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>
#include "stool/src/elias_fano_vector.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        template <typename T>
        class RankSupportVectorWrapper
        {
            T &items;

        public:
            RankSupportVectorWrapper(T &_items) : items(_items)
            {
            }
            uint64_t operator[](uint64_t index) const
            {
                return this->items[index];
            }
            uint64_t rank(uint64_t value) const
            {
                auto p2 = std::lower_bound(this->items.begin(), this->items.end(), value);
                uint64_t pos2 = std::distance(this->items.begin(), p2);

                return pos2;
            }
            uint64_t size() const
            {
                return this->items.size();
            }
        };
        class FPosDataStructure
        {
        public:
            //std::vector<uint64_t> fposArray;
            /*
    Let R be the output array of length r.
    R[i] stores the index of the L-run corresponding to i-th F-run.
    */
            static std::vector<uint64_t> construct_rle_fl_mapper(const sdsl::int_vector<> &bwt)
            {
                std::vector<uint64_t> indexes;
                uint64_t size = bwt.size();
                indexes.resize(size);
                for (uint64_t i = 0; i < size; i++)
                {
                    indexes[i] = i;
                }
                sort(indexes.begin(), indexes.end(),
                     [&](const uint64_t &x, const uint64_t &y) {
                         assert(x < bwt.size() && y < bwt.size());
                         if (bwt[x] == bwt[y])
                         {
                             return x < y;
                         }
                         else
                         {
                             return (uint64_t)bwt[x] < (uint64_t)bwt[y];
                         }
                     });
                return indexes;
            }
            /*
    Let R be the output array of length r.
    R[i] stores the starting position of the F-run corresponding to i-th L-run.
    */
            template <typename LPOSDS>
            static std::vector<uint64_t> construct_fpos_array(const sdsl::int_vector<> &bwt, const LPOSDS &lpos_vec)
            {
                std::vector<uint64_t> fvec = construct_rle_fl_mapper(bwt);
                std::vector<uint64_t> output;
                output.resize(fvec.size(), 0);
                uint64_t fsum = 0;
                for (uint64_t i = 0; i < fvec.size(); i++)
                {
                    uint64_t x = fvec[i];
                    assert(x < output.size());
                    output[x] = fsum;
                    uint64_t run = lpos_vec[x + 1] - lpos_vec[x];
                    /*
                    if(bwt[x] == 97){
                    std::cout << "[" << fsum << ", " << (fsum + run -1) << "], " << std::flush;

                    }
                    */
                    fsum += run;
                }
                return output;
            }
            template <typename LPOSDS>
            static void check(const sdsl::int_vector<> &bwt, const LPOSDS &lpos_vec, const std::vector<uint64_t> &r)
            {
                std::vector<bool> checker;

                uint64_t size = lpos_vec[lpos_vec.size() - 1];
                checker.resize(size, false);
                for (uint64_t i = 0; i < 256; i++)
                {
                    for (uint64_t x = 0; x < bwt.size(); x++)
                    {

                        if (bwt[x] == i)
                        {
                            uint64_t run = lpos_vec[x + 1] - lpos_vec[x];
                            for (uint64_t p = 0; p < run; p++)
                            {
                                uint64_t y = r[x] + p;
                                assert(y < checker.size());

                                if (y > 0 && !checker[y - 1])
                                {
                                    //std::cout << i << "/" << x << std::endl;
                                    assert(false);
                                }
                                checker[y] = true;
                            }
                        }
                    }
                }

                std::cout << "OK!" << std::endl;
            }
            template <typename LPOSDS>
            static std::vector<uint64_t> construct(const sdsl::int_vector<> &bwt, const LPOSDS &lpos_vec)
            {
                std::cout << "Building Fpos_vec(std::vector<uint64_t>) ..." << std::flush;
                std::vector<uint64_t> v1 = construct_fpos_array(bwt, lpos_vec);
                std::cout << "[Finished]" << std::endl;

                //check(bwt, lpos_vec, v1);
                return v1;
            }
            /*
            uint64_t &operator[](uint64_t i)
            {
                return fposArray[i];
            }
            */
        };
        class LightFPosDataStructure
        {
        public:
            const sdsl::int_vector<> *bwt;
            stool::WT *wt;

            std::vector<stool::EliasFanoVector> efv_vec;
            std::vector<uint64_t> C2;

            LightFPosDataStructure()
            {

            }
            template <typename LPOSVEC>
            void build(const sdsl::int_vector<> *_bwt, const LPOSVEC &_lposvec, stool::WT *_wt){
                this->bwt = _bwt;
                this->wt = _wt;
                std::cout << "Building Fpos_vec(EliasFanoVector) ..." << std::flush;
                /*
                LightFPosDataStructure::construct_C(bwt, this->C);
                LightFPosDataStructure::construct_sorted_fpos_array(_bwt, _lposvec, this->fposSortedArray);
                */
#if DEBUG
                this->rank_test();
#endif
                this->build(_lposvec);
                std::cout << "[Finished]" << std::endl;

#if DEBUG
                this->check(_lposvec);
#endif

            }
            template <typename LPOSVEC>
            void check(const LPOSVEC &_lposvec)
            {

                std::vector<uint64_t> collectVec = FPosDataStructure::construct(*this->bwt, _lposvec);
                /*
                std::vector<uint64_t> C;
                stool::EliasFanoVector fposSortedArray;
                LightFPosDataStructure::construct_C(bwt, C);
                LightFPosDataStructure::construct_sorted_fpos_array(this->bwt, _lposvec, fposSortedArray);
                */

                for (uint64_t i = 0; i < this->bwt->size(); i++)
                {
                    uint8_t c = (*this->bwt)[i];
                    uint64_t rank1 = wt->rank(i + 1, c);
                    assert(rank1 < this->efv_vec[c].size());
                    uint64_t q = C2[c] + this->efv_vec[c][rank1];

                    uint64_t p = collectVec[i];

                    if (p != q)
                    {
                        std::cout << "LightFPosDataStructure Error, i = " << i << ", Collect = " << p << ", Test = " << q << std::endl;
                        throw -1;
                    }
                }
            }

            static void construct_C(const sdsl::int_vector<> &bwt_head_chars, std::vector<uint64_t> &C)
            {
                uint64_t CHARMAX = UINT8_MAX + 1;
                std::vector<uint64_t> CK;
                CK.resize(CHARMAX, 0);
                C.resize(CHARMAX, 0);
                uint64_t rle = bwt_head_chars.size();
                for (uint64_t i = 0; i < rle; i++)
                {
                    uint8_t c = bwt_head_chars[i];
                    //uint64_t l = rlbwt.get_run(i);
                    CK[c] += 1;
                }
                for (uint64_t i = 1; i < CHARMAX; i++)
                {
                    C[i] = C[i - 1] + CK[i - 1];
                }
            }
#if DEBUG
            void rank_test()
            {
                uint64_t CHARMAX = UINT8_MAX + 1;
                std::vector<uint64_t> C_run_sum;
                C_run_sum.resize(CHARMAX, 0);
                uint64_t rle = this->bwt->size();
                for (uint64_t i = 0; i < rle; i++)
                {
                    uint8_t c = (*this->bwt)[i];
                    uint64_t rank1 = wt->rank(i + 1, c);
                    assert(C_run_sum[c] == rank1);
                    C_run_sum[c]++;
                }
            }
#endif
            template <typename LPOSVEC>
            void build(const LPOSVEC &lposvec)
            {
                uint64_t CHARMAX = UINT8_MAX + 1;
                std::vector<uint64_t> C_run_sum;
                std::vector<uint64_t> numVec;
                std::vector<uint64_t> numC;

                std::vector<stool::EliasFanoVectorBuilder> builders;

                C_run_sum.resize(CHARMAX, 0);
                numVec.resize(CHARMAX, 0);
                this->C2.resize(CHARMAX, 0);
                this->efv_vec.resize(CHARMAX);
                builders.resize(CHARMAX);
                uint64_t rle = this->bwt->size();



                for (uint64_t i = 0; i < rle; i++)
                {

                    uint8_t c = (*this->bwt)[i];
                    uint64_t l = lposvec[i + 1] - lposvec[i];
                    C_run_sum[c] += l;
                    numVec[c]++;
                }

                for (uint64_t i = 1; i < CHARMAX; i++)
                {
                    this->C2[i] = this->C2[i - 1] + C_run_sum[i - 1];
                }

                for (uint64_t i = 0; i < CHARMAX; i++)
                {
                    builders[i].initialize(C_run_sum[i] + 1, numVec[i]);
                }

                std::vector<uint64_t> tmp_sum;
                tmp_sum.resize(CHARMAX, 0);

                for (uint64_t i = 0; i < rle; i++)
                {
                    uint8_t c = (*this->bwt)[i];
                    //std::cout << (uint64_t)c << "/" << numVec[c]<< std::endl;
                    builders[c].push(tmp_sum[c]);
                    //std::cout << "a" << std::endl;
                    uint64_t l = lposvec[i + 1] - lposvec[i];
                    tmp_sum[c] += l;
                }

                for (uint64_t i = 0; i < CHARMAX; i++)
                {
                    builders[i].finish();
                    this->efv_vec[i].build_from_builder(builders[i]);
                }
            }
            /*
            template <typename LPOSVEC>
            static void construct_sorted_fpos_array(const sdsl::int_vector<> &bwt_head_chars, const LPOSVEC &lposvec, stool::EliasFanoVector &output)
            {
                uint64_t CHARMAX = UINT8_MAX + 1;
                std::vector<uint64_t> C_run_sum;
                std::vector<uint64_t> C;
                std::vector<uint64_t> numVec;
                std::vector<uint64_t> numC;

                std::vector<uint64_t> CK3;

                C_run_sum.resize(CHARMAX, 0);
                numVec.resize(CHARMAX, 0);
                numC.resize(CHARMAX, 0);
                C.resize(CHARMAX, 0);
                CK3.resize(CHARMAX, 0);

                uint64_t rle = bwt_head_chars.size();
                for (uint64_t i = 0; i < rle; i++)
                {
                    uint8_t c = bwt_head_chars[i];
                    uint64_t l = lposvec[i + 1] - lposvec[i];
                    C_run_sum[c] += l;
                    numVec[c]++;
                }


                for (uint64_t i = 1; i < CHARMAX; i++)
                {
                    C[i] = C[i - 1] + C_run_sum[i - 1];
                    numC[i] = numC[i - 1] + numVec[i - 1];
                }

                std::vector<uint64_t> X;
                X.resize(rle);
                for (uint64_t i = 0; i < rle; i++)
                {
                    uint8_t c = bwt_head_chars[i];
                    uint64_t l = lposvec[i + 1] - lposvec[i];
                    X[numC[c]] = C[c];
                    numC[c]++;
                    C[c] += l;
                }
                output.construct(&X);
            }
            */

            uint64_t operator[](uint64_t i) const
            {

                uint8_t c = (*this->bwt)[i];
                uint64_t rank1 = (*wt).rank(i + 1, c);

                assert(rank1 < this->efv_vec[c].size());
                uint64_t q = C2[c] + this->efv_vec[c][rank1];
                return q;
            }
            uint64_t size() const
            {
                return this->bwt->size();
            }
            uint64_t get_using_memory() const
            {
                uint64_t x = 0;

                for (auto &it : this->efv_vec)
                {
                    x += it.get_using_memory();
                }
                return x;
            }
        };
    } // namespace lcp_on_rlbwt
} // namespace stool