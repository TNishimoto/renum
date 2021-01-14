#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include "./bwt_analysis_result.hpp"
#include "stool/src/elias_fano_vector.hpp"

namespace stool
{
    namespace stnode_on_rlbwt
    {

        template <typename CHAR = uint8_t>
        class RLE
        {
        public:
            using char_type = CHAR;
            using INDEX = uint64_t;

        private:
            sdsl::int_vector<> head_char_vec;
            stool::EliasFanoVector lpos_vec;
            INDEX smallest_character = 0;

        public:
            RLE()
            {
            }

            const sdsl::int_vector<> *get_head_char_vec() const
            {
                return &this->head_char_vec;
            }
            const stool::EliasFanoVector *get_lpos_vec() const
            {
                return &this->lpos_vec;
            }
            CHAR get_char_by_run_index(INDEX _run_index) const
            {
                return (this->head_char_vec)[_run_index];
            }
            /*
            void set(CHARVEC &&__head_char_vec, POWVEC &&__lpos_vec)
            {
                this->head_char_vec = new CHARVEC(std::move(__head_char_vec));
                this->lpos_vec = new POWVEC(std::move(__lpos_vec));
                deleteFlag = true;
            }
            */
            /*
            POWVEC *get_lpos_vec_no_const() const
            {
                return this->lpos_vec;
            }

            

            void clear()
            {
                std::cout << "clear!" << std::endl;
                this->head_char_vec->resize(0);
                this->head_char_vec->shrink_to_fit();
                this->lpos_vec->resize(0);
                this->lpos_vec->shrink_to_fit();
            }
            */
            /*
            std::vector<INDEX> construct_pows() const
            {
                std::vector<INDEX> r;
                for (INDEX i = 0; i < lpos_vec->size(); i++)
                {
                    r.push_back((lpos_vec)[i + 1] - (lpos_vec)[i]);
                }
                return r;
            }
            void write(std::ofstream &out) const
            {
                out.write((const char *)(&(*this->head_char_vec)[0]), sizeof(CHAR) * this->head_char_vec->size());
                std::vector<INDEX> pows = this->construct_pows();
                out.write((const char *)(&pows[0]), sizeof(INDEX) * pows.size());
            }
            void write(std::string filename) const
            {
                std::ofstream out(filename, std::ios::out | std::ios::binary);
                this->write(out);
                out.close();
            }
            */
            void load(std::string filename, stool::rlbwt2::BWTAnalysisResult &analysisResult)
            {

                this->head_char_vec.width(8);
                stool::EliasFanoVectorBuilder run_bits;
                std::ifstream inp;
                std::vector<char> buffer;
                uint64_t bufferSize = 8192;
                buffer.resize(8192);

                //stool::rlbwt2::BWTAnalysisResult analysisResult;
                analysisResult.analyze(filename);
                this->smallest_character = analysisResult.min_char;

                this->head_char_vec.resize(analysisResult.run_count);
                run_bits.initialize(analysisResult.str_size + 1, analysisResult.run_count + 1);

                inp.open(filename, std::ios::binary);
                bool inputFileExist = inp.is_open();
                if (!inputFileExist)
                {
                    std::cout << filename << " cannot open." << std::endl;

                    throw std::runtime_error("error");
                }
                uint64_t textSize = stool::FileReader::getTextSize(inp);
                uint8_t prevChar = 255;
                uint64_t x = 0;
                uint64_t currentRunP = 0;

                while (true)
                {
                    bool b = stool::FileReader::read(inp, buffer, bufferSize, textSize);
                    if (!b)
                    {
                        break;
                    }

                    for (uint64_t i = 0; i < buffer.size(); i++)
                    {
                        uint8_t c = (uint8_t)buffer[i];
                        if (prevChar != c || x == 0)
                        {
                            run_bits.push_bit(true);
                            run_bits.push_bit(false);
                            this->head_char_vec[currentRunP++] = c;
                            prevChar = c;
                        }
                        else
                        {
                            run_bits.push_bit(false);

                            //run_bits.push_back(0);
                        }
                        x++;
                    }
                }

                run_bits.push_bit(true);
                run_bits.finish();

                std::cout << "BWT using memory = " << sdsl::size_in_bytes(head_char_vec) / 1000 << "[KB]" << std::endl;
                std::cout << "Run bits using memory = " << run_bits.get_using_memory() / 1000 << "[KB]" << std::endl;

                inp.close();

                analysisResult.print();
                this->lpos_vec.build_from_builder(run_bits);
            }
            CHAR get_smallest_character() const
            {
                return this->smallest_character;
            }

            uint64_t get_lindex_containing_the_position(uint64_t lposition) const
            {
                uint64_t x = this->lpos_vec.rank(lposition + 1) - 1;
                return x;
            }
            INDEX get_run(INDEX i) const
            {
                return (lpos_vec)[(i + 1)] - (lpos_vec)[i];
            }
            INDEX get_lpos(INDEX i) const
            {
                return (lpos_vec)[i];
            }
            INDEX rle_size() const
            {
                return (head_char_vec).size();
            }

            INDEX str_size() const
            {
                return (lpos_vec)[(lpos_vec).size() - 1];
            }
            INDEX get_end_rle_lposition() const
            {
                for (INDEX i = 0; i < head_char_vec.size(); i++)
                {
                    if ((head_char_vec)[i] == this->smalles_character)
                    {
                        return i;
                    }
                }
                return std::numeric_limits<INDEX>::max();
            }
            /*
            void print_info() const
            {
                std::string msg = "";
                for (INDEX i = 0; i < this->rle_size(); i++)
                {
                    if ((head_char_vec)[i] == 0)
                    {
                        msg += "#";
                    }
                    else
                    {
                        msg.push_back((char)(*this->head_char_vec)[i]);
                    }

                    msg += "^";
                    msg += std::to_string(this->get_run(i));
                    if (i + 1 < this->rle_size())
                        msg += ", ";
                }
                std::cout << "RLBWT: " << msg << std::endl;
                std::cout << "This instance has two vectors: head_char_vec and lpos_vec." << std::endl;
                std::cout << "The first array stores " << std::endl;
                for (auto c : *this->head_char_vec)
                {
                    std::cout << c << " ";
                }
                std::cout << std::endl;
                //stool::Printer::print(*this->head_char_vec);
                std::cout << "The second array stores " << std::endl;
                for (auto c : *this->lpos_vec)
                {
                    std::cout << c << " ";
                }
                std::cout << std::endl;

                //stool::Printer::print(*this->lpos_vec);
            }

            

            INDEX get_lindex_containing_the_position(INDEX lposition) const
            {
                auto p = std::upper_bound(this->lpos_vec->begin(), this->lpos_vec->end(), lposition);
                INDEX pos = std::distance(this->lpos_vec->begin(), p) - 1;
                return pos;
            }
            INDEX get_lindex_containing_the_position_with_linear_search(INDEX lposition, INDEX start_lindex) const
            {
                uint64_t rleSize = this->lpos_vec->size();
                for (uint64_t i = start_lindex; i < rleSize; i++)
                {
                    uint64_t pos = (lpos_vec)[i];
                    if (lposition < pos)
                    {
                        return i - 1;
                    }
                }
                throw std::logic_error("ERROR!");
            }
            */
            /*
    uint64_t get_using_memory() const {
        uint64_t m = 0;
        if(is_same<CHARVEC, std::vector<CHAR>>::value ){
            m += this->rle_size() * sizeof(CHAR);
        }else{
            m += this->head_char_vec->get_using_memory();
        }
        if(is_same<POWVEC, std::vector<INDEX>>::value){
            m += this->rle_size() * sizeof(INDEX);

        }else{
            m += this->lpos_vec->get_using_memory();            
        }
        return m;

    }
    */
            /*
    stool::rlbwt::ForwardBWT<CHAR, INDEX, CHARVEC,POWVEC> get_bwt(){
        return stool::rlbwt::ForwardBWT<CHAR, INDEX, CHARVEC,POWVEC>(this->head_char_vec, this->lpos_vec);
    }
    */
        };

    } // namespace stnode_on_rlbwt
} // namespace stool