#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <stdexcept>
#include <iostream>
#include <fstream>
//#include "common/io.h"
//#include "common/print.hpp"
//#include "other_functions.hpp"
//#include "OnlineRlbwt/online_rlbwt.hpp"
//#include "rlbwt.hpp"
//#include "stool/src/elias_fano_vector.hpp"
#include "stool/src/io.h"

namespace stool
{
    namespace rlbwt2
    {
        struct BWTAnalysisResult{
            uint64_t run_count;
            uint64_t str_size;

            uint64_t min_char = UINT64_MAX;
            uint64_t min_char_pos;
            uint64_t min_char_count;

        };
        static BWTAnalysisResult count_runs(std::string filename)
        {
            BWTAnalysisResult result;
            std::ifstream inp;
            std::vector<char> buffer;
            uint64_t bufferSize = 8192;
            buffer.resize(8192);

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
            uint64_t count_run = 0;
            while (true)
            {
                bool b = stool::FileReader::read(inp, buffer, bufferSize, textSize);
                if (!b)
                {
                    break;
                }

                for (uint64_t i = 0; i < buffer.size(); i++)
                {
                    uint8_t c = buffer[i];
                    if(c < result.min_char){
                        result.min_char = c;
                        result.min_char_pos = x;
                        result.min_char_count = 1;
                    }else if(c == result.min_char){
                        result.min_char_count++;
                    }

                    if (prevChar != c || x == 0)
                    {
                        count_run++;
                        prevChar = c;
                    }
                    x++;
                }
            }
            inp.close();

            result.run_count = count_run;
            result.str_size = x;
            return result;
        }

        static void load_RLBWT_from_file(std::string filename, sdsl::int_vector<> &diff_char_vec, std::vector<bool> &run_bits)
        {
            diff_char_vec.width(8);

            std::ifstream inp;
            std::vector<char> buffer;
            uint64_t bufferSize = 8192;
            buffer.resize(8192);

            BWTAnalysisResult analysisResult = count_runs(filename);
            diff_char_vec.resize(analysisResult.run_count);
            if(analysisResult.min_char != 0){
                std::cout << "We replace the character " << (int)analysisResult.min_char << "(" << (char)analysisResult.min_char << ") with 0" << std::endl;
            }

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
                    assert(buffer[i] >=0);
                    uint8_t c = buffer[i];
                    if(c == analysisResult.min_char){
                        c = 0;
                    }
                    if (prevChar != c || x == 0)
                    {
                        run_bits.push_back(1);
                        run_bits.push_back(0);

                        diff_char_vec[currentRunP++] = c;
                        prevChar = c;
                    }
                    else
                    {
                        run_bits.push_back(0);
                    }
                    x++;
                }
            }

            run_bits.push_back(1);

            inp.close();
        }
        static std::vector<uint64_t> construct_lpos_array(std::vector<bool> &run_bits)
        {
            std::vector<uint64_t> r;
            uint64_t prev_i = 0;
            for (uint64_t i = 0; i < run_bits.size(); i++)
            {
                if (run_bits[i])
                {
                    if(r.size() == 0){
                        r.push_back(i);
                    }else{
                        uint64_t diff = i - prev_i - 1;
                        uint64_t x = r[r.size()-1] + diff;
                        r.push_back(x);
                    }
                    prev_i = i;
                }
            }
            return r;
        }
    } // namespace rlbwt2
} // namespace stool