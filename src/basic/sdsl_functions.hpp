#pragma once

#include "./fmindex.hpp"
#include <sdsl/wt_algorithm.hpp>
#include <sdsl/wavelet_trees.hpp>
namespace stool
{
    namespace stnode_on_rlbwt
    {
        class SDSLFunction
        {
        public:
            static int deleteFile(string fileName)
            {
                std::cout << "Delete: " << fileName << std::endl;
                return !(remove(fileName));
            }
            static uint8_t get_last_char(std::string inputFile)
            {
                std::ifstream stream;
                stream.open(inputFile, std::ios::binary);

                if (!stream)
                {
                    std::cerr << "error reading file " << std::endl;
                    throw -1;
                }
                stream.seekg(0, std::ios::end);
                uint64_t n = (unsigned long)stream.tellg();
                stream.seekg(0, std::ios::beg);
                uint8_t c = 0;

                uint64_t k = 1;
                while (c == 0)
                {
                    stream.seekg(n - k, std::ios::beg);
                    stream.read((char *)&c, 1);
                    k++;
                }
                stream.close();

                return c;
            }
            static uint8_t load_wavelet_tree(string bwt_iv_file, wt_huff<> &wt, std::vector<uint64_t> &output_C_array)
            {
                uint8_t c = get_last_char(bwt_iv_file);

                //std::cout << "Constructing Wavelet Tree..." << std::endl;
                construct(wt, bwt_iv_file);
                //std::cout << "WT using memory = " << sdsl::size_in_bytes(wt) / 1000 << "[KB]" << std::endl;

                stool::FMIndex::constructCArray(wt, c, output_C_array);
                deleteFile(bwt_iv_file);
                return c;
            }
            static void constructDBitArray(std::string inputFile, sdsl::bit_vector &bv)
            {
                std::ifstream stream;
                stream.open(inputFile, std::ios::binary);

                if (!stream)
                {
                    std::cerr << "error reading file " << std::endl;
                    throw -1;
                }
                stream.seekg(0, std::ios::end);
                uint64_t n = (unsigned long)stream.tellg();
                stream.seekg(0, std::ios::beg);
                uint8_t c = 0;

                uint64_t k = 0;
                while (c == 0)
                {
                    k++;
                    stream.seekg(n - k, std::ios::beg);
                    stream.read((char *)&c, 1);
                }
                uint64_t textSize = (n - k + 1) - 9;
                stream.seekg(9, std::ios::beg);


                bv.resize(textSize);

                std::vector<char> buffer;
                uint64_t bufferSize = 81920;
                buffer.resize(8192);
                uint8_t prevChar = 255;
                uint64_t x = 0;
                bool bitChar = false;

                while (true)
                {
                    bool b = stool::FileReader::read(stream, buffer, bufferSize, textSize + 9);
                    if (!b)
                    {
                        break;
                    }

                    for (uint64_t i = 0; i < buffer.size(); i++)
                    {
                        uint8_t c = (uint8_t)buffer[i];

                        if (prevChar != c || x == 0)
                        {
                            prevChar = c;
                            bitChar = !bitChar;
                        }
                        bv[x++] = (uint8_t)bitChar;
                    }
                }

                stream.close();
            }
        };

    } // namespace stnode_on_rlbwt

} // namespace stool