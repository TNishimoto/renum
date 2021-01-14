#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>
//#include "next_rinterval_storage_constructor.hpp"
#include "../rlbwt/rle_wavelet_tree.hpp"

#include "suffix_tree_nodes.hpp"
#include "explicit_weiner_link_computer_on_rlbwt.hpp"
#include <thread>
#include "../debug/stnode_checker.hpp"

namespace stool
{
    namespace stnode_on_rlbwt
    {
        class STTreeAnalysisResult
        {
        public:
            uint64_t max_nodes_at_level;
            uint64_t print_interval = 100;

            uint64_t position_counter = 0;
            uint64_t input_text_length = 0;
            int64_t current_lcp = -1;
            uint64_t current_child_count = 0;
            uint64_t print_interval_counter = 0;
            std::chrono::system_clock::time_point start_time;

            void print_info()
            {
                if (this->current_lcp == 0)
                {
                    std::cout << "["
                              << "x"
                              << "/" << this->print_interval << "] ";
                    std::cout << "Current LCP";
                    std::cout << ", The peak number of children";
                    std::cout << ", The number of children with the current depth";
                    std::cout << ", Processed characters: ["
                              << "y"
                              << "/" << this->input_text_length << "]";
                    std::cout << ", Elapsed time [s]";
                    std::cout << ", The number of processed characters per second [KB/s]";

                    std::cout << std::endl;
                }
                std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
                double elapsed = std::chrono::duration_cast<std::chrono::seconds>(end - start_time).count();
                double aveg = 0;
                if(elapsed > 0){
                    aveg = ((double)this->position_counter / elapsed) / (1000);
                }

                std::cout << "[" << (this->print_interval_counter) << "/" << this->print_interval << "] ";
                std::cout << "LCP = " << this->current_lcp;
                std::cout << ", " << this->max_nodes_at_level;
                std::cout << ", " << this->current_child_count;
                std::cout << ", [" << this->position_counter << "/" << this->input_text_length << "]";
                std::cout << ", " << elapsed << "[s]";
                std::cout << ", " << aveg << "[KB/s]";

                std::cout << std::endl;
            }
            void start()
            {
                this->start_time = std::chrono::system_clock::now();
            }

            template <typename STNODES>
            void analyze(STNODES &stnodeSequencer)
            {

                this->current_lcp++;
                this->position_counter += stnodeSequencer.child_count() - stnodeSequencer.node_count();
                this->current_child_count = stnodeSequencer.child_count();
                if (this->position_counter > 0)
                {
                    bool b = false;
                    while (true)
                    {
                        uint64_t x = (this->input_text_length / this->print_interval) + 1;
                        uint64_t pp_num = x * this->print_interval_counter;
                        if (this->position_counter >= pp_num)
                        {
                            b = true;
                            this->print_interval_counter++;
                        }
                        else
                        {
                            break;
                        }
                    }
                    if (b)
                    {
                        this->print_info();
                    }
                }

                if (this->max_nodes_at_level < stnodeSequencer.child_count())
                {
                    this->max_nodes_at_level = stnodeSequencer.child_count();
                }
            }
        };

        class Application
        {
        public:
            template <typename STNODES>
            static uint64_t outputMaximalSubstrings(std::ofstream &out, STNODES &stnodeSequencer, STTreeAnalysisResult &analysis)
            {
                analysis.start();
                
                using INDEX_SIZE = typename STNODES::index_type;
                uint8_t print_type = 0;
                out.write((char *)(&print_type), sizeof(print_type));
                uint8_t index_bits = sizeof(INDEX_SIZE);
                out.write((char *)(&index_bits), sizeof(index_bits));

                uint64_t count = 0;

                analysis.input_text_length = stnodeSequencer.get_input_text_length();

                std::vector<stool::LCPInterval<INDEX_SIZE>> buffer;

                auto depth_iterator = stnodeSequencer.begin();
                while (depth_iterator != stnodeSequencer.end())
                {
                    analysis.analyze(stnodeSequencer);
                    std::vector<stool::LCPInterval<uint64_t>> r2;

                    for (auto node_it = depth_iterator.begin(); node_it != depth_iterator.end(); node_it++)
                    {
                        bool b = stnodeSequencer.check_maximal_repeat(node_it);
                        if (b)
                        {
                            INDEX_SIZE left = stnodeSequencer.get_left(node_it);
                            INDEX_SIZE right = stnodeSequencer.get_right(node_it);
                            stool::LCPInterval<INDEX_SIZE> intv(left, right, stnodeSequencer.get_current_lcp());
                            buffer.push_back(intv);
                            count++;
                            if (buffer.size() >= 8000)
                            {
                                out.write((char *)(&buffer[0]), sizeof(stool::LCPInterval<INDEX_SIZE>) * buffer.size());
                                buffer.clear();
                            }
                        }
                    }
                    depth_iterator++;
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

    } // namespace stnode_on_rlbwt
} // namespace stool