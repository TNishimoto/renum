#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <deque>
#include <vector>
#include <type_traits>
#include "../weiner_link_emulator.hpp"
#include "fast_stnode_sub_traverser.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        template <typename INDEX_SIZE, typename RLBWTDS>
        class FastSTNodeTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using STNODE_SUB_TRAVERSER = FastSTNodeSubTraverser<INDEX_SIZE, RLBWTDS>;

            std::vector<STNODE_SUB_TRAVERSER *> sub_trees;
            std::vector<std::vector<STNODE_SUB_TRAVERSER *>> new_trees;
            

            std::vector<ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> ems;
            uint64_t minimum_child_count = 1000;
            uint64_t sub_tree_limit_size = 2000;

            uint64_t current_lcp = 0;
            uint64_t _child_count = 0;
            uint64_t _node_count = 0;

            uint64_t thread_count = 1;
            std::stack<uint64_t> position_stack;
            uint64_t kk = 0;
#if DEBUG
            uint64_t prev_child_count = 0;
#endif

        public:
            RLBWTDS *_RLBWTDS;

            uint64_t get_current_lcp() const
            {
                return current_lcp;
            }
            uint64_t child_count() const
            {
                return this->_child_count;
            }

            uint64_t node_count() const
            {
                return this->_node_count;
            }

            void initialize(uint64_t size, RLBWTDS &__RLBWTDS)
            {
                this->_RLBWTDS = &__RLBWTDS;

                this->thread_count = size;

                if (this->thread_count == 1)
                {
                    this->sub_tree_limit_size = UINT64_MAX;
                }

                //assert(size == 1);

                //this->sub_trees.resize(size);

                auto st = new STNODE_SUB_TRAVERSER(this->_RLBWTDS);
                sub_trees.push_back(st);
                //sub_trees.resize(256);

                this->new_trees.resize(size);


                for (uint64_t i = 0; i < this->thread_count; i++)
                {
                    ems.push_back(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>());
                    ems[ems.size() - 1].initialize(&__RLBWTDS);
                }
            }
            uint64_t get_stnode(uint64_t L, stool::LCPInterval<uint64_t> &output)
            {
                uint64_t p = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    uint64_t csize = this->sub_trees[i]->children_count();

                    if (p <= L && L < p + csize)
                    {
                        return p + this->sub_trees[i]->get_stnode(L - p, output, this->current_lcp - 1);
                    }
                    else
                    {
                        p += csize;
                    }
                }
                return UINT64_MAX;
            }

            void process()
            {
                //this->heavyEnumerate();
                this->recompute_node_counter();
            }
            bool isStop()
            {
                return this->current_lcp > 0 && this->child_count() == 0;
                //return total_counter == strSize - 1;
            }

            void print()
            {
                std::cout << "PRINT PTREE" << std::endl;
                //this->sub_tree.print_info();

                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    this->sub_trees[i]->print_info();
                }

                std::cout << "[END]" << std::endl;
            }

            uint64_t get_using_memory()
            {
                uint64_t k = 0;

                for (auto &it : this->sub_trees)
                {
                    k += it->get_using_memory();
                }

                return k;
            }

        private:
            void recompute_node_counter()
            {

                uint64_t current_child_count = 0;
                uint64_t current_node_count = 0;
                for (auto &it : this->sub_trees)
                {
                    current_child_count += it->children_count();
                    current_node_count += it->node_count();
                }
                _node_count = current_node_count;
                _child_count = current_child_count;

                current_lcp++;
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool