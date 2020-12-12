#pragma once
#include <cassert>
#include <chrono>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <queue>
#include <vector>
#include <type_traits>
#include "weiner_link_emulator.hpp"
//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {
        template <typename INDEX_SIZE, typename RLBWTDS>
        class STNodeWTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            std::vector<RINTERVAL> stnodeVec;
            std::vector<RINTERVAL> childVec;
            std::vector<uint8_t> widthVec;

        public:
            RLBWTDS *_RLBWTDS;

            std::vector<bool> maximal_repeat_check_vec;

            STNodeWTraverser()
            {
            }
            RINTERVAL &get_stnode(uint64_t i)
            {
                return this->stnodeVec[i];
            }

            std::vector<RINTERVAL> *get_stnode_vec()
            {
                return &this->stnodeVec;
            }

            std::vector<RINTERVAL> *get_child_vec()
            {
                return &this->childVec;
            }

            /*
            uint64_t size(){
                return this->stnodeVec.size();
            }
            */
            uint64_t children_count() const
            {
                return this->childVec.size();
            }
            uint64_t node_count() const
            {
                return this->stnodeVec.size();
            }
            uint64_t get_width(uint64_t i) const
            {
                return this->widthVec[i];
            }
            /*
            void swap(STNodeWTraverser &copy)
            {
                this->stnodeVec.swap(copy.stnodeVec);
                this->childVec.swap(copy.childVec);
                this->widthVec.swap(copy.widthVec);
            }
            */
            void clear()
            {
                this->stnodeVec.clear();
                this->childVec.clear();
                this->widthVec.clear();
                this->maximal_repeat_check_vec.clear();
            }
            void first_compute(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                em.computeFirstLCPIntervalSet();
                em.move_st_internal_nodes(this->stnodeVec, this->childVec, this->widthVec);

                #if DEBUG
                    this->_RLBWTDS->checkLCPInterval(this->stnodeVec[0]);
                #endif

                this->maximal_repeat_check_vec.resize(1);
                this->maximal_repeat_check_vec[0] = true;

                //this->move_from(em);
            }
            void computeNextLCPIntervalSet(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &inputSet, uint64_t start_index, uint64_t width, uint64_t rank, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                assert(inputSet.node_count() > 0);

                this->clear();
                uint64_t end_index = start_index + width - 1;
                //uint64_t k = 0;
                assert(end_index  < inputSet.stnodeVec.size());

                for (uint64_t i = start_index; i <= end_index; i++)
                {

                    em.computeNextLCPIntervalSet(inputSet.stnodeVec[i], inputSet.childVec, rank, inputSet.widthVec[i]);
                    //this->move_from(em);
                    em.move_st_internal_nodes(this->stnodeVec, this->childVec, this->widthVec);
                    //assert(k == inputSet.widthVec[i]);

                    rank += inputSet.widthVec[i];
                }
            }
            /*
            void spill(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &item, uint64_t limit_child_count)
            {
                assert(item.children_count() <= limit_child_count);

                uint64_t capacity = limit_child_count - item.children_count();
                uint64_t k = 0;
                uint64_t k2 = 0;
                int64_t x = this->node_count() - 1;
                while (x >= 0 && k <= capacity && (this->children_count() - k) > limit_child_count)
                {
                    k += this->widthVec[x--];
                    k2++;
                }
                uint64_t current_children_count = this->children_count();
                for (uint64_t i = current_children_count - k; i < current_children_count; i++)
                {
                    item.childVec.push_back(this->childVec[i]);
                }
                for (uint64_t i = current_children_count - k; i < current_children_count; i++)
                {
                    this->childVec.pop_back();
                }

                uint64_t current_node_count = this->node_count();

                for (uint64_t i = current_node_count - k2; i < current_node_count; i++)
                {
                    item.stnodeVec.push_back(this->stnodeVec[i]);
                    item.widthVec.push_back(this->widthVec[i]);
                }
                for (uint64_t i = current_node_count - k2; i < current_node_count; i++)
                {
                    this->stnodeVec.pop_back();
                    this->widthVec.pop_back();
                }
            }
            */
            void print()
            {
                std::cout << "[" << this->node_count() << ", " << this->children_count() << "]" << std::endl;
            }
            void print_info()
            {
                std::cout << "PRINT TREE" << std::endl;
                for (uint64_t i = 0; i < this->node_count(); i++)
                {
                    this->stnodeVec[i].print2(this->_RLBWTDS->_fposDS);
                }
                std::cout << std::endl;
            }

            void add(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &item)
            {
                assert(item.stnodeVec.size() == item.widthVec.size());
                assert(item.stnodeVec.size() == item.maximal_repeat_check_vec.size());

                assert(this->stnodeVec.size() == this->widthVec.size());
                assert(this->stnodeVec.size() == this->maximal_repeat_check_vec.size());

                for (uint64_t i = 0; i < item.node_count(); i++)
                {
                    this->stnodeVec.push_back(item.stnodeVec[i]);
#if DEBUG
                    if (this->_RLBWTDS->str_size() < 100)
                    {
                        std::cout << "Add: " << std::endl;
                        item.stnodeVec[i].print();
                    }
#endif
                    this->widthVec.push_back(item.widthVec[i]);
                    this->maximal_repeat_check_vec.push_back(item.maximal_repeat_check_vec[i]);
                }
                for (uint64_t i = 0; i < item.children_count(); i++)
                {
                    this->childVec.push_back(item.childVec[i]);
                }
                item.clear();
            }
            /*
            void check()
            {
                uint64_t x = 0;
                for (uint64_t i = 0; i < this->node_count(); i++)
                {
                    uint64_t lcpIntvBeginPos = this->_RLBWTDS->get_fpos(this->stnodeVec[i].beginIndex, this->stnodeVec[i].beginDiff);
                    uint64_t lcpIntvEndPos = this->_RLBWTDS->get_fpos(this->stnodeVec[i].endIndex, this->stnodeVec[i].endDiff);

                    for (uint64_t j = 0; j < this->widthVec[i]; j++)
                    {
                        uint64_t weinerBeginPos = this->_RLBWTDS->get_fpos(this->childVec[x].beginIndex, this->childVec[x].beginDiff);
                        uint64_t weinerEndPos = this->_RLBWTDS->get_fpos(this->childVec[x].endIndex, this->childVec[x].endDiff);
                        assert(lcpIntvBeginPos <= weinerBeginPos && weinerEndPos <= lcpIntvEndPos);

                        x++;
                    }
                }
            }
            */
            /*
            void check(uint64_t t)
            {
                uint64_t x = 0;
                for (uint64_t i = 0; i < t; i++)
                {
                    uint64_t lcpIntvBeginPos = this->_RLBWTDS->get_fpos(this->stnodeVec[i].beginIndex, this->stnodeVec[i].beginDiff);
                    uint64_t lcpIntvEndPos = this->_RLBWTDS->get_fpos(this->stnodeVec[i].endIndex, this->stnodeVec[i].endDiff);

                    for (uint64_t j = 0; j < this->widthVec[i]; j++)
                    {
                        uint64_t weinerBeginPos = this->_RLBWTDS->get_fpos(this->childVec[x].beginIndex, this->childVec[x].beginDiff);
                        uint64_t weinerEndPos = this->_RLBWTDS->get_fpos(this->childVec[x].endIndex, this->childVec[x].endDiff);
                        assert(lcpIntvBeginPos <= weinerBeginPos && weinerEndPos <= lcpIntvEndPos);

                        x++;
                    }
                }
            }
            */

        private:
        };

    } // namespace lcp_on_rlbwt
} // namespace stool