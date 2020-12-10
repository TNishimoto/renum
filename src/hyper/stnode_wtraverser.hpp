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
            std::vector<bool> maximal_repeat_check_vec;

            STNodeWTraverser()
            {
            }
            RINTERVAL& get_stnode(uint64_t i){
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
            uint64_t children_count() const {
                return this->childVec.size();
            }
            uint64_t node_count() const {
                return this->stnodeVec.size();
            }

            void swap(STNodeWTraverser &copy)
            {
                this->stnodeVec.swap(copy.stnodeVec);
                this->childVec.swap(copy.childVec);
                this->widthVec.swap(copy.widthVec);

            }
            void clear()
            {
                this->stnodeVec.clear();
                this->childVec.clear();
                this->widthVec.clear();

            }
            void first_compute(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                em.computeFirstLCPIntervalSet();
                em.move_st_internal_nodes(this->stnodeVec, this->childVec, this->widthVec);

                this->maximal_repeat_check_vec.resize(1);
                this->maximal_repeat_check_vec[0] = true;

                //this->move_from(em);
            }
            void computeNextLCPIntervalSet(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &inputSet, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS> &em)
            {
                assert(inputSet.node_count() > 0);

                this->clear();
                uint64_t rank = 0;
                for (uint64_t i = 0; i < inputSet.node_count(); i++)
                {
                    em.computeNextLCPIntervalSet(inputSet.stnodeVec[i], inputSet.childVec, inputSet.children_count(), rank);
                    //this->move_from(em);
                    em.move_st_internal_nodes(this->stnodeVec, this->childVec, this->widthVec);

                    rank += inputSet.widthVec[i];
                }
            }
            void spill(STNodeWTraverser<INDEX_SIZE, RLBWTDS> &item, uint64_t limit_child_count){
                assert(item.children_count() <= limit_child_count);

                uint64_t capacity = limit_child_count - item.children_count();
                uint64_t k = 0;
                uint64_t k2 = 0;
                int64_t x = this->node_count()-1;
                while(x >= 0 && k <= capacity && (this->children_count() - k) > limit_child_count ){
                    k += this->widthVec[x--];
                    k2++;
                }
                uint64_t current_children_count = this->children_count();
                for(uint64_t i = current_children_count - k;i < current_children_count;i++){
                    item.childVec.push_back(this->childVec[i]);
                }
                for(uint64_t i = current_children_count - k;i < current_children_count;i++){
                    this->childVec.pop_back();
                }

                uint64_t current_node_count = this->node_count();

                for(uint64_t i = current_node_count - k2;i < current_node_count;i++){
                    item.stnodeVec.push_back(this->stnodeVec[i]);
                    item.widthVec.push_back(this->widthVec[i]);
                }
                for(uint64_t i = current_node_count - k2;i < current_node_count;i++){
                    this->stnodeVec.pop_back();
                    this->widthVec.pop_back();
                }

            }
            void print(){
                std::cout << "[" << this->node_count() << ", " << this->children_count() << "]" << std::endl;
                
            }
            void print2(const RLBWTDS &ds){
                for(uint64_t i=0;i<this->node_count();i++){
                    this->stnodeVec[i].print2(ds._fposDS);
                }
            }


        private:
        };

    } // namespace lcp_on_rlbwt
} // namespace stool