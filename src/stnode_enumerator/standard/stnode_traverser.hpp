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
#include "stnode_sub_traverser.hpp"
#include "../single/single_stnode_traverser.hpp"
#include "multi_thread.hpp"

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace lcp_on_rlbwt
    {

        template <typename INDEX_SIZE, typename RLBWTDS>
        class STNodeTraverser
        {
            using RINTERVAL = RInterval<INDEX_SIZE>;
            using STNODE_SUB_TRAVERSER = STNodeSubTraverser<INDEX_SIZE, RLBWTDS>;
            using ITERATOR = STNodeIterator<STNodeTraverser>;
            using DEPTH_ITERATOR = STDepthIterator<STNodeTraverser>;
            using CHAR = typename RLBWTDS::CHAR;

            std::vector<STNODE_SUB_TRAVERSER *> sub_trees;

            //std::vector<std::vector<STNODE_SUB_TRAVERSER *>> new_trees;

            std::vector<ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> ems;
            uint64_t minimum_child_count = 1000;
            uint64_t sub_tree_limit_size = 2000;

            int64_t current_lcp = -1;
            uint64_t _child_count = 0;
            uint64_t _node_count = 0;

            uint64_t thread_count = 1;
            bool store_edge_chars = false;
            std::stack<uint64_t> position_stack;
#if DEBUG
            uint64_t prev_child_count = 0;
#endif

        public:
            RLBWTDS *_RLBWTDS;
            using index_type = INDEX_SIZE;
            bool has_edge_characters() const
            {
                return this->store_edge_chars;
            }
            int64_t get_current_lcp() const
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
            void clear()
            {

                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    this->sub_trees[i]->clear();
                    delete this->sub_trees[i];
                    this->sub_trees[i] = nullptr;
                }
                this->sub_trees.clear();
                this->current_lcp = -1;
                this->_child_count = 0;
                this->_node_count = 0;
            }

            void initialize(uint64_t size, RLBWTDS &__RLBWTDS, bool _store_edge_chars)
            {
                this->store_edge_chars = _store_edge_chars;
                this->_RLBWTDS = &__RLBWTDS;

                this->thread_count = size;

                for (uint64_t i = 0; i < this->thread_count; i++)
                {
                    ems.push_back(ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>());
                    ems[ems.size() - 1].initialize(&__RLBWTDS);
                }
            }

            void load(std::ifstream &file)
            {
                assert(false);
                throw -1;
                /*
                this->clear();

                std::cout << "\033[33m";

                std::cout << "LOAD" << std::endl;
                file.seekg(0, std::ios::end);
                uint64_t textSize = (uint64_t)file.tellg() / sizeof(char);
                file.seekg(0, std::ios::beg);
                std::cout << "Data_SIZE " << textSize << std::endl;

                uint64_t sub_tree_count = 0;

                file.read((char *)&this->current_lcp, sizeof(uint64_t));
                file.read((char *)&sub_tree_count, sizeof(uint64_t));

                std::cout << "lcp = " << this->current_lcp << std::endl;
                std::cout << "sub_tree_count = " << sub_tree_count << std::endl;

                for (uint64_t i = 0; i < sub_tree_count; i++)
                {

                    auto st = new STNODE_SUB_TRAVERSER(this->sub_tree_limit_size, this->_RLBWTDS, this->store_edge_chars);
                    st->load(file);
                    this->sub_trees.push_back(st);
                    this->_child_count += st->children_count();
                    this->_node_count += st->node_count();
                }
                //this->print();
                std::cout << "\033[39m" << std::endl;
                */
            }
            void write(std::ofstream &out)
            {
                uint64_t sub_tree_size = this->sub_trees.size();
                uint64_t lcp = this->current_lcp;
                std::cout << "\033[32m";
                std::cout << "WRITE " << std::endl;

                std::cout << "Subtree count : " << sub_tree_size << std::endl;
                std::cout << "LCP : " << lcp << std::endl;

                out.write((char *)&lcp, sizeof(uint64_t));
                out.write((char *)&sub_tree_size, sizeof(uint64_t));

                for (auto &it : this->sub_trees)
                {
                    it->write(out);
                }
                out.close();
                //this->print();
                std::cout << "\033[39m" << std::endl;
            }

            uint64_t parallel_count()
            {

                std::vector<uint64_t> fst_pos_vec;
                std::vector<uint64_t> children_count_vec;

                fst_pos_vec.resize(this->thread_count, UINT64_MAX);
                children_count_vec.resize(this->thread_count, 0);

                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    if (i < this->thread_count)
                    {
                        fst_pos_vec[i] = i;
                    }
                    else
                    {
                        position_stack.push(i);
                    }
                }

                //auto start = std::chrono::system_clock::now();
                std::vector<thread> threads;

                for (uint64_t i = 0; i < this->thread_count; i++)
                {
                    threads.push_back(thread(parallel_count_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees), fst_pos_vec[i], ref(position_stack), ref(ems[i]), ref(children_count_vec[i])));
                }

                for (thread &t : threads)
                    t.join();

                assert(position_stack.size() == 0);
                uint64_t count = 0;
                for (auto &it : children_count_vec)
                {
                    count += it;
                }
                return count;
            }
            void get_lcp_intervals(std::vector<stool::LCPInterval<uint64_t>> &output)
            {
                for (auto &it : this->sub_trees)
                {
                    it->get_lcp_intervals(this->current_lcp, output);
                }
            }
            /*
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
            */

            bool succ()
            {
                if (this->is_finished())
                {
                    return false;
                }
                else
                {
                    if (this->current_lcp >= 0)
                    {
                        this->heavyEnumerate();
                        this->recompute_node_counter();
                    }
                    else
                    {
                        SingleSTNodeTraverser<INDEX_SIZE, ExplicitWeinerLinkEmulator<INDEX_SIZE, RLBWTDS>> tmp_traverser;
                        tmp_traverser.initialize(&ems[0], this->store_edge_chars);
                        tmp_traverser.succ();
                        STNodeVector<INDEX_SIZE> tmp;
                        tmp_traverser.convert_to_vector(tmp);
                        this->import(tmp_traverser.get_current_lcp(), tmp);
                    }
                    return true;
                }
            }
            bool is_finished() const
            {
                return this->current_lcp >= 0 && this->child_count() == 0;
                //return total_counter == strSize - 1;
            }

            void print()
            {
                std::cout << "PRINT STNODE Traverser" << std::endl;
                std::cout << "[LCP, STNODE_COUNT, CHILDREN_COUNT] = [" << this->current_lcp << ", " << this->_node_count << ", " << this->_child_count << "]" << std::endl;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    this->sub_trees[i]->print();
                }

                std::cout << "[END]" << std::endl;
            }

            uint64_t get_using_memory()
            {
                uint64_t k = 0;

                uint64_t p = this->sub_trees.size() * 32;

                for (auto &it : this->sub_trees)
                {
                    k += it->get_using_memory();
                }

                return p + k;
            }

            void to_stnode_vector(STNodeVector<INDEX_SIZE> &item)
            {
                for (auto &it : this->sub_trees)
                {
                    it->to_stnode_vector(item);
                }
            }
            void import(uint64_t lcp, STNodeVector<INDEX_SIZE> &item)
            {

                assert(this->sub_trees.size() == 0);
                auto st = new STNODE_SUB_TRAVERSER(this->sub_tree_limit_size, this->store_edge_chars);
                this->sub_trees.push_back(st);

                this->current_lcp = lcp;
                this->sub_trees[0]->import(item);
                this->_child_count = this->sub_trees[0]->children_count();
                this->_node_count = this->sub_trees[0]->node_count();
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

            void parallel_succ()
            {
#if DEBUG
                std::cout << "PARALLEL PROCESS" << std::endl;
#endif
                std::vector<uint64_t> fst_pos_vec;
                fst_pos_vec.resize(this->thread_count, UINT64_MAX);
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    if (i < this->thread_count)
                    {
                        fst_pos_vec[i] = i;
                    }
                    else
                    {
                        position_stack.push(i);
                    }
                }

                //auto start = std::chrono::system_clock::now();
                std::vector<thread> threads;

                for (uint64_t i = 0; i < this->thread_count; i++)
                {
                    threads.push_back(thread(parallel_succ_stnodes<INDEX_SIZE, RLBWTDS>, ref(sub_trees), fst_pos_vec[i], ref(position_stack), ref(ems[i]), this->sub_tree_limit_size));
                }

                for (thread &t : threads)
                    t.join();
                    //auto end = std::chrono::system_clock::now();
#if DEBUG
                std::cout << "PARALLEL PROCESS[END]" << std::endl;
#endif

                assert(position_stack.size() == 0);
            }
            void single_succ()
            {
                STNodeVector<INDEX_SIZE> tmp;
                std::queue<STNODE_SUB_TRAVERSER *> uqueue;
                for (auto &it : this->sub_trees)
                {
                    it->computeNextSTNodes(ems[0], tmp);
                    it->clear();
                    uqueue.push(it);
                    distribute(this->sub_trees, tmp, uqueue, this->sub_tree_limit_size);
                }
            }
            void heavyEnumerate()
            {
                bool isSingleProcess = false;
                /*
                if (current_lcp == 12)
                {
                    throw -1;
                }
                */

                if (current_lcp >= 0)
                {

                    isSingleProcess = this->child_count() < minimum_child_count || this->thread_count == 1;
                    //bool b = true;
                    if (isSingleProcess)
                    {
                        this->single_succ();
                    }
                    else
                    {

                        this->parallel_succ();
                    }
                }
                else
                {
                    assert(false);
                    throw -1;
                }

                this->remove_empty_trees();

                if ((double)this->sub_trees.size() * 2 < (double)this->sub_trees.capacity())
                {
                    this->sub_trees.shrink_to_fit();
                }
            }
            void remove_empty_trees()
            {
                uint64_t nonEmptyCount = 0;
                for (uint64_t i = 0; i < this->sub_trees.size(); i++)
                {
                    if (this->sub_trees[i]->node_count() > 0)
                    {
                        if (i != nonEmptyCount)
                        {

                            this->sub_trees[nonEmptyCount]->swap(*this->sub_trees[i]);

                            //this->sub_trees[i] = STNODE_SUB_TRAVERSER();
                            //sub_trees[i]._RLBWTDS = this->_RLBWTDS;
                        }
                        nonEmptyCount++;
                    }
                }
                for (uint64_t i = nonEmptyCount; i < this->sub_trees.size(); i++)
                {
                    this->sub_trees[i]->clear();
                    //this->recycle_sub_trees.push(this->sub_trees[i]);
                    delete this->sub_trees[i];
                    this->sub_trees[i] = nullptr;
                }
                this->sub_trees.resize(nonEmptyCount);
            }

        public:
            DEPTH_ITERATOR begin()
            {
                this->clear();

                this->succ();

                return DEPTH_ITERATOR(this, true);
            }
            DEPTH_ITERATOR end()
            {
                return DEPTH_ITERATOR(this, false);
            }

            void increment(ITERATOR &iter) const
            {
                this->increment(iter.child_index, iter.node_index, iter.array_index);
            }

            void increment(INDEX_SIZE &child_index, INDEX_SIZE &node_index, INDEX_SIZE &array_index) const
            {
                uint64_t left = 0, right = 0;
                child_index = this->sub_trees[array_index]->increment(child_index, left, right);

                if (child_index >= this->sub_trees[array_index]->get_integer_array_size())
                {
                    if (array_index + 1 < this->sub_trees.size())
                    {
                        child_index = 1;
                        node_index = 0;
                        array_index++;
                    }
                    else
                    {
                        child_index = std::numeric_limits<INDEX_SIZE>::max();
                        node_index = std::numeric_limits<INDEX_SIZE>::max();
                        array_index = std::numeric_limits<INDEX_SIZE>::max();
                    }
                }
                else
                {
                    node_index++;
                }
            }

            INDEX_SIZE get_children_count(const ITERATOR &iter) const
            {
                return this->get_children_count(iter.child_index, iter.array_index);
            }
            INDEX_SIZE get_children_count(INDEX_SIZE child_index, INDEX_SIZE array_index) const
            {
                uint64_t left = 0, right = 0;
                uint64_t R = this->sub_trees[array_index]->increment(child_index, left, right);
                return R - child_index - 1;
            }

            INDEX_SIZE get_left(const ITERATOR &iter) const
            {
                return this->get_left(iter.child_index, iter.array_index);
            }
            INDEX_SIZE get_left(INDEX_SIZE child_index, INDEX_SIZE array_index) const
            {

                uint64_t left = 0, right = 0;
                this->sub_trees[array_index]->increment(child_index, left, right);
                return left;
            }

            INDEX_SIZE get_right(const ITERATOR &iter) const
            {
                return this->get_right(iter.child_index, iter.array_index);
            }
            INDEX_SIZE get_right(INDEX_SIZE child_index, INDEX_SIZE array_index) const
            {

                uint64_t left = 0, right = 0;
                this->sub_trees[array_index]->increment(child_index, left, right);
                return right;
            }
            CHAR get_edge_character(uint64_t child_index, INDEX_SIZE array_index) const
            {
                assert(this->store_edge_chars);
                CHAR c = this->sub_trees[array_index]->get_edge_character(child_index);
                return c;
            }
            
            INDEX_SIZE get_edge_character(const ITERATOR &iter, uint64_t ith_child) const
            {
                return this->get_edge_character(iter.child_index + ith_child, iter.array_index);
            }
            INDEX_SIZE get_child_left_boundary(INDEX_SIZE child_index, INDEX_SIZE array_index) const
            {
                return this->sub_trees[array_index]->get_child_left_boundary(child_index);
            }
            INDEX_SIZE get_child_left_boundary(const ITERATOR &iter, uint64_t ith_child) const
            {
                return this->get_child_left_boundary(iter.child_index + ith_child, iter.array_index);
            }

            INDEX_SIZE get_child_right_boundary(INDEX_SIZE child_index, INDEX_SIZE array_index) const
            {
                return this->sub_trees[array_index]->get_child_left_boundary(child_index);
            }
            INDEX_SIZE get_child_right_boundary(const ITERATOR &iter, uint64_t ith_child) const
            {
                return this->get_child_right_boundary(iter.child_index + ith_child, iter.array_index);
            }

            void set_current_first_iterator(ITERATOR &iter) const
            {
                this->set_current_first_iterator(iter.child_index, iter.node_index, iter.array_index);
            }
            void set_current_first_iterator(INDEX_SIZE &child_index, INDEX_SIZE &node_index, INDEX_SIZE &array_index) const
            {
                if (this->is_finished())
                {
                    node_index = std::numeric_limits<INDEX_SIZE>::max();
                    child_index = std::numeric_limits<INDEX_SIZE>::max();
                    array_index = std::numeric_limits<INDEX_SIZE>::max();
                }
                else
                {
                    child_index = 1;
                    node_index = 0;
                    array_index = 0;
                }
            }

            bool check_maximal_repeat(INDEX_SIZE node_index, INDEX_SIZE array_index) const
            {
                return this->sub_trees[array_index]->check_maximal_repeat(node_index);
            }

            bool check_maximal_repeat(ITERATOR &iter) const
            {
                return this->check_maximal_repeat(iter.node_index, iter.array_index);
            }
        };

    } // namespace lcp_on_rlbwt
} // namespace stool