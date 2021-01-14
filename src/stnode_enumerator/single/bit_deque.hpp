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

//#include "range_distinct/range_distinct_on_rlbwt.hpp"

namespace stool
{
    namespace stnode_on_rlbwt
    {
        class BitDeque
        {
            std::deque<uint8_t> deq;
            uint8_t fst_value = 0;
            uint8_t last_value = 0;

            uint64_t fst_unused_size = 8;
            uint64_t last_unused_size = 8;

        private:
            void back_to_front()
            {
                assert(this->fst_value == 0);
                assert(this->deq.size() == 0);
                assert(this->fst_unused_size == 8);
                this->fst_value = this->last_value << this->last_unused_size;
                this->fst_unused_size = this->last_unused_size;
                this->last_value = 0;
                this->last_unused_size = 8;

            }
            void front_to_back()
            {
                assert(this->last_value == 0);
                assert(this->deq.size() == 0);
                assert(this->last_unused_size == 8);
                this->last_value = this->fst_value >> this->fst_unused_size;
                this->last_unused_size = this->fst_unused_size;
                this->fst_value = 0;
                this->fst_unused_size = 8;
            }

        public:
            void clear()
            {
                this->deq.clear();
                this->fst_value = 0;
                this->last_value = 0;
                this->fst_unused_size = 8;
                this->last_unused_size = 8;
            }
            bool operator[](uint64_t i) const
            {
                //uint64_t fstw = 8 - start;
                uint64_t j = (i + fst_unused_size) / 8;
                uint64_t diff = (i + fst_unused_size) - (j * 8);
                uint8_t v = 0;
                if (j == 0)
                {
                    v = this->fst_value;
                }
                else if ((j - 1) == this->deq.size())
                {
                    v = this->last_value;
                }
                else
                {
                    v = this->deq[j - 1];
                }
                //std::cout << "i: " << i << "[" << j << ", " << diff << "]" << ((v >> diff) & 1) << "/" << ((uint64_t)v) << std::endl;
                return (v >> diff) & 1;
            }
            uint64_t size() const
            {
                return this->deq.size() * 8 + (8 - fst_unused_size) + (8 - last_unused_size);
            }
            void push_back(bool value)
            {
#if DEBUG
                uint64_t prevSize = this->size();
#endif
                if (last_unused_size == 0)
                {
                    this->deq.push_back(this->last_value);
                    this->last_value = 0;
                    last_unused_size = 8;
                }
                if (value)
                {
                    uint8_t p = this->last_value;
                    p = p | (1 << (8 - last_unused_size));
                    this->last_value = p;
                }
                last_unused_size--;
                assert(this->size() == prevSize + 1);
            }
            void pop_back()
            {
#if DEBUG
                uint64_t prevSize = this->size();
#endif

                if (last_unused_size == 8)
                {
                    if (this->deq.size() > 0)
                    {
                        this->last_value = this->deq[this->deq.size() - 1];
                        this->deq.pop_back();
                        last_unused_size = 0;
                    }
                    else
                    {
                        if (this->fst_unused_size != 8)
                        {
                            this->front_to_back();
                        }
                        else
                        {
                            std::cout << "Empty Error!" << std::endl;
                            assert(false);
                            throw -1;
                        }
                    }
                }
                uint8_t p = this->last_value;
                p = p << (last_unused_size + 1);
                p = p >> (last_unused_size + 1);
                this->last_value = p;
                last_unused_size++;
                assert(this->size() == prevSize - 1);
            }
            void push_front(bool value)
            {
#if DEBUG
                uint64_t prevSize = this->size();
#endif

                if (fst_unused_size == 0)
                {
                    this->deq.push_front(this->fst_value);
                    this->fst_value = 0;
                    fst_unused_size = 8;
                }
                if (value)
                {
                    uint8_t p = this->fst_value;
                    p = p | (1 << (fst_unused_size - 1));
                    this->fst_value = p;
                }
                fst_unused_size--;
                assert(this->size() == prevSize + 1);
            }
            void pop_front()
            {
#if DEBUG
                uint64_t prevSize = this->size();
#endif

                if (fst_unused_size == 8)
                {
                    if (this->deq.size() > 0)
                    {
                        this->fst_value = this->deq[0];
                        this->deq.pop_front();
                        fst_unused_size = 0;
                    }
                    else
                    {
                        if (this->last_unused_size != 8)
                        {
                            this->back_to_front();
                        }
                        else
                        {
                            std::cout << "Empty Error!" << std::endl;
                            assert(false);
                            throw -1;
                        }
                    }
                }
                uint8_t p = this->fst_value;
                p = p >> (fst_unused_size + 1);
                p = p << (fst_unused_size + 1);
                this->fst_value = p;
                fst_unused_size++;

                assert(this->size() == prevSize - 1);
            }
            std::vector<bool> to_vector()
            {
                std::vector<bool> output;
                for (uint64_t i = 0; i < this->size(); i++)
                {
                    output.push_back((*this)[i]);
                }
                return output;
            }
            std::vector<uint64_t> to_vector2()
            {
                std::vector<uint64_t> output;
                for (uint64_t i = 0; i < this->deq.size(); i++)
                {
                    output.push_back(this->deq[i]);
                }
                return output;
            }
        };
    } // namespace stnode_on_rlbwt
} // namespace stool
