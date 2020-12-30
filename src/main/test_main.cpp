// License: MIT http://opensource.org/licenses/MIT

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "../../module/stool/src/cmdline.h"
#include <thread>
#include <mutex>
#include <deque>
#include <queue>

//#include "../minimal_substrings/naive_minimal_substrings.hpp"

using namespace std;

using INDEXTYPE = int64_t;
using CHAR = char;

std::mutex mtx;
uint64_t u_value = 0;

struct Hoge
{
    std::deque<uint64_t> vec;
};

void parallel_process_deque(std::deque<uint64_t> *vec, uint64_t size)
{
    for (uint64_t i = 0; i < size; i++)
    {
        vec->push_back(i);
    }
}
void parallel_process_queue(std::queue<uint64_t> *vec, uint64_t size)
{
    for (uint64_t i = 0; i < size; i++)
    {
        vec->push(i);
    }
}
void parallel_process_vector(std::vector<uint64_t> *vec, uint64_t size)
{
    
    for (uint64_t i = 0; i < size; i++)
    {
        (*vec).push_back(i);
    }    
}
void parallel_process_structure(Hoge *vec, uint64_t size)
{
    
    for (uint64_t i = 0; i < size; i++)
    {
        vec->vec.push_back(i);
    }    
}

void test_deque(uint64_t size, uint64_t thread_count, bool is_parallel)
{
    std::cout << "TEST DEQUE " << (is_parallel ? "PARALLEL" : "SINGLE") << std::endl;

    std::vector<std::deque<uint64_t> *> tmp;
    for (uint64_t i = 0; i < thread_count; i++)
    {
        auto vec = new std::deque<uint64_t>();
        tmp.push_back(vec);
    }

    if (is_parallel)
    {
        std::vector<thread> threads;
        for (uint64_t i = 0; i < thread_count; i++)
        {
            threads.push_back(thread(parallel_process_deque, tmp[i], size));
        }

        for (thread &t : threads)
            t.join();
    }
    else
    {
        for (uint64_t i = 0; i < thread_count; i++)
        {
            parallel_process_deque(tmp[i], size);
        }
    }

    uint64_t y = 0;
    uint64_t z = 0;

    for (uint64_t i = 0; i < thread_count; i++)
    {
        y += tmp[i]->size() * 8;
        z += tmp[i]->size();
    }
    std::cout << z << "/" << (y / 1000) << " KB" << std::endl;
}
void test_vector(uint64_t size, uint64_t thread_count, bool is_parallel)
{
    std::cout << "TEST VECTOR " << (is_parallel ? "PARALLEL" : "SINGLE") << std::endl;

    std::vector<std::vector<uint64_t> *> tmp;
    for (uint64_t i = 0; i < thread_count; i++)
    {
        auto vec = new std::vector<uint64_t>();
    
        tmp.push_back(vec);
    }

    if (is_parallel)
    {
        std::vector<thread> threads;
        for (uint64_t i = 0; i < thread_count; i++)
        {
            threads.push_back(thread(parallel_process_vector, tmp[i], size));
        }
        
        for (thread &t : threads)
            t.join();
            
    }
    else
    {
        for (uint64_t i = 0; i < thread_count; i++)
        {
            parallel_process_vector(tmp[i], size);
        }
    }

    uint64_t y = 0;
    uint64_t z = 0;

    for (uint64_t i = 0; i < thread_count; i++)
    {
        y += tmp[i]->size() * 8;
        z += tmp[i]->size();
    }
    std::cout << z << "/" << (y / 1000) << " KB" << std::endl;
}
void test_queue(uint64_t size, uint64_t thread_count, bool is_parallel)
{
    std::cout << "TEST QUEUE " << (is_parallel ? "PARALLEL" : "SINGLE") << std::endl;
    std::vector<std::queue<uint64_t> *> tmp;

    for (uint64_t i = 0; i < thread_count; i++)
    {
        auto vec = new std::queue<uint64_t>();

        tmp.push_back(vec);
    }
    if (is_parallel)
    {
        std::vector<thread> threads;
        for (uint64_t i = 0; i < thread_count; i++)
        {
            threads.push_back(thread(parallel_process_queue, tmp[i], size));
        }

        for (thread &t : threads)
            t.join();
    }
    else
    {
        for (uint64_t i = 0; i < thread_count; i++)
        {
            parallel_process_queue(tmp[i], size);
        }
    }

    uint64_t y = 0;
    uint64_t z = 0;

    for (uint64_t i = 0; i < thread_count; i++)
    {
        y += tmp[i]->size() * 8;
        z += tmp[i]->size();
    }
    std::cout << z << "/" << (y / 1000) << " KB" << std::endl;
}
void test_hoge(uint64_t size, uint64_t thread_count, bool is_parallel)
{
    std::cout << "TEST HOGE " << (is_parallel ? "PARALLEL" : "SINGLE") << std::endl;
    std::vector<Hoge*> tmp;

    for (uint64_t i = 0; i < thread_count; i++)
    {
        auto vec = new Hoge();

        tmp.push_back(vec);
    }
    if (is_parallel)
    {
        std::vector<thread> threads;
        for (uint64_t i = 0; i < thread_count; i++)
        {
            threads.push_back(thread(parallel_process_structure, tmp[i], size));
        }

        for (thread &t : threads)
            t.join();
    }
    else
    {
        for (uint64_t i = 0; i < thread_count; i++)
        {
            parallel_process_structure(tmp[i], size);
        }
    }

    uint64_t y = 0;
    uint64_t z = 0;

    for (uint64_t i = 0; i < thread_count; i++)
    {
        y += tmp[i]->vec.size() * 8;
        z += tmp[i]->vec.size();
    }
    std::cout << z << "/" << (y / 1000) << " KB" << std::endl;
}
int main(int argc, char *argv[])
{

    cmdline::parser p;
    p.add<uint64_t>("size", 'i', "size", true, 1000);
    p.add<uint64_t>("thread", 't', "thread", true, 1);

    p.add<string>("mode", 'm', "mode", true, "1");

    //p.add<string>("tree_file", 't', "file type", false, "NULL");

    p.parse_check(argc, argv);
    uint64_t size = p.get<uint64_t>("size");
    string mode = p.get<string>("mode");
    uint64_t thread = p.get<uint64_t>("thread");

    uint64_t usize = size / thread;
    std::cout << "Size: " << size << std::endl;
    if (mode == "DS")
    {
        test_deque(usize, thread, false);
    }
    else if (mode == "DP")
    {
        test_deque(usize, thread, true);
    }
    else if (mode == "QS")
    {
        test_queue(usize, thread, false);
    }
    else if (mode == "QP")
    {
        test_queue(usize, thread, true);
    }
    else if (mode == "VS")
    {
        test_vector(usize, thread, false);
    }
    else if (mode == "VP")
    {
        test_vector(usize, thread, true);
    }
    else if (mode == "HS")
    {
        test_hoge(usize, thread, false);
    }
    else if (mode == "HP")
    {
        test_hoge(usize, thread, true);
    }
    else
    {
        throw -1;
    }

    return 0;
}