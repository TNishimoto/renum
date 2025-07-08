#include <cassert>
#include <chrono>
#include <string>

#include "stool/include/stool.hpp"
#include "../stnode_enumerator/single/bit_deque.hpp"

std::vector<bool> to_vector(std::deque<bool> &vec)
{
    std::vector<bool> r;
    for (auto &it : vec)
    {
        r.push_back(it);
    }
    return r;
}

bool check_test(stool::renum::BitDeque &vec1, std::deque<bool> &vec2)
{
    /*
    auto r = vec1.to_vector();
    stool::Printer::print(r);
    auto r2 = to_vector(vec2);
    stool::Printer::print(r2);
    */

    if (vec1.size() != vec2.size())
    {
        std::string s = std::string("String sizes are different!") + ", collect = " + std::to_string(vec1.size()) + ", test = " + std::to_string(vec2.size());

        throw std::logic_error(s);
    }
    for (uint64_t i = 0; i < vec1.size(); i++)
    {
        if (vec1[i] != vec2[i])
        {
            std::string msg = "collect_vec[" + std::to_string(i) + "] != test_vec[" + std::to_string(i) + "]";

            throw std::logic_error("Values are different! " + msg);
        }
    }
}

int main(int argc, char *argv[])
{
    /*
  cmdline::parser p;
  p.add<std::string>("input_file", 'i', "input file name", true);
  p.add<std::string>("mode", 'm', "mode", false, "xx");
  p.add<std::string>("output_file", 'o', "output file name", false, "");

  p.parse_check(argc, argv);
  std::string inputFile = p.get<std::string>("input_file");
  std::string mode = p.get<std::string>("mode");
  std::string outputFile = p.get<std::string>("output_file");
  std::string format = "binary";
  */

    stool::renum::BitDeque bdeq;

    /*
    bdeq.push_back(true);
    bdeq.push_back(false);
    bdeq.pop_front
    */

    std::cout << "Test[push_back]" << std::endl;
    for (uint64_t i = 0; i < 20; i++)
    {
        bdeq.push_back(i % 2 == 0);
        auto r = bdeq.to_vector();
        stool::Printer::print(r);
    }
    std::cout << "Test[pop_back]" << std::endl;
    for (uint64_t i = 0; i < 20; i++)
    {
        bdeq.pop_back();
        auto r = bdeq.to_vector();
        stool::Printer::print(r);
    }
    std::cout << "Test[push_front]" << std::endl;
    for (uint64_t i = 0; i < 20; i++)
    {
        bdeq.push_front(i % 2 == 0);
        auto r = bdeq.to_vector();
        stool::Printer::print(r);
    }
    std::cout << "Test[pop_front]" << std::endl;
    for (uint64_t i = 0; i < 20; i++)
    {
        bdeq.pop_front();
        auto r = bdeq.to_vector();
        stool::Printer::print(r);
    }

    /*
    std::cout << "Test[pop_back]" << std::endl;
    for(uint64_t i=0;i<20;i++){
        bdeq.pop_back();
        auto r = bdeq.to_vector();
        stool::Printer::print(r);        
    }
    */

    stool::renum::BitDeque deq1;
    std::deque<bool> deq2;
    for (uint64_t i = 0; i < 10; i++)
    {
        deq1.push_back(i % 2 == 0);
        deq2.push_back(i % 2 == 0);
    }
    auto seq = stool::create_deterministic_integers<int64_t>(100, 4, -4, 100);

    for (auto &it : seq)
    {
        auto r2 = to_vector(deq2);
        stool::Printer::print(r2);
        if (it >= 0)
        {
            if (it % 2 == 1)
            {
                deq1.push_back(it <= 2);
                deq2.push_back(it <= 2);
            }
            else
            {
                deq1.push_front(it <= 2);
                deq2.push_front(it <= 2);
            }
        }
        else if (deq2.size() > 0)
        {
            if (it % 2 == 0)
            {
                deq1.pop_front();
                deq2.pop_front();
            }
            else
            {
                deq1.pop_back();
                deq2.pop_back();
            }
        }
        check_test(deq1, deq2);
    }
    std::cout << "OK!" << std::endl;
}
