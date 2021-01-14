#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <vector>
#include <random>
#include <chrono>
#include "divsufsort64.h"
#include "stool/src/cmdline.h"


using namespace std;

std::vector<uint8_t> construct_bwt(const std::vector<uint8_t> &text)
{
    std::vector<uint8_t> bwt;
    std::vector<uint64_t> sa;

    uint64_t n = text.size();
    sa.resize(n);

    divsufsort64((const unsigned char *)&text[0], (int64_t *)&sa[0], n);

    bwt.resize(n);
    for (uint64_t i = 0; i < bwt.size(); i++)
    {
        if (sa[i] == 0)
        {
            bwt[i] = text[n - 1];
        }
        else
        {
            bwt[i] = text[sa[i] - 1];
        }
    }

    return bwt;
}

bool load(string &filename, std::vector<uint8_t> &output)
{
    std::ifstream file;
    file.open(filename, std::ios::binary);

    if (!file)
    {
        std::cerr << "error reading file " << endl;
        return false;
    }
    file.seekg(0, ios::end);
    auto n = (unsigned long)file.tellg();
    file.seekg(0, ios::beg);

    output.resize(n / sizeof(char));

    file.read((char *)&(output)[0], n);
    file.close();
    file.clear();

    return true;
}
bool write(string filename, std::vector<uint8_t> &text)
{
    std::cout << "writing: " << filename << std::endl;
    auto start = std::chrono::system_clock::now();

    ofstream os(filename, ios::out | ios::binary);
    if (!os)
        return 1;
    os.write((const char *)(&text[0]), sizeof(char) * text.size());
    os.close();

    auto end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "wrote: " << filename << ", time = " << elapsed << std::endl;

    return true;
}

int main(int argc, char *argv[])
{
    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file path", true);
    p.add<string>("output_file", 'o', "output bwt file path", false, "");
    p.add<string>("special_character", 's', "special character", false, "");

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    string outputFile = p.get<string>("output_file");
    string specialCharacter = p.get<string>("special_character");

    uint8_t sc = 0;
    if(specialCharacter.size() > 0){
        sc = specialCharacter[0];
    }
    std::cout << "Special character: " << sc << "(" << (uint)sc << ")" << std::endl;
    if (outputFile.size() == 0)
    {
            outputFile = inputFile + ".bwt";
    }
    std::vector<uint8_t> text;
    load(inputFile, text);
    bool xb = false;
    for (uint64_t i = 0; i < text.size(); i++)
    {
        if (text[i] == 0)
        {
            xb = true;
        }
    }
    if (xb)
    {
        std::cout << "Character Error!" << std::endl;
        throw -1;
    }

    text.push_back(sc);
    std::vector<uint8_t> bwt = construct_bwt(text);

    if(bwt.size() < 100){
        std::cout << "Input Text: ";
        for(auto c : text){
            std::cout << c;
        }
        std::cout << std::endl;

        std::cout << "Output BWT: ";
        for(auto c : bwt){
            std::cout << c;
        }
        std::cout << std::endl;
    }

    write(outputFile, bwt);
}
