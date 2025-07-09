#include <iostream>
#include <fstream>
#include <cstdlib>
#include <stdio.h>
#include <vector>
#include <random>
#include <chrono>
#include "libdivsufsort/sa.hpp"
#include "stool/include/stool.hpp"
#include <divsufsort64.h>





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


bool load(std::string &filename, std::vector<uint8_t> &output)
{
    std::ifstream file;
    file.open(filename, std::ios::binary);

    if (!file)
    {
        std::cerr << "error reading file " << std::endl;
        return false;
    }
    file.seekg(0, std::ios::end);
    auto n = (unsigned long)file.tellg();
    file.seekg(0, std::ios::beg);

    output.resize(n / sizeof(char));

    file.read((char *)&(output)[0], n);
    file.close();
    file.clear();

    return true;
}
bool write(std::string filename, std::vector<uint8_t> &text)
{
    std::cout << "writing: " << filename << std::endl;
    auto start = std::chrono::system_clock::now();

    std::ofstream os(filename, std::ios::out | std::ios::binary);
    if (!os)
        return 1;
    os.write((const char *)(&text[0]), sizeof(char) * text.size());
    os.close();

    auto end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "wrote: " << filename << ", time = " << elapsed << std::endl;

    return true;
}
void constructOccCharVec(const std::vector<uint8_t> &text, std::vector<bool> &output){
    output.resize(256, false);
    for(uint64_t i=0;i<text.size();i++){
        output[text[i]] = true;
    }
}
void replace(std::vector<uint8_t> &text, uint8_t oldChar, uint8_t newChar){
    for(uint64_t i=0;i<text.size();i++){
        if(text[i] == oldChar){
            text[i] = newChar;
        }
    }
}
void sanityze(std::vector<uint8_t> &text){
    std::vector<bool> occVec;
    constructOccCharVec(text, occVec);
    if(occVec[0]){
        std::cout << "This text contains character 0" << std::endl;
        uint8_t replaceChar = 0;
        for(uint64_t i=1;i<occVec.size();i++){
            if(i != 8 && !occVec[i]){
                replaceChar = i;
                break;
            }
        }
        if(replaceChar != 0){
            std::cout << "We replace the character 0 with " << (int)replaceChar << "." << std::endl;
            replace(text, 0, replaceChar);
        }else{
            std::cout << "We cannot replace the character 0."<< std::endl;
            throw -1;
        }

    }

    if(occVec[8]){
        std::cout << "This text contains character 8" << std::endl;
        uint8_t replaceChar = 8;
        for(uint64_t i=9;i<occVec.size();i++){
            if(!occVec[i]){
                replaceChar = i;
                break;
            }
        }
        if(replaceChar != 8){
            std::cout << "We replace the character 8 with " << (int)replaceChar << "." << std::endl;
            replace(text, 8, replaceChar);
        }else{
            std::cout << "We cannot replace the character 8."<< std::endl;
            throw -1;
        }

    }

}


int main(int argc, char *argv[])
{
    cmdline::parser p;
    p.add<std::string>("input_file", 'i', "input file path", true);
    p.add<std::string>("output_file", 'o', "output bwt file path", false, "");
    p.add<std::string>("special_character", 's', "special character", false, "");

    p.parse_check(argc, argv);
    std::string inputFile = p.get<std::string>("input_file");
    std::string outputFile = p.get<std::string>("output_file");
    std::string specialCharacter = p.get<std::string>("special_character");

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
    sanityze(text);
    /*
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
    */

    text.push_back(sc);
    //std::vector<uint64_t> sa = libdivsufsort::construct_suffix_array(text);
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
