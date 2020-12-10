#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "stool/src/cmdline.h"
#include "common.hpp"
#include "libdivsufsort/sa.hpp"
#include "../postorder_maximal_substring_intervals.hpp"
#include "../forward_bwt.hpp"
using namespace std;
using CHAR = char;
using INDEX = uint64_t;

uint64_t input_text_size = 0;
std::vector<std::pair<std::string, uint64_t>> execution_time_messages;

std::vector<char> getInputText(string filepath)
{
    vector<char> text;

    std::ifstream ifs(filepath);
    bool inputFileExist = ifs.is_open();
    if (!inputFileExist)
    {
        std::cout << filepath << " cannot open." << std::endl;
        throw - 1;
    }

    std::string tmp;
    //vector<char> result;
    std::cout << "Loading the input text..." << std::endl;

    while (getline(ifs, tmp))
    {
        if(tmp.size() >= 2){
            std::string pref = tmp.substr(0,2);
            uint64_t startPos = 0; 
            if(pref == "1 "){
                startPos = 2;
            }else if(pref == "-1"){
                startPos = 3;
            }
            for(uint64_t i=startPos;i<tmp.size();i++ ){
                text.push_back(tmp[i]);
            }
        }
        text.push_back((char)1);

    }
    text.pop_back();
    text.push_back((char)0);

    return text;
}

std::vector<INDEX> createLineOccList(std::vector<CHAR> &text){
    std::vector<INDEX> r;
    r.push_back(0);
    for(uint64_t i=0;i<text.size();i++){
        if(text[i] == (char)1){
            r.push_back(i);
        }
    }
    return r;
}
std::vector<INDEX> createSampledLineOccList(std::vector<CHAR> &text){
    uint64_t counter=0;
    std::vector<INDEX> r;
    for(uint64_t i=0;i<text.size();i++){
        if(i % 64 == 0){
            r.push_back(counter);
        }
        if(text[i] == 1){
            counter++;
        }
    }
    return r;
}


std::pair<uint64_t, uint64_t> getLineOcc(uint64_t pos, std::vector<INDEX> &lineOccList, std::vector<INDEX> &sampledLineOccList){
        assert(pos != 0);
        uint64_t index = (pos / 64);
        uint64_t pindex = sampledLineOccList[index];

            assert(lineOccList[pindex] < pos);
            for(uint64_t x=pindex;x<lineOccList.size()-1;x++){
                if(lineOccList[x] < pos && pos < lineOccList[x+1]){
                    return std::pair<uint64_t, uint64_t>(x, pos - (lineOccList[x]+1) );
                }
            }
            return std::pair<uint64_t, uint64_t>(lineOccList.size(), pos - (lineOccList[lineOccList.size()-1]+1) );

        /*
        auto p = std::lower_bound(lineOccList.begin(), lineOccList.end(), pos);
        INDEX r = distance(lineOccList.begin(), p);
        uint64_t lineStartPos = 0;
        if(r != 0){
            lineStartPos = lineOccList[r-1]+1;
        }
        */
}

std::vector<pair<uint64_t,uint64_t>> getOccurrences(stool::LCPInterval<INDEX> &interval, std::vector<INDEX> &sa, std::vector<INDEX> &lineOccList, std::vector<INDEX> &sampledLineOccList){
    std::vector<pair<uint64_t,uint64_t>> r;
    for(uint64_t x=interval.i;x<=interval.j;x++){
        std::pair<uint64_t, uint64_t> p = getLineOcc(sa[x], lineOccList, sampledLineOccList);
        r.push_back(p);
    }
    return r;
}


bool lineFilter(stool::LCPInterval<INDEX> &interval, std::vector<INDEX> &sa, std::vector<CHAR> &text){
    for(uint64_t i=0;i<interval.lcp;i++){
        if(text[sa[interval.i]+i] == (char)1){
            return false;            
        }
    }
    return true;
}


std::string getInfo(stool::LCPInterval<INDEX> &interval, std::vector<INDEX> &sa, std::vector<CHAR> &text, std::vector<INDEX> &lineOccList){
    std::string ms = "";
    for(uint64_t i=0;i<interval.lcp;i++){
        ms.push_back(text[sa[interval.i]+i]);
    }
    return ms;
}

uint64_t iterateMS(std::vector<CHAR> &T, std::ofstream &out)
{
    //vector<CHAR> T = stool::load_char_vec_from_file(filename, true); // input text
    input_text_size = T.size();

    auto start_sa = std::chrono::system_clock::now();
    std::vector<INDEX> sa = stool::construct_suffix_array(T);
    auto end_sa = std::chrono::system_clock::now();
    double sa_construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_sa - start_sa).count();
    execution_time_messages.push_back(std::pair<std::string, uint64_t>("SA construction time\t\t", sa_construction_time));

    using BWT = stool::esaxx::ForwardBWT<CHAR, std::vector<CHAR>, std::vector<INDEX>>;
    BWT bwt(&T, &sa);

    auto start_lcp = std::chrono::system_clock::now();
    std::vector<INDEX> lcpArray = stool::constructLCP<CHAR, INDEX>(T, sa);
    auto end_lcp = std::chrono::system_clock::now();
    double lcp_array_construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_lcp - start_lcp).count();
    execution_time_messages.push_back(std::pair<std::string, uint64_t>("LCP array construction time\t", lcp_array_construction_time));

    auto start_ms = std::chrono::system_clock::now();
    stool::esaxx::PostorderMaximalSubstringIntervals<CHAR, INDEX, std::vector<INDEX>, BWT> pmsi;
    pmsi.construct(&lcpArray, &bwt);

    std::vector<INDEX> lineOccList = createLineOccList(T);
    std::vector<INDEX> sampledLineOccList = createSampledLineOccList(T);

    uint64_t count = 0;
    for (auto it : pmsi)
    {
        //out.write(reinterpret_cast<const char *>(&it), sizeof(stool::LCPInterval<INDEX>));
        if(lineFilter(it, sa, T)){
            if(it.lcp == 0) continue;
            std::string msstr = getInfo(it, sa, T, lineOccList);

            std::vector<pair<uint64_t,uint64_t>> occList = getOccurrences(it,sa,lineOccList, sampledLineOccList);
            std::string listStr="";
            std::sort(occList.begin(), occList.end(), [](const pair<uint64_t,uint64_t>& lhs, const pair<uint64_t,uint64_t>& rhs) {
                return lhs.first < rhs.first;});
            for(auto p : occList){
                std::string p1 = std::to_string(p.first+1);
                //std::string p2 = std::to_string(p.second);
                for(auto c : p1){
                    listStr.push_back(c);
                }
                //listStr.push_back('#');
                //for(auto c : p2){
                //    listStr.push_back(c);
                //}

                listStr.push_back(' ');
            }

            //out.write(msstr);
            out << msstr << std::endl;
            out << listStr << std::endl;

		    //out.write((const char *)&msstr[0], sizeof(char) * msstr.size());

            ++count;
            if(count % 10000 == 0) std::cout << "+" << std::flush;
        }
    }
            if(count % 10000 == 0) std::cout << std::endl;

    auto end_ms = std::chrono::system_clock::now();
    double ms_construction_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_ms - start_ms).count();
    execution_time_messages.push_back(std::pair<std::string, uint64_t>("MS Construction time\t\t", ms_construction_time));
    return count;
}

int main(int argc, char *argv[])
{

    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file name", true);
    p.add<string>("output_file", 'o', "output file name", false, "");
    //p.add<bool>("print", 'p', "print info", false, true);
    //p.add<string>("format", 'f', "output format (binary or csv)", false, "binary");
    //p.add<string>("mode", 'm', "mode(esaxx, rlbwt or succinct)", false, "esaxx");
    //p.add<bool>("memory", 'u', "using only main memory (0 or 1)", false, 1);

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    string outputFile = p.get<string>("output_file");
    //string format = p.get<string>("format");
    //string mode = p.get<string>("mode");
    //bool usingMemory = p.get<bool>("memory");


    //bool isPrint = p.get<bool>("print");

    if (outputFile.size() == 0)
    {
        outputFile = inputFile + ".max.log";
    }

    auto start = std::chrono::system_clock::now();

    std::ofstream out(outputFile, std::ios::out | std::ios::binary);
    if (!out)
    {
        throw std::runtime_error("Cannot open the output file!");
    }
    uint64_t ms_count = 0;
    std::vector<stool::LCPInterval<INDEX>> intervals;
    //mode = "non-compressed";


    vector<char> text = getInputText(inputFile); // input text
    /*
    for(uint64_t i=0;i<text.size();i++){
        std::cout << text[i] << std::flush;
        if(text[i] == 1){
            std::cout << std::endl;
        }
    }
    */


    ms_count = iterateMS(text, out);

    auto end = std::chrono::system_clock::now();
    double elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    std::cout << "\033[31m";
    std::cout << "______________________RESULT______________________" << std::endl;
    std::cout << "File \t\t\t\t\t : " << inputFile << std::endl;
    std::cout << "Output \t\t\t\t\t : " << outputFile << std::endl;
    std::cout << "The length of the input text \t\t : " << input_text_size << std::endl;
    std::cout << "The number of maximum substrings \t : " << ms_count << std::endl;
    std::cout << "Excecution time \t\t\t : " << elapsed << "[ms]" << std::endl;
    for (auto it : execution_time_messages)
    {
        std::cout << "|\t " << it.first << " : " << it.second << "[ms]" << std::endl;
    }

    std::cout << "_______________________________________________________" << std::endl;
    std::cout << "\033[39m" << std::endl;
}