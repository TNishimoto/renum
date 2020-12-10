// License: MIT http://opensource.org/licenses/MIT
/*
  This code was copied from https://takeda25.hatenablog.jp/entry/20101202/1291269994 and I modified it.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "stool/src/cmdline.h"
//#include "../esa.hxx"
#include "../test/old_postorder_maximal_substrings.hpp"

using namespace std;
using INDEXTYPE = int64_t;

using LCPPair = std::pair<stool::LCPInterval<INDEXTYPE>, std::pair<uint64_t, uint64_t>>;

// Computes the kai squared value of the input values.
// The first return value is the kai squared value. 
// The second return value is true if a is larger than or equal to the expected value of a; otherwise it is false. 
std::pair<double, bool> get_kai_squared_value(double a, double b, double c, double d)
{
    double sum_ab = a + b;
    double sum_cd = c + d;
    double sum_ac = a + c;
    double sum_bd = b + d;
    double sum_abcd = a + b + c + d;

    double expected_a = sum_ab * (sum_ac / sum_abcd);
    double expected_b = sum_ab * (sum_bd / sum_abcd);
    double expected_c = sum_cd * (sum_ac / sum_abcd);
    double expected_d = sum_cd * (sum_bd / sum_abcd);

    double kai_a = ((a - expected_a) * (a - expected_a)) / expected_a;
    double kai_b = ((b - expected_b) * (b - expected_b)) / expected_b;
    double kai_c = ((c - expected_c) * (c - expected_c)) / expected_c;
    double kai_d = ((d - expected_d) * (d - expected_d)) / expected_d;
    double result = kai_a + kai_b + kai_c + kai_d;
    return std::pair<double, bool>(result, a >= expected_a);
}


template <typename sa_type, typename index_type>
string toCSVLine(LCPPair &item, vector<char> &text, sa_type &sa, index_type positive_sum, index_type negative_sum, string delimiter)
{
    string result = "";
    //index_type sum_y = item.second.first + item.second.second;
    //index_type plus_count_y = item.second.first;
    //double y_value = (double)plus_count_y / (double)sum_y;
    index_type plus_value = item.second.first;
    index_type minus_value = item.second.second;
    string line = item.first.getText(text, sa);

    index_type other_positive_sum = positive_sum - plus_value;
    index_type other_negative_sum = negative_sum - minus_value;

    std::pair<double, bool> r = get_kai_squared_value(plus_value, minus_value, other_positive_sum, other_negative_sum);
    string b1 = (r.second ? "+" : "-");

    result += std::to_string(r.first) + delimiter + b1 + delimiter + std::to_string(plus_value) + delimiter + std::to_string(minus_value) + delimiter + std::to_string(other_positive_sum) + delimiter + std::to_string(other_negative_sum) + delimiter + line;

    return result;
}

template <typename index_type>
std::pair<vector<char>, index_type> getInputText(string filepath)
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
        for (char c : tmp)
        {
            text.push_back(c);
        }
        text.push_back((char)1);
    }
    text.pop_back();
    text.push_back((char)0);

    index_type border = text.size();
    for (uint64_t i = 0; i < text.size(); i++)
    {
        if (text[i] == '-')
        {
            border = i;
            break;
        }
    }
    std::pair<vector<char>, index_type> r(std::move(text), border);
    return r;
}
INDEXTYPE getBorderPosition(vector<char> &text)
{
    vector<INDEXTYPE> result;
    for (uint64_t i = 0; i < text.size(); i++)
    {
        if (text[i] == '-')
        {
            return i;
        }
    }
    return text.size();
}
template <typename index_type>
vector<bool> constructForbiddenIndexesForMaximalSubstrings(vector<char> &text, vector<index_type> &sa, stool::PostorderMaximalSubstrings<index_type> &ms)
{
    vector<index_type> excludedPositions;
    excludedPositions.resize(text.size(), std::numeric_limits<index_type>::max());
    uint64_t prev = 0;
    uint64_t rank = 0;
    for (uint64_t i = 0; i < text.size(); i++)
    {
        if (text[i] == (char)1 || text[i] == '+' || text[i] == '-')
        {
            for (uint64_t j = prev; j < i; j++)
            {
                excludedPositions[j] = rank;
            }
            rank++;
            prev = i;
        }
    }

    vector<bool> result_vec;
    result_vec.resize(ms.size(), true);

    for(index_type i=0;i<ms.size();i++)
    {
        stool::LCPInterval<index_type> interval = ms[i];
        index_type start_pos = sa[interval.i];
        index_type end_pos = sa[interval.i] + interval.lcp - 1;

        bool b1, b2;

        if (start_pos == 0)
        {
            if (excludedPositions[0] == 0)
            {
                b1 = true;
            }
            else
            {
                b1 = false;
            }
        }
        else
        {
            b1 = (excludedPositions[start_pos] - excludedPositions[start_pos - 1]) == 0;
        }
        b2 = (excludedPositions[end_pos] - excludedPositions[start_pos]) == 0;
        result_vec[i] = b1 && b2 && (interval.lcp > 0);
    }
    return result_vec;
}

template <typename sa_type>
std::pair<uint64_t, uint64_t> getFrequency(stool::LCPInterval<INDEXTYPE> &interval, sa_type &sa, INDEXTYPE border)
{
    //vector<uint64_t> occurrences;
    INDEXTYPE plus = 0;
    INDEXTYPE minus = 0;
    for (INDEXTYPE i = interval.i; i <= interval.j; i++)
    {
        if (sa[i] < border)
        {
            plus++;
        }
        else
        {
            minus++;
        }
        //occurrences.push_back(sa[i]);
    }
    return std::pair<uint64_t, uint64_t>(plus, minus);
}

int main(int argc, char *argv[])
{

    cmdline::parser p;
    p.add<string>("input_file", 'i', "input file name", true);
    p.add<string>("output_file", 'o', "output file name", false, "");
    p.add<bool>("print", 'p', "print info", false, true);

    p.add<int64_t>("occ_threshold", 'a', "the threshold of occurrences", false, 0);
    p.add<int64_t>("len_threshold", 'b', "the threshold of length", false, 0);
    p.add<int64_t>("occlen_threshold", 'c', "the threshold of occ length", false, 0);

    p.parse_check(argc, argv);
    string inputFile = p.get<string>("input_file");
    string outputFile = p.get<string>("output_file");
    //string format = p.get<string>("format");
    int64_t occ_threshold = p.get<int64_t>("occ_threshold");
    int64_t len_threshold = p.get<int64_t>("len_threshold");
    int64_t occlen_threshold = p.get<int64_t>("occlen_threshold");

    //bool isPrint = p.get<bool>("print");

    if (outputFile.size() == 0)
    {
        outputFile = inputFile + "." + std::to_string(occ_threshold) + "." + std::to_string(len_threshold) + "." + std::to_string(occlen_threshold) + ".csv";
    }

    std::pair<vector<char>,INDEXTYPE> pr = getInputText<INDEXTYPE>(inputFile); // input text
    vector<char> &T = pr.first;
    //INDEXTYPE n = T.size();
    vector<INDEXTYPE> SA; // suffix array
    stool::PostorderMaximalSubstrings<INDEXTYPE> maximal_substrings = stool::PostorderMaximalSubstrings<INDEXTYPE>::construct(T, SA);
    vector<bool> filterVec = constructForbiddenIndexesForMaximalSubstrings(T, SA, maximal_substrings);

    ofstream os(outputFile, ios::out | ios::binary);
    if (!os)
        return 1;
    vector<LCPPair> buffer;
    INDEXTYPE positive_sum = 0;
    INDEXTYPE negative_sum = 0;

    for(INDEXTYPE i=0;i<maximal_substrings.size();i++)
    {
        stool::LCPInterval<INDEXTYPE> interval = maximal_substrings[i];
        if (!filterVec[i])
        {
            continue;
        }

        std::pair<int64_t, int64_t> freq = getFrequency(interval, SA, pr.second);
        positive_sum += freq.first;
        negative_sum += freq.second;
        if (freq.first + freq.second >= occ_threshold && interval.lcp >= len_threshold && (interval.lcp * (freq.first + freq.second) >= occlen_threshold))
        {
        }
        else
        {
            continue;
        }
        buffer.push_back(LCPPair(interval, freq));
    }

    std::sort(buffer.begin(), buffer.end(),
              [&](const LCPPair &x, const LCPPair &y) {
                  std::pair<double, bool> kai_x = get_kai_squared_value(x.second.first, x.second.second, positive_sum - x.second.first, negative_sum - x.second.second);
                  std::pair<double, bool> kai_y = get_kai_squared_value(y.second.first, y.second.second, positive_sum - y.second.first, negative_sum - y.second.second);
                  return kai_y.first < kai_x.first;
              });
    string column0 = "Kai_squared";
    string column1 = "Bias";
    string column2 = "MS(TRUE)";
    string column3 = "MS(FALSE)";
    string column4 = "Others(TRUE)";
    string column5 = "Others(FALSE)";
    string column6 = "Maximal_substring";

    os << column0 << "," << column1 << ", " << column2 << ", " << column3 << ", " << column4 << ", " << column5 << "," << column6 << std::endl;
    if(buffer.size() < 1000){
        std::cout << column0 << "\t" << column1 << "\t" << column2 << "\t" << column3 << "\t" << column4 << "\t" << column5 << "\t" << column6 << std::endl;
    }
    for (uint64_t z = 0; z < buffer.size(); z++)
    {
        LCPPair &item = buffer[z];
        string csvLine = toCSVLine(item, T, SA, positive_sum, negative_sum, ",");
        os << csvLine;
        if(buffer.size() < 1000){
            std::cout << toCSVLine(item, T, SA, positive_sum, negative_sum, "\t") << std::endl;
        }            
    }
    os << std::endl;
    buffer.clear();
    os.close();

    std::cout << "\033[36m";
    std::cout << "___________RESULT___________" << std::endl;
    std::cout << "File: " << inputFile << std::endl;
    std::cout << "Output: " << outputFile << std::endl;
    std::cout << "Occurrence threshold: " << occ_threshold << std::endl;
    std::cout << "Length threshold: " << len_threshold << std::endl;
    std::cout << "Occurrence*Length threshold: " << occlen_threshold << std::endl;

    std::cout << "The length of the input text: " << T.size() << std::endl;
    //std::cout << "The number of maximum substrings: " << maximumSubstringCount << std::endl;
    std::cout << "_________________________________" << std::endl;
    std::cout << "\033[39m" << std::endl;
}