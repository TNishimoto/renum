# Suffix Tree Node Iterator on RLBWT

//This library provides the implementation of a uni-directional iterator for suffix tree nodes.  
//The iterator uses $xr$ bytes and runs in $O(log (n/r))$ time per node,  
//where $n$ is the length of the input BWT, $r$ is the number of runs in the BWT, and $x$ is approximately 7~20.  
//This library also provides a program enumerating maximal repeats.  

## Download
The source codes in 'module' directory are maintained in different repositories. 
So, to download all the necessary source codes, do the following:

> git clone https://github.com/TNishimoto/stnode_iterator_on_rlbwt.git  
> cd stnode_iterator_on_rlbwt  
> git submodule init  
> git submodule update  

You need sdsl-lite(https://github.com/simongog/sdsl-lite) to excecute this program. Please edit CMakeLists.txt to set SDSL library and include directory paths.

## Compile
> mkdir build  
> cd build  
> cmake -DCMAKE_BUILD_TYPE=Release ..  
> make  

## Executions && Examples

### text_to_bwt.out  
This program outputs the BWT of a given text using libdivsufsort.  

usage: ./text_to_bwt.out --input_file=string [options] ...  
options:  
  -i, --input_file           input file path (string)  
  -o, --output_file          output bwt file path (string [=])  
  -s, --special_character    special character (string [=])  
  -?, --help                 print this message  
  
$ echo -n GATCAATGAGGTGGACACCAGAGGCGGGGACTTGT > sample.txt  
$ ./text_to_bwt.out -i sample.txt -s $  
Special character: $(36)  
Input Text: GATCAATGAGGTGGACACCAGAGGCGGGGACTTGT$  
Output BWT: TCGCGCGGGATACAGAGGAT$GTGAGCATGGAAGTC  
writing: sample.txt.bwt  
wrote: sample.txt.bwt, time = 2  

### maximal_repeat.out  

This program outputs the LCP intervals for all the maximal repeats 
in the string represented by a given BWT.  
The output file is econded in binary format, and print.out can decode the encoded file.  
The program can use multithreading for computing maximal repeats.  
It works with a single thread if you use "-p 1".  
It uses all the processors in your computer if you use "-p -1" or you do not use the -p option.  

usage: ./maximal_repeat.out --input_file=string [options] ...  
options:  
  -i, --input_file     input bwt file path (string)  
  -o, --output_file    output file path (default: input_file_path.max) (string [=])  
  -p, --thread_num     thread number for parallel processing (int [=-1])  
  -?, --help           print this message  

$ ./maximal_repeat.out -i ./sample.txt.bwt  
______________________RESULT______________________  
RLBWT File                               : ./sample.txt.bwt  
Output                                   : ./sample.txt.bwt.max  
Peak children count                      : 22  
Thread number                            : 8  
Integer Type                             : UINT32_t  
The length of the input text             : 36  
The number of runs on BWT                : 31  
The number of maximum substrings         : 16  
Excecution time                          : 0 [s]  
Character per second                     : inf [KB/s]  
         Preprocessing time              : 0 [s]  
         Enumeration time                : 0 [s]  
_______________________________________________________  

### print.out  
This program converts the LCP intervals stored in a given binary file into a new file in text format.  
It also outputs the string represented by each LCP interval if you use option -s 1.

usage: ./print.out --input_file=string --lcp_interval_file=string [options] ...   
options:  
  -i, --input_file           input bwt file path (string)  
  -l, --lcp_interval_file    LCP interval file path (string)  
  -o, --output_file          output file path (string [=])  
  -s, --string_flag          Output the string represented by each interval if this flag is 1 (bool [=1])  
  -?, --help                 print this message  

$ ./print.out -i ./sample.txt.bwt -l ./sample.txt.bwt.max -s 1  
(i, j, LCP, substring)  
0, 35, 0,  
30, 35, 1, T  
10, 15, 1, C  
16, 29, 1, G  
1, 9, 1, A  
28, 29, 2, GT  
8, 9, 2, AT  
2, 4, 2, AC  
22, 27, 2, GG  
5, 7, 2, AG  
32, 34, 2, TG  
10, 12, 2, CA  
16, 20, 2, GA  
25, 26, 3, GGG  
22, 23, 4, GGAC  
18, 19, 4, GAGG  
______________________RESULT______________________  
BWT File                                 : ./sample.txt.bwt  
Interval File                            : ./sample.txt.bwt.max  
Output File                              : ./sample.txt.bwt.max.interval.log  
File Type: LCPInterval  
Integer Type: uint32t  
Count: 16  
_______________________________________________________  