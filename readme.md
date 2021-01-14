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

This program outputs the suffix-tree intervals (a.k.a suffix-array interval) for all the maximal repeats 
in the string represented by a given BWT.  
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