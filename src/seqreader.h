/**
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#
# The reader sequence class is referred from Kraken2
# Copyright 2013-2019, Derrick Wood <dwood@cs.jhu.edu>
*/

#ifndef SEQREADER_H_
#define SEQREADER_H_

#include <bits/stdc++.h>
using namespace std;

class Reader{
    /// 0 = FASTA, 1 = FASTQ
private:
      stringstream ss_;
      string str_buffer_;  // used to prevent realloc upon every load/parse
      int FASTA = 0, FASTQ = 1;
      char *block_buffer_;
      size_t block_buffer_size_;

public:
    int file_format_;
    Reader();
    ~Reader();
    Reader(const Reader &rhs) = delete;
    Reader& operator=(const Reader &rhs) = delete;


    bool LoadBlock(istream &ifs, size_t block_size);
    bool NextSequence(string &headID, string &seq);
    bool ReadNextSequence(std::istream &is, string &headID, string &seq, std::string &str_buffer, int file_format);
};

#endif
