/**
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#
# The reader sequence class is referred from Kraken2
# Copyright 2013-2019, Derrick Wood <dwood@cs.jhu.edu>
*/

#include "seqreader.h"
#include <sysexits.h>
#include <err.h>
#include <fcntl.h>

void StripString(string &str) {
  while (isspace(str.back()))
    str.pop_back();
}

Reader::Reader(){
    str_buffer_.reserve(8192);
    block_buffer_ = new char[8192];
    block_buffer_size_ = 8192;
};

Reader::~Reader(){
    delete[] block_buffer_;
};


bool Reader::LoadBlock(istream &ifs, size_t block_size){
    ss_.clear();
    ss_.str("");
    if (block_buffer_size_ < block_size) {
        delete[] block_buffer_;
        block_buffer_ = new char[block_size];
        block_buffer_size_ = block_size;
    }
    ifs.read(block_buffer_, block_size);
    if (! ifs && ifs.gcount() <= 0)
        return false;

    str_buffer_.assign(block_buffer_, ifs.gcount());
    ss_ << str_buffer_;
    if (getline(ifs, str_buffer_)){
        ss_ << str_buffer_ << "\n";
    }
    if (file_format_ == FASTQ) {
        while (getline(ifs, str_buffer_)) {
        ss_ << str_buffer_ << "\n";
        if (str_buffer_[0] == '@')
            break;
        }
        int lines_to_read = 0;
        if (getline(ifs, str_buffer_)) {
            ss_ << str_buffer_ << "\n";
            lines_to_read = str_buffer_[0] == '@' ? 3 : 2;
            while (lines_to_read-- > 0 && getline(ifs, str_buffer_))
                ss_ << str_buffer_ << "\n";
        }
    }
    else {
        while (ifs) {
            if (ifs.peek() == '>')
                break;
            if (getline(ifs, str_buffer_)){
                ss_ << str_buffer_ << "\n";
            }
        }
    }
    return true;
}

bool Reader::NextSequence(string &headID, string &seq){
    return ReadNextSequence(ss_, headID, seq, str_buffer_, file_format_);
}

bool Reader::ReadNextSequence(std::istream &is, string &headID, string &seq, std::string &str_buffer, int file_format){
    if (! getline(is, str_buffer))
        return false;
    StripString(str_buffer);

    auto first_whitespace_ch = str_buffer.find_first_of(" \t\r", 1);
    auto substr_len = first_whitespace_ch;
    if (substr_len != std::string::npos)
        substr_len--;
    if (str_buffer.size() > 1)
        headID.assign(str_buffer, 1, substr_len);
    else
        return false;


    if (file_format == FASTQ) {
        if (! getline(is, str_buffer))
            return false;
        StripString(str_buffer);
        seq.assign(str_buffer);
        if (! getline(is, str_buffer))  //  + line, discard
            return false;
        if (! getline(is, str_buffer))
            return false;
    }
    else if (file_format == FASTA) {
        seq.assign("");
        while (is && is.peek() != '>') {
            if (! getline(is, str_buffer))
                return ! seq.empty();
        StripString(str_buffer);
        seq.append(str_buffer);
        }
    }
    return true;
}
