#!/bin/bash
#
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#

if [ $# -ne 1 ]; then
	echo "Usage: $0 <Directory: directory to store taxonomy data> "
	echo "Note: if the chosen directory is not empty, then its content will be erased."
	exit
fi

if  [ -d $1/taxonomy ] ; then
	exit
fi

mkdir -m 775 $1/taxonomy
cd $1/taxonomy
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

# Extrat downloaded data
if  [ -s taxdump.tar.gz ] ; then
	tar -zxf taxdump.tar.gz
	if [ -s nodes.dmp ]; then
		exit
	else
		echo "Failed to uncompress taxonomy data."
	fi
else
	echo "Failed to download taxonomy data!"
	exit
fi

