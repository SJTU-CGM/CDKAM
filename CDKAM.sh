#!/bin/bash
#
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#

DBDR=$1
RANK=0
FSCRPT=$(readlink -f "$0")
LDIR=$(dirname "$FSCRPT")

if [ $# -lt 4 ]; then
	echo "Usage: $0 <DTB path> <input> <output> <--fasta/--fastq>"
	exit
fi

$LDIR/classify $LDIR/$1/database $LDIR/$1/nameFamily.txt $2 $3 $4

