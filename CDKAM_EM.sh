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
	echo "Usage: $0 DBname input output --fasta/--fastq"
	echo "Or"
	echo "Usage: $0 DBname input output --fasta/--fastq nthread N"
	exit
fi

for library_name in "database_Taxo"  "database_Suffix"  "database_Size"
do
	if [ ! -e $LDIR/$1/$library_name ]; then
		echo "Database does not contain necessary file $library_name"
		exit 0
	fi
done

$LDIR/classifyEM $LDIR/$1/database $LDIR/$1/nameFamily.txt $2 $3 $4 $5 $6

