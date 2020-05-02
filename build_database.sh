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

if [ $# -lt 1 ]; then
	echo "Usage: $0 <DTB path> "
	exit
fi

$LDIR/buildDB $LDIR/$1/targets.txt $LDIR/$1/nameFamily.txt $LDIR/$1/database_full.txt
$LDIR/uniqueDB $LDIR/$1/database_full.txt $LDIR/$1/database
rm $LDIR/$1/database_full.txt