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

touch $LDIR/$1/targets.txt
for library_name in "archaea"  "bacteria"  "viral"  "fungi"  "plant"  "human"  "protozoa"
do
	if [ -e $LDIR/$1/$library_name.txt ]; then
		echo "Loading $library_name.txt"
		cat $LDIR/$1/$library_name.txt >> $LDIR/$1/targets.txt
	fi
done


$LDIR/buildDB $LDIR/$1/targets.txt $LDIR/$1/nameFamily.txt $LDIR/$1/database_full.txt
$LDIR/uniqueDB $LDIR/$1/database_full.txt $LDIR/$1/database
rm $LDIR/$1/database_full.txt
