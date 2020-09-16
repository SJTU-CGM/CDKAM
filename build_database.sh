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
	echo "Building the database."
	echo "Usage: $0 <DTB path> "
	exit
fi

if [ -e $LDIR/$1/targets.txt ]; then
	rm $LDIR/$1/targets.txt
fi
touch $LDIR/$1/targets.txt

for library_name in "archaea"  "bacteria"  "viral"  "fungi"  "plant"  "human"  "protozoa"
do
	if [ -e $LDIR/$1/$library_name.txt ]; then
		echo "Loading $library_name.txt"
		cat $LDIR/$1/$library_name.txt >> $LDIR/$1/targets.txt
	fi
done

if [ ! -e $LDIR/$1/taxanomy/nodes.dmp ]; then
	echo "The taxonomy tree is missing. Please run the command   ./download_taxonomy.sh $1"
	exit 0
fi

echo "The first step: collecting kmers."
$LDIR/buildDB $LDIR/$1/targets.txt $LDIR/$1/nameFamily.txt $LDIR/$1/database_full.txt

if [ ! -e $LDIR/$1/database_full.txt ]; then
	echo "The first step is failed."
	exit 0
fi

echo "The second step: solving kmers collision."
$LDIR/uniqueDB $LDIR/$1/database_full.txt $LDIR/$1/database
rm $LDIR/$1/database_full.txt
echo "Building database: DONE."
