#!/bin/bash
#
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#

FSCRPT=$(readlink -f "$0")
LDIR=$(dirname "$FSCRPT")

echo |cpp -fopenmp -dM |grep -i open > $LDIR/.tmp
NB=`wc -l < $LDIR/.tmp`



g++ -o $LDIR/src/genSpecies $LDIR/src/genSpecies.cpp -O3
g++ -o $LDIR/src/buildDB  $LDIR/src/compress.cpp -O3
g++ -o $LDIR/src/uniqueDB  $LDIR/src/DTB_unique.cpp -O3
g++ -o $LDIR/src/translate  $LDIR/src/translate.cpp -O3
g++ -o $LDIR/src/classify  $LDIR/src/classify.cpp $LDIR/src/seqreader.cpp -O3

mv $LDIR/src/buildDB $LDIR/
mv $LDIR/src/uniqueDB $LDIR/
mv $LDIR/src/translate $LDIR/
mv $LDIR/src/classify $LDIR/

