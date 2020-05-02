#!/usr/bin/env perl

# Copyright 2013-2019, Derrick Wood <dwood@cs.jhu.edu>
#
# This file is part of the Kraken 2 taxonomic sequence classification system.

# Modified by Bui Van-Kien (buikien.dp@sjtu.edu.cn) 
# used for CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Function: Print the header of each genome

use strict;
use warnings;
use File::Basename;
use Getopt::Long;

my $PROG = basename $0;

my $lenient = 0;
GetOptions("lenient" => \$lenient)
  or die "Usage: $PROG [--lenient] <fasta filename(s)>\n";

while (<>) {
  next unless /^>/;
  # while (/.../g) needed because non-redundant DBs sometimes have multiple
  #   sequence IDs in the header; extra sequence IDs are prefixed by
  #   '\x01' characters (if downloaded in FASTA format from NCBI FTP directly).
  while (/(?:^>|\x01)(\S+)/g) {
    my $seqid = $1;
     print "$seqid\n";
  }
}
