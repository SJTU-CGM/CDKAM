#!/bin/bash
#
# CDKAM: a metagenomic classification tool using discriminative k-mers and approximate matching strategy
# Copyright 2019-2020
# Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University
# Contact information: buikien.dp@sjtu.edu.cn, ccwei@sjtu.edu.cn
#

set -u  # Protect against uninitialized vars.
set -e  # Stop on error
set -o pipefail  # Stop on failures in non-final pipeline commands

masking_flag=""
if [ -z "$CDKAM_MASK_LC" ]; then
  masking_flag="--no-mask"
fi

ftp_flag=""
if [ -n "$CDKAM_USE_FTP" ]; then
  ftp_flag="--use-ftp"
fi


if [ -f $CDKAM_DB_NAME/targets.txt ]; then
	rm -f $CDKAM_DB_NAME/targets.txt
fi

if [ ! -d $CDKAM_DB_NAME ]; then
	mkdir $CDKAM_DB_NAME
fi

touch $CDKAM_DB_NAME/targets.txt

./download_taxonomy.sh $CDKAM_DB_NAME

download --db $CDKAM_DB_NAME --download-library archaea $masking_flag $ftp_flag
download --db $CDKAM_DB_NAME --download-library bacteria $masking_flag $ftp_flag
download --db $CDKAM_DB_NAME --download-library fungi $masking_flag $ftp_flag
download --db $CDKAM_DB_NAME --download-library viral $masking_flag $ftp_flag
#download --db $CDKAM_DB_NAME --download-library human --no-mask $ftp_flag

