# CDKAM
Copyright 2019-2020\
Author: Van-Kien Bui, Chaochun Wei\
Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University

**1) Introduction**

Classification tool using Discriminative K-mers and Approximate Matching algorithm (CDKAM) is a new metagenome sequence classification tool for the third generation sequencing data with high error rates. 

**2) Requirements**

- Linux operation system
- Memory: 70 GB
- Disk space: 200 GB
- Perl 5.8.5 (or up) and GCC 4.8.5 (or up).
- Dustmasker https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/. 

It is suggested to install BLAST+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), which has already included dustmasker.
Low-complexity sequences, e.g. "ACACACACACACACACA", "ATATATATATATATATAT" are known to occur in many different organisms and are typically less informative for use in alignments. The masked regions are not processed further by CDKAM.


**3) Datasets** 

Datasets can be found at OneDrive: 
- The first simulated dataset
https://1drv.ms/u/s!AvI5WFKEnJrGeQlkB-KTexns4m8?e=lxVWKy
- The second simulated dataset
https://1drv.ms/u/s!AvI5WFKEnJrGeu0OzwT1556LlG0?e=ThGos7
- The third simulated dataset
https://1drv.ms/u/s!AvI5WFKEnJrGe6q7R76aKHVx29k?e=bg2ogf
- The fourth simulated dataset
https://1drv.ms/u/s!AvI5WFKEnJrGfDYHCWoOfBN06zs?e=g3m7Zz
- A sample of sequencing Nanopore MinION data
https://www.st-va.ncbi.nlm.nih.gov/bioproject/PRJNA493153
- Zymo mock dataset:
https://github.com/LomanLab/mockcommunity

**4) Installation**

- Firstly, download the package of the latest CDKAM release: https://github.com/SJTU-CGM/CDKAM

- Then, go in the extracted sub-directory "CDKAM". 
Then:\
$ ./install.sh

**5) Running CDKAM**

It might take 6-8 hours for downloading the reference genomes (about 85 GB), and 24-30 hours for building the database.\
To build CDKAM with different value of X, please change the value of variable RANGE in /src/compress.cpp file as follows:\
RANGE = 5 for X = 20%\
RANGE = 7 for X = 15%\
RANGE = 10 for X = 10%\
RANGE = 20 for X = 5%\
The default version of CDKAM selects X = 15% with RANGE = 7.


- Downloading database:\
*Standard installation with archaea, bacteria and viral reference genomes*\
./download --standard --db $DBname\
*Or custom installation*\
mkdir $DBname\
./download_taxonomy.sh $DBname\
./download --download-library archaea --db $DBname\
./download --download-library bacteria --db $DBname\
./download --download-library fungi --db $DBname\
./download --download-library viral --db $DBname\
./download --download-library human --db $DBname

- Building database:\
./build_database.sh $DBname

- Running classification by default (using approximate matching strategies):\
./CDKAM.sh $DBname input output --fasta/--fastq \
Using --fasta if the input is FASTA file, --fastq if the input is FASTQ file.

- Multi-threading:\
./CDKAM.sh $DBname input output --fasta/--fastq nthread N\
where N is the number of threads.

- CDKAM also supports classification in Exact Matching mode:\
./CDKAM_EM.sh $DBname input output --fasta/--fastq 

- Running translation:\
./translate $DBname input output\
,where input is the result of the previous classification process.

**6) Output format**

Normal mode\
(Read ID) (Length of read) (Taxonomy ID)

Example:
- *1	985	-1*
- *2	733	28116*
- *3	886	590*

Translation mode\
(Read ID) (Genus taxonomy ID) (Genus taxonomy ID)  (AG) Full taxomomic path to Genus level      
or \
(Read ID) (Genus taxonomy ID) (Species taxonomy ID)  (AS) Full taxomomic path to Species level     

Example:
- *1	-1	-1*
- *2	816	28116	 (AS)	(P) Bacteroidetes | (C) Bacteroidia | (O) Bacteroidales | (F) Bacteroidaceae | (G) Bacteroides | (S) Bacteroides ovatus*
- *3	590	590	 (AG)	(P) Proteobacteria | (C) Gammaproteobacteria | (O) Enterobacterales | (F) Enterobacteriaceae | (G) Salmonella*

(AS) indicates that the read is assigned to a taxonomy ID at Species level\
(AG) indicates that the read is assigned to a taxonomy ID at Genus level

**7) Testing tool**

Compare the classification results between CDKAM and Kraken2. The database contains archaea, viral and fungi.

https://github.com/buikiendp/TestCDKAM



**Contact:**

If you have any questions, feel free to contact us:\
   buikien.dp@sjtu.edu.cn\
   ccwei@sjtu.edu.cn
   
   
