# CDKAM
Copyright 2019-2020\
Author: Van-Kien Bui, Chaochun Wei\
Department of Bioinformatics and Biostatistics, Shanghai Jiao Tong University

**1) Introduction**

Classification tool using Discriminative K-mers and Approximate Matching algorithm (CDKAM) is a new metagenome sequence classification tool for the third generation sequencing data with high error rates. 

**2) Requirements**

- Linux operation system
- Memory: 50 GB
- Disk space: 150 GB
- Perl 5.8.5 (or up) and GCC 6.3.1 (or up).
- Dustmasker https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/dustmasker/. 

It is suggested to install BLAST+ (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), which has already included dustmasker.
Low-complexity sequences, e.g. "ACACACACACACACACA", "ATATATATATATATATAT" are known to occur in many different organisms and are typically less informative for use in alignments. The masked regions are not processed further by CDKAM.

**3) Datasets** 

Datasets can be found at OneDrive: 
- First simulated dataset
https://1drv.ms/u/s!AvI5WFKEnJrGcNO2PkmiFgBz3lk?e=Ipgwac
- Second simulated dataset
https://1drv.ms/u/s!AvI5WFKEnJrGcv0BBP0jM72Vibs?e=j4tOX3
- A sample of sequencing Nanopore data
https://www.st-va.ncbi.nlm.nih.gov/bioproject/PRJNA493153

**4) Installation**

- First, download the package of the latest CDKAM release: https://github.com/SJTU-CGM/CDKAM

- Go in the extracted sub-directory "CDKAM". 
Then:\
$ ./install.sh

**5) Running CDKAM**
- Downloading database:\
./download --standard --db $DBname

- Building database:\
./build_database.sh $DBname

- Running classification:\
./CDKAM.sh $DBname input output --fasta/--fastq \
Using --fasta if the input is FASTA file, --fastq if the input is FASTQ file.

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
(Read ID) (Species/Genus taxonomy ID) (Genus taxonomy ID) (S/G) Scientific Name

Example:
- *1	-1	-1*
- *2	28116	816	 (S) Bacteroides*
- *3	590	590	 (G) Salmonella*

(S) indicates that the read is assigned to a taxonomy ID at Species level\
(G) indicates that the read is assigned to a taxonomy ID at Genus level



**Contact:**

If you have any questions, feel free to contact us:\
   buikien.dp@sjtu.edu.cn\
   ccwei@sjtu.edu.cn
   
   
