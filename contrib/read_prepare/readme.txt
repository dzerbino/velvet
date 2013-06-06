
------------------------------------------------------------------------------------------

PE_READ_PREPARE.pl: Preprocessing of paired-end data for de novo genome assembly
<Manojkumar Sumathiselvaraju>

Introduction:
   Quality filtering of raw data forms the primary step in NGS data analysis.
This script combines the features of popular quality filtering tools like
FASTX-toolkits. Commonly, reads from a paired-end data are pre-processed
separately, and paired later for downstream analysis. This kind of anlaysis has 
the following pros and cons.

Pros: 1) One can process read 1 and read 2 seperately at the same time.
Cons: 1) Read pairing will become a memory intensive process.
      2) Many times one has to re-write the script depending on header format for pairing. 

Here, this script attempts to circumvent the above mentioned two issues by
processing both the read pairs at one go. By doing so, the header information is
merely used for printing in the quality passed reads, and also allows one to
process data of any volume even in a laptop provided there is enough space to
store the output.

Disclaimer:
   Though I've tested this script several time, there is no guarantee that the
results will be perfect. Hence, it's up to the user to optimize the paramters
depending on his/her data quality profile.

Any Bugs:
   If you find any bugs in this script or need more information, please
reach me at manojsbiotech [at] gmail [dot] com.

------------------------------------------------------------------------------------------

How To Use:
	To know how to use the script, in your terminal just run the script as
described below:

[user@ngs-server]$ perl ~/pe_read_prepare.0.1.pl 

Usage:

perl pe_read_prepare.0.01.pl [Options] <READ1> <READ2> <Output_Prefix>

Optional Parameters:

-s      <integer>       <read pair trim start position>         default (1)
-s_r1   <integer>       <read 1 trim start position>            default (1)
-s_r2   <integer>       <read 2 trim start position>            default (1)

-l      <integer>       <read pair trim length>                 default (no trimming)
-l_r1   <integer>       <read 1 trim length>                    default (no trimming)
-l_r2   <integer>       <read 2 trim length>                    default (no trimming)

-e      <33 / 64>       <quality encoding>                      default (33)
-q      <integer>       <minimum base quality score>            default (20)
-p      <integer>       <% of base with Qscore >= q>            default (100)
-n      <integer>       <No. of N's to tolerate per read>       default (0)

List Of Output Files:

1) <Output_File_Name>_01_PE_R1.fastq.gz
2) <Output_File_Name>_01_PE_R2.fastq.gz

[Files (1) & (2) will contain read pairs that passed your quality filter. The
order of the read pairs are conserved.]

3) <Output_File_Name>_02_SR_R1.fastq.gz
4) <Output_File_Name>_02_SR_R2.fastq.gz

[Files (3) & (4) will contain quality filter passed reads whos mate got
discarded during the process.]

5) <Output_File_Name>_03_PAIRS.fasta.gz

[File (5) will contain quality filter passed shuffled read pairs in fasta 
format ready for Velvet assembly. It can be used for velvet assembly with
-shortPaired parameter.]

6) <Output_File_Name>_03_SR_R1.fasta.gz
7) <Output_File_Name>_03_SR_R2.fasta.gz

[Files (6) & (7) will contain quality filter passed reads whos mate got
discarded during the process in fasta format. These files can be used for
velvet assembly with -short2 and -short3 parameters respectively.]

Note:

1) Use either -s Or (-s_r1 and -s_r2); either -l Or (-l_r1 and -l_r2)
2) The input files must be in fastq format & they must be paired
3) Gzipped fastq files are also accepted as inputs
4) Values of -l Or (-l_r1 and -l_r2) must be <= Read Length

------------------------------------------------------------------------------------------

EXAMPLE PE_READ_PREPARE: 

[user@ngs-server]$ sh ./run_example_read_prepare.sh

perl pe_read_prepare.0.1.pl -e 64 -q 10 -p 100 -s 1 -n 2 -l_r1 15 -l_r2 5 \ 
			    read1.fq.gz read2.fq.gz READ >& log.read_example.txt &

INPUT FILES

1) read1.fq.gz
2) read2.fq.gz

OUTPUT FILES

1) READ_01_PE_R1.fastq.gz
2) READ_01_PE_R2.fastq.gz

3) READ_02_SR_R1.fastq.gz
4) READ_02_SR_R2.fastq.gz

5) READ_03_PAIRS.fasta.gz
6) READ_03_SR_R1.fasta.gz
7) READ_03_SR_R2.fasta.gz

------------------------------------------------------------------------------------------
