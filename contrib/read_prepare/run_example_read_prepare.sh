
perl pe_read_prepare.0.1.pl -e 64 -q 10 -p 100 -n 2 -s 1 -l_r1 15 -l_r2 5 read1.fq.gz read2.fq.gz READ >& log.read_example.txt &

#-------------------------------------------------------------------------------------------------------------------------------    
# Parameter Explanation
#-------------------------------------------------------------------------------------------------------------------------------    
#
#	Mandatory Inputs
#	
#	read1.fq.gz	=>	Input read - first in pair (R1)
#	read2.fq.gz	=>	Input read - second in pair (R2)
#	READ		=>	Output prefix
#
#	Optinal Parameters
#
#	-e		=>	Quality encoding (+64)
#	-q		=>	Minimum base quality for filtering (10)
#	-p		=>	Percentage of bases with minimum base quality (100)
#	-n		=>	Number of N's to tolerate in the preprocessed read (2)
#	-s		=>	Start position for both read1 and read2 (1)
#	-l_r1		=>	Length of read1 (15)
#	-l_r2		=>	Length of read2 (5)
#
#	Tip
#
#	The parameters for quality filtering were determined after infering the data quality. Here, FastQC was used to decide on 
#	the parameters. You can run FastQC on the read1.fq.gz and read2.fq.gz to understand the above statements.
#
#	[user@ngs-server]$ ~/FastQC/./fastqc -k 7 --nogroup -t 2 read1.fq.gz read2.fq.gz & 
#
#--------------------------------------------------------------------------------------------------------------------------------	
