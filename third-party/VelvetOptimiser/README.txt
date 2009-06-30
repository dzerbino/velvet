Usage: VelvetOptimiser.pl <-f 'velveth parameters'> [-s <hash start>] [-e <hash end>] [-a yes]

Where:  <-f 'velveth parameters'> is the file section of the parameter line normally passed to velveth in quotes.
        -s <hash start> The hash value you want velvet to start looking from. Default: 19. MUST BE ODD > 0 & <=31!
        -e <hash end> The hash value you want velvet to stop looking at. Default: 31. MUST BE ODD AND > START & <= 31!
        -a <yes> The final optimised assembly will include read tracking and amos file outputs (however, intermediate assemblies won't.)

Examples:

VelvetOptimiser.pl -f '-short -fastq test_reads.fastq -longPaired test_long.fa' -s 23 -e 27 -a yes

This example will run velvet using illumina fastq reads in test_reads.fastq and some paired end long reads in test_long.fa.
It will run using hash values 23, 25 and 27 and will use read tracking on the final optimised assembly.

VelvetOptimiser.pl -f '-shortPaired short_paired.fa -short -fastq short.fastq' -s 17 -e 21

This example will run velvet using Illumina fasta short paired reads and another set of single ended fastq reads.  
It will use k values of 17, 19 and 21 and then optimise what it thinks is the best one.

The script looks at the velveth parameters and decides what to do to get the best assembly.  It will ask for an insert length about halfway through its run if there are any specified paired end reads.
It uses a better search method for the optimisation of coverage cutoff and has a max that is 0.8 * exp_cov.  

It automatically uses -long_mult_cutoff 1. This can be altered only by editing the script at this stage.  

Note that the velveth parameters have to be in ' marks.

