NAME
====

AssemblyAssembler.py

VERSION
=======

Version 1.1

LICENCE
=======

Copyright 2010 - Jacob Crawford - Cornell University.
	
jc598@cornell.edu

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
    
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
        
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.


INTRODUCTION
============

The AssemblyAssembler (AA) is designed to automate a directed series of 
assemblies using the Velvet assembler (Daniel Zerbino, EBI UK).  The 
assembly routine was initially developed in an attempt to improve de novo 
transcriptome assemblies using Velvet (prior to the release of Oases 
(Daniel Zerbino, EBI UK)).  The AA conducts Velvet assemblies using 
default parameter values across a user specified range of kmer values.  
It then determines which kmer values produced the best assembly (largest 
max contig) and conducts additional assemblies in that region of the 
kmer-parameter space.  Lastly, the AA takes contigs from all previous 
assemblies and uses them as input for a final assembly.  The final contig 
set is stored in Finalcontigs.fa and all Velvet operations recorded in 
GrandVelvetLog.txt, both stored in the directory FinalDir.  

DISCLAIMER
==========
This script comes with NO GUARANTEE. Misassemblies are possible in Velvet
assemblies and are possible using this script.  It is up to the user to 
determine the quality of the assembly.  

PREREQUISITES
=============

Python >= 2.5
Python modules sys and os

COMMAND LINE
============
	
	AssemblyAssembler.py -s INT -e INT -v /path/to/velvetdir -f 'velveth input line'
  
  	-s		The starting (lowest) hash value (no default). Enter integer value (no quotes).
  			It is not recommended that you attempt small values (< k=15) for this parameter.
  	
  	-e		The end (highest) hash value (no default). Enter integer value (no quotes).  
  			The difference between the highest and lowest hash values must be at least 16.
  			I found it useful to set this parameter to (length of read-1).
  	
  	-v		The full path to the Velvet directory containing velveth and velvetg binaries. 
  			Example: /programs/velvet_0.7.61      (no quotes).
  			
  	-f		The velveth parameters you would enter for a standard velveth run.  This entry MUST
  			be INSIDE QUOTES. You MUST enter 1) the input file type (e.g. -fasta or -fastq), 2) 
  			the read type (e.g. -short or =shortPaired) and 3) the full path to the sequence 
  			file(s). Example: Ô-fastq -short /myfiles/Illumina/s_8_1_sequence.txtÕ).
  			
  	*** At this time, this script is restricted to very basic velveth and velvetg runs.  In other
  		words, advanced parameters changes, and other tweaks such as, for example, 
  		-cov_cutoff and -min_contig_lgth are not supported.  This is simply due to the fact that 
  		such parameter adjustments do not seem to benefit de novo transcriptome assemblies.  
  		There is no restriction on the number of categories and the MAXKMERLENGTH as they are simply
  		entered in the velveth input line (-f), provided you have compiled Velvet accordingly.

EXAMPLES
========

Assemble a lane of Illumina single-end reads by searching across the range k = 19 to k = 35:

% AssemblyAssembler.py -s 19 -e 35 -v /programs/velvet_0.7.61 -f '-short -fastq s_8_1_sequence.txt'

Assemble two lanes of Illumina paired-end reads by searching across the range k = 19 to k = 59:

% AssemblyAssembler.py -s 19 -e 35 -v /programs/velvet_0.7.61 -f Ô-fastq -shortPaired s_7_shuffled.txt -shortPaired2 s_8_shuffled.txtÕ

BUGS
====

* None that I am aware of.

CONTACT
=======

Jacob Crawford <jc598@cornell.edu>




