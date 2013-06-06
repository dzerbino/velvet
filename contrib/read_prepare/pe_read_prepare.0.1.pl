#---------------------------------------[SCRIPT : PAIRED-END READ PREPARE]----------------------------------------

# Script        : pe_read_prepare.0.1.pl
# Date          : June 03, 2013
# Usage         : perl ~/pe_read_prepare.0.1.pl [options] <READ1> <READ2> <Output_Prefix>
# Author	: Manojkumar Sumathiselvaraju
# Contact	: manojsbiotech [at] gmail [dot] com
# Institute	: Next Generation Genomics Facility, C-CAMP, NCBS-TIFR, Bangalore, India

#---------------------------------------[ASSIGNING & VERIFYING PARAMETERS]----------------------------------------

use Getopt::Long;

if ($#ARGV < 0){
	how_to_use();
}

my ($s, $s_r1, $s_r2, $l, $l_r1, $l_r2, $Q, $q, $p, $n);

GetOptions ('s=i'	=> \$s,
	    's_r1=i'	=> \$s_r1,
	    's_r2=i'	=> \$s_r2,
	    'l=i'	=> \$l,
	    'l_r1=i'	=> \$l_r1,
	    'l_r2=i'	=> \$l_r2,
	    'q=i'	=> \$q,
	    'e=i'       => \$e,
	    'p=i'	=> \$p,
	    'n=i'	=> \$n);

$e = 33	 if (!defined $e);
$q = 20	 if (!defined $q);
$p = 100 if (!defined $p);
$n = 0   if (!defined $n);

if(defined $s && (defined $s_r1 || defined $s_r2 )){
	print STDERR "\nError: Wrong Parameter Combination\n";
	print STDERR "Parameters -s and (s_r1 or s_r2) cannot be used together\n";
	how_to_use();
}

if(defined $l && (defined $l_r1 || defined $l_r2 )){
	print STDERR "\nError: Wrong Parameter Combination\n";
	print STDERR "Parameters -l and (l_r1 or l_r2) cannot be used together\n";
        how_to_use();
}	

if(!defined $s || (!defined $s_r1 && defined $s_r2)){
	$s = 1;
}

if(defined $s){
	$s_r1 = $s, $s_r2 = $s;
}

if(defined $l){
	$l_r1 = $l-$s+1, $l_r2 = $l-$s+1;
}
elsif(defined $l_r1 && defined $l_r2){
	$l_r1 = $l_r1-$s+1, $l_r2 = $l_r2-$s+1;
}

#--------------------------------------------[OPENING INPUT FASTQ FILES]------------------------------------------

$T0=localtime();
print "$T0\tPaired-End Read Preprocessing Has Begum\n";

$gzip_R1 = grep(/\.gz/, $ARGV[0]);
$gzip_R2 = grep(/\.gz/, $ARGV[1]);

if(($gzip_R1 == 1) && ($gzip_R2 == 1)){
	open("R1", "gunzip -c $ARGV[0] |");
        open("R2", "gunzip -c $ARGV[1] |");
}
elsif(($gzip_R1 == 1) && ($gzip_R2 == 0)){
	open("R1", "gunzip -c $ARGV[0] |");
	open("R2", "$ARGV[1]");
}
elsif(($gzip_R1 == 0) && ($gzip_R2 == 1)){
	open("R1", "$ARGV[0]");
	open("R2", "gunzip -c $ARGV[1] |");
}
else{	
	open("R1", "$ARGV[0]");
	open("R2", "$ARGV[1]");
}

#---------------------------------------[READING INPUT & OPENIG OUTPUT FILES]-------------------------------------

open(O1A, "| gzip -c > $ARGV[2]_01_PE_R1.fastq.gz");
open(O1B, "| gzip -c > $ARGV[2]_01_PE_R2.fastq.gz");

open(O2A, "| gzip -c > $ARGV[2]_02_SR_R1.fastq.gz");
open(O2B, "| gzip -c > $ARGV[2]_02_SR_R2.fastq.gz");

open(O3A, "| gzip -c > $ARGV[2]_03_PAIRS.fasta.gz");
open(O3B, "| gzip -c > $ARGV[2]_03_SR_R1.fasta.gz");
open(O3C, "| gzip -c > $ARGV[2]_03_SR_R2.fasta.gz");

if(!defined $l && !defined $l_r1 &&  !defined $l_r2){
	while($_=<R1>){
		$p1_H1=$_;   $p1_S=<R1>; $p1_H2=<R1>; $p1_Q=<R1>;
        	$p2_H1=<R2>; $p2_S=<R2>; $p2_H2=<R2>; $p2_Q=<R2>;
		
		chomp ($p1_H1, $p1_S, $p1_Q, $p2_H1, $p2_S, $p2_Q);

		$p1_S = substr($p1_S, $s_r1-1,); $p1_Q = substr($p1_Q, $s_r1-1,);
                $p2_S = substr($p2_S, $s_r2-1,); $p2_Q = substr($p2_Q, $s_r2-1,);
 
		@all_counts = read_filter ($q, $p, $n, $e, $p1_H1, $p1_S, $p1_Q, $p2_H1, $p2_S, $p2_Q, @all_counts);
	}
}
else{
	while($_=<R1>){
		$p1_H1=$_;   $p1_S=<R1>; $p1_H2=<R1>; $p1_Q=<R1>;
		$p2_H1=<R2>; $p2_S=<R2>; $p2_H2=<R2>; $p2_Q=<R2>;

		chomp ($p1_H1, $p1_S, $p1_Q, $p2_H1, $p2_S, $p2_Q);

		$p1_S = substr($p1_S, $s_r1-1, $l_r1); $p1_Q = substr($p1_Q, $s_r1-1, $l_r1);
		$p2_S = substr($p2_S, $s_r2-1, $l_r2); $p2_Q = substr($p2_Q, $s_r2-1, $l_r2);
	
		@all_counts = read_filter ($q, $p, $n, $e, $p1_H1, $p1_S, $p1_Q, $p2_H1, $p2_S, $p2_Q, @all_counts);
	}
}
close(R1, R2, O1A, O1B, O2A, O2B, O3A, O3B, O3C);

print_summary_stats(@all_counts);

#---------------------------------------------[SUBROUTINE : READ_FILTER]------------------------------------------

sub read_filter() {
	($a_q, $b_p, $c_n, $Q_a, $r1_h, $r1_s, $r1_q, $r2_h, $r2_s, $r2_q, @tempa) = @_;
	($t_pair, $t_se_r1, $t_se_r2, $t_qf_r1, $t_qf_r2, $t_nn_r1, $t_nn_r2) = @tempa;

	if($c_n > 0){
		$n_r1 = count_ns($r1_s);	
		$n_r2 = count_ns($r2_s);
	}
	else{
		$n_r1 = grep (/N/, $r1_s);
		$n_r2 = grep (/N/, $r2_s);
	}
	if(($n_r1 <= $c_n) && ($n_r2 <= $c_n)){
		$q_r1 = quality_filter ($a_q, $r1_q, $Q_a, $n_r1);
		$q_r2 = quality_filter ($a_q, $r2_q, $Q_a, $n_r2);
		if(($q_r1 >= $b_p) && ($q_r2 >= $b_p)){
			print O1A "$r1_h\n$r1_s\n+\n$r1_q\n";
			print O1B "$r2_h\n$r2_s\n+\n$r2_q\n";
			$r1_h = substr($r1_h, 1,);
			$r2_h = substr($r2_h, 1,);
			print O3A ">$r1_h\n$r1_s\n";
                        print O3A ">$r2_h\n$r2_s\n";
			$t_pair += 1;
		}
		elsif(($q_r1 >= $b_p) && ($q_r2 <= $b_p)){
			$tr1_h = substr($r1_h, 1,);
			print O2A "$r1_h\n$r1_s\n+\n$r1_q\n";
			print O3B ">$tr1_h\n$r1_s\n";
			$t_se_r1 += 1; $t_qf_r2 += 1;
		}
		elsif(($q_r1 <= $b_p) && ($q_r2 >= $b_p)){
			$tr2_h = substr($r2_h, 1,);
			print O2B "$r2_h\n$r2_s\n+\n$r2_q\n";
			print O3C ">$tr2_h\n$r2_s\n";
			$t_se_r2 += 1; $t_qf_r1 += 1;
                }
		else{
			$t_qf_r2 += 1; $t_qf_r1 += 1;
		}
	}
	elsif (($n_r1 <= $c_n) && ($n_r2 >= $c_n)){
		$t_nn_r2 += 1;
		$q_r1 = quality_filter ($a_q, $r1_q, $Q_a, $n_r1);
		if($q_r1 >= $b_p){
			$tr1_h = substr($r1_h, 1,);
			print O2A "$r1_h\n$r1_s\n+\n$r1_q\n";
			print O3B ">$tr1_h\n$r1_s\n";
			$t_se_r1 += 1;
		}
		else{
			$t_qf_r1 += 1;
		}
	}
	elsif (($n_r1 >= $c_n) && ($n_r2 <= $c_n)){
		$t_nn_r1 += 1;
		$q_r2 = quality_filter ($a_q, $r2_q, $Q_a, $n_r2);
		if($q_r2 >= $b_p){
			$tr2_h = substr($r2_h, 1,);
			print O2B "$r2_h\n$r2_s\n+\n$r2_q\n";
			print O3C ">$tr2_h\n$r2_s\n";
			$t_se_r2 += 1;
                }
                else{
			$t_qf_r2 += 1;
                }
	}
	else{
		$t_nn_r1 += 1; $t_nn_r2 += 1;
	}
	return ($t_pair, $t_se_r1, $t_se_r2, $t_qf_r1, $t_qf_r2, $t_nn_r1, $t_nn_r2);
}

#-----------------------------------------[SUBROUTINE : COUNT N'S PER READ]---------------------------------------

sub count_ns() {
	($nseq) = @_; $ncount = 0;
	for($i=0; $i < length $nseq; $i += 1){
		$ncount += 1 if(substr($nseq, $i, 1) =~ m/N/);
	}
	return($ncount);
}

#-----------------------------------------[SUBROUTINE : READ QUALITY FILTER]---------------------------------------

sub quality_filter() {
	($mq, $rq, $Q_b, $n_b) = @_; $high_qual = 0;
	for($i=0; $i < length $rq; $i+=1){		
		$high_qual += 1 if((ord(substr($rq,$i,1))-$Q_b) >= $mq);			
	}
	$qpercent = ($high_qual/((length $rq)-$n_b))*100;	
	return($qpercent);	
}

#-----------------------------------------[SUBROUTINE : PRINT SUMMARY STATS]---------------------------------------

sub print_summary_stats(){
	@a = @_;
	$total_rp = $a[0]*2;
	$total_se = $a[1] + $a[2];
	$total_qf = $a[3] + $a[4];
	$total_nn = $a[5] + $a[6];
	$total_reads += $_ for @a;
	$total_reads += $a[0]; $thalf = $total_reads/2;

	print "\nInput Read Stats:\n";
	print "Total No. Of Reads		= $total_reads\n";
	print "Total No. Of R1 Reads		= $thalf\n";
	print "Total No. Of R2 Reads		= $thalf\n\n";
	
	print "Output Read Stats:\n";
	print "Total No. Of Paired Reads	= $total_rp ($a[0] + $a[0])\n";
	print "Total No. Of Singleton Reads	= $total_se ($a[1] + $a[2])\n";
	print "Total No. Of NNNN Reads		= $total_nn ($a[5] + $a[6])\n";
	print "Total No. Of Quality Fail Reads	= $total_qf ($a[3] + $a[4])\n\n";
	
	$T=localtime();
	print "$T\tProcessing Completed. Reads Are Ready For Assembly.\n";
}

#----------------------------------------------[SUBROUTINE : HOW TO USE]------------------------------------------

sub how_to_use(){
	print STDERR <<EOF;

Usage:

perl pe_read_prepare.0.1.pl [Options] <READ1> <READ2> <Output_Prefix>

Optional Parameters:

-s      <integer>       <read pair trim start position>		default (1)
-s_r1   <integer>       <read 1 trim start position>		default (1)
-s_r2   <integer>       <read 2 trim start position>		default (1)

-l      <integer>       <read pair trim length>			default (no trimming)
-l_r1   <integer>       <read 1 trim length>         		default (no trimming)
-l_r2   <integer>       <read 2 trim length>			default (no trimming)

-e	<33 / 64>	<quality encoding>			default (33)
-q      <integer>       <minimum base quality score>		default (20)
-p      <integer>       <% of base with Qscore >= q>		default (100)
-n      <integer>       <No. of N's to tolerate per read>	default (0)

List Of Output Files:

1) <Output_File_Name>_01_PE_R1.fastq.gz
2) <Output_File_Name>_01_PE_R2.fastq.gz
3) <Output_File_Name>_02_SR_R1.fastq.gz
4) <Output_File_Name>_02_SR_R2.fastq.gz
5) <Output_File_Name>_03_PAIRS.fasta.gz
6) <Output_File_Name>_03_SR_R1.fasta.gz
7) <Output_File_Name>_03_SR_R2.fasta.gz

Note:

1) Either -s Or (-s_r1 and -s_r2), Either -l Or (-l_r1 and -l_r2)
2) The input files must be in fastq format & they must be paired
3) Gzipped files are also accepted as inputs
4) Values of -l Or (-l_r1 and -l_r2) must be <= Read Length

EOF
        exit(1);
}

#------------------------------------------------------------------------------------------
