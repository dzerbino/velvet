#!/usr/bin/perl -w
#
#       VelvetOptimiser.pl
#
#       Copyright 2008, 2009 Simon Gladman <simon.gladman@csiro.au>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

#
#   pragmas
#
use strict;

#
#   includes
#
use Getopt::Std;
use POSIX qw(strftime);
use FindBin;
use lib "$FindBin::Bin";
use VelvetOpt::Assembly;
use VelvetOpt::hwrap;
use VelvetOpt::gwrap;
use VelvetOpt::Utils;

#
#   global var decs
#
my %opts;
my @hashvals;
my %assemblies;
my $readfile;
my $logfile = "logfile.txt";
my $hashs = 19;
my $hashe = 31;
my $ass_num = 1;
my $interested = 0;

#
#
#   subroutines
#
#

#
#   isOdd
#
sub isOdd {
    my $x = shift;
    if($x % 2 == 1){
        return 1;
    }
    else {
        return 0;
    }
}

#
#   by_num - sorter
#
sub by_num { $a <=> $b }

#
#   getOptRoutine
#
sub getOptRoutine {

    my $readfile = shift;

    #   Choose the optimisation path depending on the types of read files sent to velveth
    #       For short only:                 shortOpt routine
    #       For short and long:             shortLong routine
    #       For short paired:               shortPaired routine
    #       For short and long paired:      longPaired routine
    #       For short paired and long:      shortPaired routine
    #       For short paired & long paired: shortlongPaired routine

    #look at velveth string ($readfile) and look for keywords from velvet manual...
    my $long = 0;
    my $longPaired = 0;
    my $shortPaired = 0;
    my $short = 0;

    #standard cases..
    if($readfile =~ /-short.? /) { $short = 1; }
    if($readfile =~ /-long /) { $long = 1; }
    if($readfile =~ /-shortPaired /) { $shortPaired = 1; }
    if($readfile =~ /-longPaired /) { $longPaired = 1; }

    #weird cases to cover the non-use of the short keyword (since its the default.)
    if(!($readfile =~ /(-short.? )|(-long )|(-shortPaired )|(-longPaired )/)) { $short = 1; } #if nothing is specified, assume short.
    if(!($readfile =~ /-short.? /) && ($readfile =~ /(-long )|(-longPaired )/)) { $short = 1; } #if long or longPaired is specified, also assum short since very unlikely to only have long...

    if($short && !($long || $longPaired || $shortPaired)){
        return "shortOpt";
    }
    elsif($short && $long && !($longPaired || $shortPaired)){
        return "shortLong";
    }
    elsif($short && $longPaired && !$shortPaired){
        return "longPaired";
    }
    elsif($short && $shortPaired && !$longPaired){
        return "shortPaired";
    }
    elsif($short && $shortPaired && $longPaired){
        return "shortLongPaired";
    }
    elsif($shortPaired && !$short && !$long && !$longPaired){
        return "shortPaired";
    }
    else {
        return "Unknown";
    }
}

#
#   covCutoff - the coverage cutoff optimisation routine.
#
sub covCutoff{

    my $ass = shift;
    #get the assembly score and set the current cutoff score.
    my $ass_score = $ass->{assmscore};
    print "In covCutOff and assembly score is: $ass_score..\n" if $interested;
	


	sub func {
		my $ass = shift;
		my $cutoff = shift;
		my $ass_score = $ass->{assmscore};
		my $ps = $ass->{pstringg};
        if($ps =~ /cov_cutoff/){
            $ps =~ s/cov_cutoff\s+\d+(\.\d+)?/cov_cutoff $cutoff/;
        }
        else {
            $ps .= " -cov_cutoff $cutoff";
        }
        $ass->{pstringg} = $ps;

        print STDERR strftime("%b %e %H:%M:%S", localtime), "\t\tSetting cov_cutoff to $cutoff.\n";
        print OUT strftime("%b %e %H:%M:%S", localtime), "\t\tSetting cov_cutoff to $cutoff.\n";

        my $worked = VelvetOpt::gwrap::objectVelvetg($ass);
        if($worked){
            $ass->getAssemblyDetails();
        }
        else {
            die "Velvet Error in covCutoff!\n";
        }
        $ass_score = $ass->{assmscore};
		print OUT $ass->toStringNoV();
		
		return $ass_score;
		
	}
	
	print STDERR strftime("%b %e %H:%M:%S", localtime), " Beginning coverage cutoff optimisation\n";
    print OUT strftime("%b %e %H:%M:%S", localtime), " Beginning coverage cutoff optimisation\n";

	my $dir = $ass->{ass_dir};
    $dir .= "/stats.txt";
    print "\tLooking for exp_cov in $dir\n";
    my $expCov = VelvetOpt::Utils::getExpCov($dir);
	
	my $a = 0;
	my $b = 0.8 * $expCov;
	my $t = 0.618;
	my $c = $a + $t * ($b - $a);
	my $d = $b + $t * ($a - $b);
	my $fc = func($ass, $c);
	my $fd = func($ass, $d);

	my $iters = 1;
	
	print STDERR "\t\tLooking for best cutoff score between $a and $b\n";
	print OUT "\t\tLooking for best cutoff score between $a and $b\n";
	
	while(abs($a -$b) > 1){
		
		if($fc > $fd){
			print STDERR "\t\tMax cutoff lies between $d & $b\n";
			print OUT "\t\tMax cutoff lies between $d & $b\n";
			$a = $d;
			$d = $c;
			$fd = $fc;
			$c = $a + $t * ($b - $a);
			$fc = func($ass, $c);
		}
		else {
			print STDERR "\t\tMax cutoff lies between $a & $c\n";
			print OUT "\t\tMax cutoff lies between $a & $c\n";
			$b = $c;
			$c = $d;
			$fc = $fd;
			$d = $b + $t * ($a - $b);
			$fd = func($ass, $d);
		}
		$iters ++;
	}

	print STDERR "\t\tOptimum value of cutoff is " . int($b) . "\n";
	print STDERR "\t\tTook $iters iterations\n";
	print OUT "\t\tOptimum value of cutoff is " . int($b) . "\n";
	print OUT "\t\tTook $iters iterations\n";

    return 1;

}

#
#   expCov - find the expected coverage for the assembly and run velvetg with that exp_cov.
#
sub expCov {

    print STDERR strftime("%b %e %H:%M:%S", localtime), " Looking for the expected coverage\n";
    print OUT strftime("%b %e %H:%M:%S", localtime), " Looking for the expected coverage\n";

    my $ass = shift;

    #need to get the directory of the assembly and add "stats.txt" to it and then send it to
    #the histogram methods in SlugsUtils.pm...
    my $dir = $ass->{ass_dir};
    $dir .= "/stats.txt";
    print "Looking for exp_cov in $dir\n";
    my $expCov = VelvetOpt::Utils::getExpCov($dir);

    print STDERR strftime("%b %e %H:%M:%S", localtime), "\t\tExpected coverage set to $expCov\n";
    print OUT strftime("%b %e %H:%M:%S", localtime), "\t\tExpected coverage set to $expCov\n";

    #re-write the pstringg with the new velvetg command..
    my $vg = $ass->{pstringg};
    if($vg =~ /exp_cov/){
        $vg =~ s/exp_cov\s+\d+/exp_cov $expCov/;
    }
    else {
        $vg .= " -long_mult_cutoff 1 -exp_cov $expCov";
    }

    $ass->{pstringg} = $vg;
    my $worked = VelvetOpt::gwrap::objectVelvetg($ass);
    if($worked){
        $ass->getAssemblyDetails();
    }
    else {
        die "Velvet Error in expCov!\n";
    }
    print OUT $ass->toStringNoV();

}

#
#   insLengthLong - get the Long insert length and use it in the assembly..
#
sub insLengthLong {
    print STDERR strftime("%b %e %H:%M:%S", localtime), " Getting the long insert length\n";
    print OUT strftime("%b %e %H:%M:%S", localtime), " Getting the long insert length\n";
    my $ass = shift;
    print STDERR "/tPlease type in the insert length for the long reads: ";
    my $len = <>;
    chomp($len);
    while($len =~ /\D+/){
        print STDERR "\tThe length needs to be a number, please re-enter: ";
        $len = <>;
        chomp($len);
    }
    print STDERR strftime("%b %e %H:%M:%S", localtime), " Running assembly with insert length $len\n";
    print OUT strftime("%b %e %H:%M:%S", localtime), " Running assembly with insert length $len\n";

    #re-write the pstringg with the new velvetg command..
    my $vg = $ass->{pstringg};
    if($vg =~ /ins_length_long/){
        $vg =~ s/ins_length_long\s+\d+/ins_length_long $len/;
    }
    else {
        $vg .= " -ins_length_long $len";
    }

    $ass->{pstringg} = $vg;
    my $worked = VelvetOpt::gwrap::objectVelvetg($ass);
    if($worked){
        $ass->getAssemblyDetails();
    }
    else {
        die "Velvet Error in insLengthLong!\n";
    }
    print OUT $ass->toStringNoV();
}

#
#   insLengthShort - get the short insert length and use it in the assembly..
#
sub insLengthShort {
    print STDERR strftime("%b %e %H:%M:%S", localtime), " Getting the short insert length\n";
    print OUT strftime("%b %e %H:%M:%S", localtime), " Getting the short insert length\n";
    my $ass = shift;
    print STDERR "\tPlease type in the insert length for the short reads: ";
    my $len = <>;
    chomp($len);
    while($len =~ /\D+/){
        print STDERR "\tThe length needs to be a number, please re-enter: ";
        $len = <>;
        chomp($len);
    }
    print STDERR strftime("%b %e %H:%M:%S", localtime), " Running assembly with short insert length $len\n";
    print OUT strftime("%b %e %H:%M:%S", localtime), " Running assembly with short insert length $len\n";

    #re-write the pstringg with the new velvetg command..
    my $vg = $ass->{pstringg};
    if($vg =~ /ins_length /){
        $vg =~ s/ins_length\s+\d+/ins_length $len/;
    }
    else {
        $vg .= " -ins_length $len";
    }

    $ass->{pstringg} = $vg;
    my $worked = VelvetOpt::gwrap::objectVelvetg($ass);
    if($worked){
        $ass->getAssemblyDetails();
    }
    else {
        die "Velvet Error in insLengthShort!\n";
    }
    print OUT $ass->toStringNoV();
}

#
#   usage stuff
#
my $usage = "\nVelvetOptimiser.pl: A script to run the Velvet assembler and optimise its output. Simon Gladman - CSIRO 2008, 2009.\n\n";
$usage .= "Usage: VelvetOptimiser.pl <-f 'velveth parameters'> [-s <hash start>] [-e <hash end>] [-a <yes>]\n\n";
$usage .= "Where:\t<-f 'velveth parameters'> is the parameter line normally passed to velveth in quotes.\n";
$usage .= "\t-s <hash start> The hash value you want velvet to start looking from. Default: 19. MUST BE ODD > 0 & <=31!\n";
$usage .= "\t-e <hash end> The hash value you want velvet to stop looking at. Default: 31. MUST BE ODD AND > START & <= 31!\n";
$usage .= "\t-a <yes> The final optimised assembly will include read tracking and amos file outputs (however, intermediate assemblies won't.)\n";
$usage .= "\nIf the optimizer requires an insert length for some paired end data, it will ask for it when it gets to the optimization step.\n";

#
#
#   get all the input parameters.
#
#
getopt('as:e:f:', \%opts);

#foreach my $key (keys %opts){
#    print "$key: " . $opts{$key} . "\n" if($opts{$key});
#}

print STDERR "
****************************************************

                 VelvetOptimiser.pl

         Simon Gladman - CSIRO 2008, 2009

****************************************************\n";

#
#
#check all input paramters
#
#
print STDERR strftime("%b %e %H:%M:%S", localtime), " Starting to check input parameters.\n";

unless($opts{'f'}){
    print STDERR "\tYou must supply the velveth parameter line.\n";
    die $usage;
}

$readfile = $opts{'f'};

if($opts{'s'}){
    $hashs = $opts{'s'};
    unless($hashs =~ /^\d+$/){ die "\tFatal error! Start hash not a number!\n$usage";}
    if($hashs > 31){
        print STDERR "\tStart hash value too high.  New start hash value is 31.\n";
        $hashs = 31;
    }
    if(!&isOdd($hashs)){
        $hashs = $hashs - 1;
        print STDERR "\tStart hash value not odd.  Subtracting one. New start hash value = $hashs\n";
    }

}

if($opts{'e'}){
    $hashe = $opts{'e'};
    unless($hashe =~ /^\d+$/){ die "\tFatal error! End hash not a number!\n$usage";}
    if($hashe > 31 || $hashe < 1){
        print STDERR "\tEnd hash value not in workable range.  New end hash value is 31.\n";
        $hashe = 19;
    }
    if($hashe < $hashs){
        print STDERR "\tEnd hash value lower than start hash value.  New end hash value = $hashs.\n";
        $hashe = $hashs;
    }
    if(!&isOdd($hashe)){
        $hashe = $hashe - 1;
        print STDERR "\tEnd hash value not odd.  Subtracting one. New end hash value = $hashe\n";
    }

}

#check the velveth parameter string..
my $vh_ok = VelvetOpt::hwrap::_checkVHString("check 21 $readfile");

unless($vh_ok){ die "Please re-start with a corrected velveth parameter string." }

print STDERR strftime("%b %e %H:%M:%S", localtime), " Finished checking input parameters.\n";

#
#
#   Perform common tasks - write details to log file and screen, run velveth and vanilla velvetg
#
#

#let user know about parameters to run with.
print STDERR "Will run velvet optimiser with the following paramters:\n";
print STDERR "\tVelveth parameter string:\n\t\t$readfile\n";
print STDERR "\tVelveth start hash values:\t$hashs\n";
print STDERR "\tVelveth end hash value:\t\t$hashe\n";
if($opts{'a'}){
    print STDERR "\tRead tracking for final assembly on.\n";
} else {
    print STDERR "\tRead tracking for final assembly off.\n";
}

#open the log file..
open OUT, ">$logfile" or die "Couldn't open $logfile for writing.";

print OUT strftime("%b %e %H:%M:%S", localtime), "\n";

#send run parameters to log file.
print OUT "Will run velvet optimiser with the following paramters:\n";
print OUT "\tVelveth parameter string:\n\t\t$readfile\n";
print OUT "\tVelveth start hash values:\t$hashs\n";
print OUT "\tVelveth end hash value:\t\t$hashe\n\n";
if($opts{'a'}){
    print OUT "\tRead tracking for final assembly on.\n";
} else {
    print OUT "\tRead tracking for final assembly off.\n";
}

#get the velveth and velvetg version numbers...
my $response = VelvetOpt::hwrap::_runVelveth(" ");
$response =~ /Version\s+(\d+\.\d+\.\d+)/s;
my $vhversion = $1;

$response = VelvetOpt::gwrap::_runVelvetg(" ");
$response =~ /Version\s+(\d+\.\d+\.\d+)/s;
my $vgversion = $1;

print OUT "Velveth version: $vhversion.\nVelvetg version: $vgversion\n";

#build the hashval array
for(my $i = $hashs; $i <= $hashe; $i += 2){
    push @hashvals, $i;
}

print STDERR strftime("%b %e %H:%M:%S", localtime), " Beginning velveth runs.\n";
print OUT strftime("%b %e %H:%M:%S", localtime), "\n\n\tBeginning velveth runs.\n";

#now make an assembly object for all of the hash values chosen, run velveth on them and getHashDetails...
foreach my $hashval (@hashvals){

    print STDERR strftime("%b %e %H:%M:%S", localtime), "\t\tRunning velveth with hash value: $hashval.\n";

    #make the velveth command line.
    my $vhline = "data_$hashval $hashval $readfile";

    #make a new VelvetAssembly and store it in the %assemblies hash...
    $assemblies{$ass_num} = VelvetOpt::Assembly->new(ass_id => $ass_num, pstringh => $vhline, versionh => $vhversion);

    #run velveth on this assembly object
    my $vhresponse = VelvetOpt::hwrap::objectVelveth($assemblies{$ass_num});

    unless($vhresponse){ die "Velveth didn't run on hash value of $hashval.\n$!\n";}

    #run the hashdetail generation routine.
    $vhresponse = $assemblies{$ass_num}->getHashingDetails();

    #print the objects to the log file...
    print OUT $assemblies{$ass_num}->toStringNoV();

    #increment the $ass_num
    $ass_num ++;
}

print STDERR strftime("%b %e %H:%M:%S", localtime), " Finished velveth runs.\n";

print STDERR strftime("%b %e %H:%M:%S", localtime), " Beginning vanilla velvetg runs.\n";
print OUT strftime("%b %e %H:%M:%S", localtime), "\n\n\tBeginning vanilla velvetg runs.\n";

#now for each assembly object, run vanilla velvetg and get their assembly details...
foreach my $key (sort by_num keys %assemblies){
    print STDERR strftime("%b %e %H:%M:%S", localtime), "\t\tRunning velvetg on assembly for hash value: " . $assemblies{$key}->{hashval} ."\n";

    #make the velvetg commandline.
    my $vgline = "data_" . $assemblies{$key}->{hashval};

    #save the velvetg commandline in the assembly.
    $assemblies{$key}->{pstringg} = $vgline;

    #run velvetg
    my $vgresponse = VelvetOpt::gwrap::objectVelvetg($assemblies{$key});

    unless($vgresponse){ die "Velvetg didn't run on the directory $vgline.\n$!\n";}

    #run the assembly details routine..
    $assemblies{$key}->getAssemblyDetails();

    #print the objects to the log file...
    print OUT $assemblies{$key}->toStringNoV();
}


#
#
#   Now perform a velvetg optimisation based upon the file types sent to velveth
#
#

#
#   get the best assembly so far...
#

my $bestId;
my $maxScore = -100;
my $asmscorenotneg = 1;

foreach my $key (keys %assemblies){
	if(($assemblies{$key}->{assmscore} != -1) && $asmscorenotneg){
    	if($assemblies{$key}->{assmscore} > $maxScore){
        	$bestId = $key;
        	$maxScore = $assemblies{$key}->{assmscore};
    	}
	}
	elsif($assemblies{$key}->{n50} && $asmscorenotneg){
		if($assemblies{$key}->{n50} > $maxScore){
			$bestId = $key;
			$maxScore = $assemblies{$key}->{n50};
		}
	}
	else {
		$asmscorenotneg = 0;
		if($assemblies{$key}->{totalbp} > $maxScore){
        	$bestId = $key;
        	$maxScore = $assemblies{$key}->{totalbp};
    	}
	}
}
print "\n\nThe best assembly so far is:\n" if $interested;
print $assemblies{$bestId}->toStringNoV() if $interested;

#   determine the optimisation route for the assembly based on the velveth parameter string.
my $optRoute = &getOptRoutine($readfile);

print STDERR strftime("%b %e %H:%M:%S", localtime), " Hash value of best assembly by assembly score: ". $assemblies{$bestId}->{hashval} . "\n";

print OUT strftime("%b %e %H:%M:%S", localtime), " Best assembly by assembly score - assembly id: $bestId\n";

print STDERR strftime("%b %e %H:%M:%S", localtime), " Optimisation routine chosen for best assembly: $optRoute\n";
print OUT strftime("%b %e %H:%M:%S", localtime), " Optimisation routine chosen for best assembly: $optRoute\n";

#now send the best assembly so far to the appropriate optimisation routine...

if($optRoute eq "shortOpt"){

    &covCutoff($assemblies{$bestId});

}
elsif($optRoute eq "shortLong"){

    &expCov($assemblies{$bestId});
    &covCutoff($assemblies{$bestId});

}
elsif($optRoute eq "longPaired"){
    &expCov($assemblies{$bestId});
    &insLengthLong($assemblies{$bestId});
    &covCutoff($assemblies{$bestId});
}
elsif($optRoute eq "shortPaired"){
    &expCov($assemblies{$bestId});
    &insLengthShort($assemblies{$bestId});
    &covCutoff($assemblies{$bestId});
}
elsif($optRoute eq "shortLongPaired"){
    &expCov($assemblies{$bestId});
    &insLengthShort($assemblies{$bestId});
    &insLengthLong($assemblies{$bestId});
    &covCutoff($assemblies{$bestId});
}
else{
    print STDERR "There was an error choosing an optimisation routine for this assembly.  Please change the velveth parameter string and try again.\n";
    print OUT "There was an error choosing an optimisation routine for this assembly.  Please change the velveth parameter string and try again.\n";
}

#   once it comes back from the optimisation routines, we need to turn on read tracking and amos output if it was selected in the options.
#
#
#   The final assembly run!
#
#
if($opts{'a'}){
    $assemblies{$bestId}->{pstringg} .= " -amos_file yes -read_trkg yes";

    my $final = VelvetOpt::gwrap::objectVelvetg($assemblies{$bestId});
    $assemblies{$bestId}->getAssemblyDetails();
}

print STDERR strftime("%b %e %H:%M:%S", localtime), "\n\n\nFinal optimised assembly details:\n";
print OUT strftime("%b %e %H:%M:%S", localtime), "\n\n\nFinal optimised assembly details:\n";
print STDERR $assemblies{$bestId}->toStringNoV();
print OUT $assemblies{$bestId}->toStringNoV();
print STDERR "\n\nAssembly output files are in the following directory:\n" . $assemblies{$bestId}->{ass_dir} . "\n\n";
print OUT "\n\nAssembly output files are in the following directory:\n" . $assemblies{$bestId}->{ass_dir} . "\n";

#delete superfluous directories..
foreach my $key(keys %assemblies){
	unless($key == $bestId){ 
		my $dir = $assemblies{$key}->{ass_dir};
		`rm -r $dir`;
	} 
}
