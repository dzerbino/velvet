#       VelvetOpt::Assembly.pm
#
#       Copyright 2008 Simon Gladman <simon.gladman@csiro.au>
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
package VelvetOpt::Assembly;

=head1 NAME

VelvetOpt::Assembly.pm - Velvet assembly container class.

=head1 AUTHOR

Simon Gladman, CSIRO, 2007, 2008.

=head1 LICENSE

Copyright 2008 Simon Gladman <simon.gladman@csiro.au>

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

=head1 SYNOPSIS

    use VelvetOpt::Assembly;
    my $object = VelvetOpt::Assembly->new(
        timestamph => "23 November 2008 15:00:00",
        ass_id => "1",
        versionh => "0.7.04",
        ass_dir => "/home/gla048/Desktop/newVelvetOptimiser/data_1"
    );
    print $object->toString();

=head1 DESCRIPTION

A container class to hold the results of a Velvet assembly.  Includes timestamps,
version information, parameter strings and assembly output metrics.

=head2 Uses

=over 8

=item strict

=item warnings

=item Carp

=head2 Fields

=over 8

=item assmscore

The assembly score metric for this object

=item timstamph

The timestamp of the start of the velveth run for this assembly

=item timestampg

The date and time of the end of the velvetg run.

=item ass_id

The assembly id number.  Sequential for all the runs for this optimisation.

=item versionh

The version number of velveth used in this assembly

=item versiong

The version number of velvetg used in this assembly

=item readfilename

The name of the file containing all the reads (or a qw of them if more than one...)

=item pstringh

The velveth parameter string used in this assembly

=item pstringg

The velvetg parameter string used in this assembly

=item ass_dir

The assembly directory path (full)

=item hashval

The hash value used for this assembly

=item rmapfs

The roadmap file size

=item sequences

The total number of sequences in the input files

=item nconts

The number of contigs in the final assembly

=item totalbp

The total number of bases in the contigs

=item n50

The n50 of the assembly

=item maxlength

The length of the longest contig in the assembly

=item maxcont

The size of the largest contig in the assembly

=item nconts1k

The number of contigs greater than 1k in size

=item totalbp1k

the sum of the length of contigs > 1k in size

=item velvethout

The velveth output

=item velvetgout

The velvetg output

=back
=back

=head2 Methods

=over 8

=item new

Returns a new VelvetAssembly object.

=item accessor methods

Accessor methods for all fields.

=item _checkVHString

Checks the velveth parameter string for completeness..  Returns 1 if correct or 0 if not

=item calcAssemblyScore

Calculates the assembly score of the object (after velvetg has been run.) and stores it in self.

=item getHashingDetails

Gets the details of the outputs from the velveth run and stores it in self.

=item getAssemblyDetails

Gets the details of the outputs from the velvetg run and stores it in self.

=item toString

Returns a string representation of the object's contents.

=item toStringNoV

Returns a string representation of the object's contents without the velvet outputs which are large.

=back

=cut

use strict;
use Carp;
use warnings;
use Cwd;

my $usage = "Incorrect velveth parameter string: Needs to be of the form\n{[-file_format][-read_type] filename}\n";
$usage .= "Where:\n\tFile format options:
        -fasta
        -fastq
        -fasta.gz
        -fastq.gz
        -eland
        -gerald

Read type options:
        -short
        -shortPaired
        -short2
        -shortPaired2
        -long
        -longPaired\n\nThere can be more than one filename specified as long as its a different type.\nStopping run\n";

#constructor
sub new {
    my $class = shift;
    my $self = {@_};
    if($self->{pstringh}){
        my $vhcheck = &_checkVHString($self->{pstringh});
        if(!$vhcheck){
            croak $usage;
        }
    }
    bless ($self, $class);
    return $self;
}

#accessor methods
sub assmscore{ $_[0]->{assmscore}=$_[1] if defined $_[1]; $_[0]->{assmscore}}
sub timestamph{ $_[0]->{timestamph}=$_[1] if defined $_[1]; $_[0]->{timestamph}}
sub timestampg{ $_[0]->{timestampg}=$_[1] if defined $_[1]; $_[0]->{timestampg}}
sub ass_id{ $_[0]->{ass_id}=$_[1] if defined $_[1]; $_[0]->{ass_id}}
sub versionh{ $_[0]->{versionh}=$_[1] if defined $_[1]; $_[0]->{versionh}}
sub versiong{ $_[0]->{versiong}=$_[1] if defined $_[1]; $_[0]->{versiong}}
sub readfilename{ $_[0]->{readfilename}=$_[1] if defined $_[1]; $_[0]->{readfilename}}
sub pstringh{ $_[0]->{pstringh}=$_[1] if defined $_[1]; $_[0]->{pstringh}}
sub pstringg{ $_[0]->{pstringg}=$_[1] if defined $_[1]; $_[0]->{pstringg}}
sub ass_dir{ $_[0]->{ass_dir}=$_[1] if defined $_[1]; $_[0]->{ass_dir}}
sub hashval{ $_[0]->{hashval}=$_[1] if defined $_[1]; $_[0]->{hashval}}
sub rmapfs{ $_[0]->{rmapfs}=$_[1] if defined $_[1]; $_[0]->{rmapfs}}
sub nconts{ $_[0]->{nconts}=$_[1] if defined $_[1]; $_[0]->{nconts}}
sub n50{ $_[0]->{n50}=$_[1] if defined $_[1]; $_[0]->{n50}}
sub maxlength{ $_[0]->{maxlength}=$_[1] if defined $_[1]; $_[0]->{maxlength}}
sub nconts1k{ $_[0]->{nconts1k}=$_[1] if defined $_[1]; $_[0]->{nconts1k}}
sub totalbp{ $_[0]->{totalbp}=$_[1] if defined $_[1]; $_[0]->{totalbp}}
sub totalbp1k{ $_[0]->{totalbp1k}=$_[1] if defined $_[1]; $_[0]->{totalbp1k}}
sub velvethout{ $_[0]->{velvethout}=$_[1] if defined $_[1]; $_[0]->{velvethout}}
sub velvetgout{ $_[0]->{velvetgout}=$_[1] if defined $_[1]; $_[0]->{velvetgout}}
sub sequences{ $_[0]->{sequences}=$_[1] if defined $_[1]; $_[0]->{sequences}}

#velveth parameter string checking function.
sub _checkVHString {
    my %fileform = ();
    my %readform = ();

    my @Fileformats = qw(-fasta -fastq -fasta.gz -fastq.gz -eland -gerald);
    my @Readtypes = qw(-short -shortPaired -short2 -shortPaired2 -long -longPaired);

    foreach(@Fileformats){ $fileform{$_} = 1;}
    foreach(@Readtypes){ $readform{$_} = 1;}

    my $line = shift;
    my @l = split /\s+/, $line;

    #first check for a directory name as the first parameter...
    my $dir = shift @l;
    if(!($dir =~ /\w+/) || ($dir =~ /^\-/)){
        carp "**** $line\n\tNo directory name specified as first parameter in velveth string. Internal error!\n";
        return 0;
    }
    #print "VH Check passed directory..\n";
    my $hash = shift @l;
    unless($hash =~ /^\d+$/){
        carp "**** $line\n\tHash value in velveth string not a number. Internal error!\n";
        return 0;
    }

    #print "VH check passed hash value..\n";

    my $i = 0;
    my $ok = 1;
    foreach(@l){
        #print $_ . "\n";
        if(/^-/){
            #print "Got a dash..\n";
            #s/-//;
            if(!$fileform{$_} && !$readform{$_}){
                carp "**** $line\n\tIncorrect fileformat or readformat specified.\n\t$_ is an invalid velveth switch.\n";
                return 0;
            }
            elsif($fileform{$_}){
                if(($i + 1) > $#l){
                    carp "$line\n\tNo filename supplied after file format type $l[$i].\n";
                    return 0;
                }
                if($readform{$l[$i+1]}){
                    if(($i+2) > $#l){
                        carp "$line\n\tNo filename supplied after read format type $l[$i+1].\n";
                        return 0;
                    }
                    if(-e $l[$i+2]){
                        $ok = 1;
                    }
                    else{
                        carp "**** 320 $line\n\tVelveth filename " . $l[$i+2] . " doesn't exist.\n";
                        return 0;
                    }
                }
                elsif (-e $l[$i+1]){
                    $ok = 1;
                }
                else {
                   carp "**** 328 $line\n\tVelveth filename " . $l[$i+1] . " doesn't exist.\n";
                    return 0;
                }
            }
            elsif($readform{$_}){
                if(($i + 1) > $#l){
                    carp "$line\n\tNo filename supplied after read format type $l[$i].\n";
                    return 0;
                }
                #print "i + 1: " . $l[$i+1] . "\n";
                if($fileform{$l[$i+1]}){
                    if(($i+2) > $#l){
                        carp "$line\n\tNo filename supplied after file format type $l[$i+1].\n";
                        return 0;
                    }
                    if(-e $l[$i+2]){
                        $ok = 1;
                    }
                    else{
                        carp "**** 346 $line\n\tVelveth filename " . $l[$i+2] . " doesn't exist.\n";
                        return 0;
                    }
                }
                elsif (-e $l[$i+1]){
                    $ok = 1;
                }
                else {
                    carp "**** 354 $line\n\tVelveth filename " . $l[$i+1] ." doesn't exist.\n";
                    return 0;
                }
            }
        }
        elsif(!-e $_){

            carp "**** 361 $line\n\tVelveth filename $_ doesn't exist.\n";
            return 0;
        }
        $i ++;
    }
    if($ok){
        return 1;
    }
}
#assemblyScoreCalculator
sub calcAssemblyScore {
    my $self = shift;
    unless(!$self->nconts1k || !$self->n50){
        #$self->{assmscore} = $self->totalbp1k / $self->nconts1k * $self->n50 / 10000;
        $self->{assmscore} = $self->totalbp1k;
        return 1;
    }
    $self->{assmscore} = -1;
    return 0;
}

#getHashingDetails
sub getHashingDetails {
    my $self = shift;
    unless(!$self->timestamph || !$self->pstringh){
        my $programPath = cwd;
        $self->pstringh =~ /^(\w+)\s+(\d+)\s+(.*)$/;
        $self->{ass_dir} = $programPath . "/" . $1;
        $self->{rmapfs} = -s $self->ass_dir . "/Roadmaps";
        $self->{hashval} = $2;
        $self->{readfilename} = $3;
        my @t = split /\n/, $self->velvethout;
        foreach(@t){
            if(/^(\d+).*total\.$/){
                $self->{sequences} = $1;
                last;
            }
        }
        return 1;
    }
    return 0;
}

#getAssemblyDetails
sub getAssemblyDetails {
    my $self = shift;
    unless(!$self->timestampg || !$self->pstringg || !$self->velvetgout){
        my @t = split /\n/, $self->velvetgout;
        foreach(@t){
            if(/^Final graph/){
                /(\d+) nodes/;
                $self->{nconts} = $1;
                /n50 of (\d+), max (\d+)/;
                $self->{n50} = $1;
                $self->{maxlength} = $2;
                last;
            }
        }
        unless(!(-e $self->ass_dir . "/stats.txt")){
            open IN, $self->ass_dir . "/stats.txt";
            my $b = 0;
            my $l = 0;
            my $tb = 0;
            while(<IN>){
                if(/^\d+\s+(\d+)/){
                    my $c = $1;
                    if($c >= 1000){
                        $l ++;
                        $b += $c;
                    }
                    $tb += $c;
                }
            }
            $self->{nconts1k} = $l;
            $self->{totalbp} = $tb;
            $self->{totalbp1k} = $b;
            close IN;
        }
        $self->calcAssemblyScore();

        return 1;
    }
    return 0;
}

#toString method
sub toString {
    my $self = shift;
    my $tmp = $self->toStringNoV();
    if(defined $self->velvethout){
        $tmp .= "Velveth Output:\n" . $self->velvethout() . "\n";
    }
    if(defined $self->velvetgout){
        $tmp .= "Velvetg Output:\n" . $self->velvetgout() . "\n";
    }
    $tmp .= "**********************************************************\n";
    return $tmp;
}


#toStringNoV method
sub toStringNoV {
    my $self = shift;
    my $tmp = "********************************************************\n";
    if($self->ass_id()){
        $tmp .= "Assembly id: " . $self->ass_id(). "\n";
    }
    if($self->assmscore()){
        $tmp .= "Assembly score: " .$self->assmscore(). "\n";
    }
    if($self->timestamph()){
        $tmp .= "Velveth timestamp: " . $self->timestamph(). "\n";
    }
    if($self->timestampg()){
        $tmp .= "Velvetg timestamp: " . $self->timestampg(). "\n";
    }
    if(defined $self->versionh){
        $tmp .= "Velveth version: " . $self->versionh(). "\n";
    }
    if(defined $self->versiong){
        $tmp .= "Velvetg version: " . $self->versiong(). "\n";
    }
    if(defined $self->readfilename){
        $tmp .= "Readfile(s): " . $self->readfilename(). "\n";
    }
    if(defined $self->pstringh){
        $tmp .= "Velveth parameter string: " . $self->pstringh(). "\n";
    }
    if(defined $self->pstringg){
        $tmp .= "Velvetg parameter string: " . $self->pstringg(). "\n";
    }
    if(defined $self->ass_dir){
        $tmp .= "Assembly directory: " . $self->ass_dir(). "\n";
    }
    if(defined $self->hashval){
        $tmp .= "Velvet hash value: " . $self->hashval(). "\n";
    }
    if(defined $self->rmapfs){
        $tmp .= "Roadmap file size: " . $self->rmapfs(). "\n";
    }
    if(defined $self->sequences){
        $tmp .= "Total number of sequences: " . $self->sequences(). "\n";
    }
    if(defined $self->nconts){
        $tmp .= "Total number of contigs: " . $self->nconts(). "\n";
    }
    if(defined $self->n50){
        $tmp .= "n50: " . $self->n50(). "\n";
    }
    if(defined $self->maxlength){
        $tmp .= "length of longest contig: " . $self->maxlength(). "\n";
    }
    if(defined $self->totalbp){
        $tmp .= "Total bases in contigs: " . $self->totalbp(). "\n";
    }
    if(defined $self->nconts1k){
        $tmp .= "Number of contigs > 1k: " . $self->nconts1k(). "\n";
    }
    if(defined $self->totalbp1k){
        $tmp .= "Total bases in contigs > 1k: " . $self->totalbp1k(). "\n";
    }
    $tmp .= "**********************************************************\n";
    return $tmp;
}

1;
