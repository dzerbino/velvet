#
#       VelvetOpt::Utils.pm
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

package VelvetOpt::Utils;

use strict;
use warnings;
use POSIX qw(ceil floor);
use Carp;

#no longer used!
#my $upperLim = 500; #the upper limit of coverage to consider for an exp_cov determination.
#my $lowerLim = 5; #the lower limit of coverage to consider for an exp_cov determination.
#gets rid of values that are too high in the instance of the 16s etc and too low in the
#instance of long reads with no short read coverage to avoid the whole exp_cov = 0 thing...


my $min_len = 100; #the minimum length of the contigs to consider for the expected coverage
#estimation.

#   histogramMaxBin:
#   based on a version of histogram.pl posted 2007 by Christian at
#   http://snippets.aktagon.com/snippets/62-How-to-generate-a-histogram-with-Perl

#   it returns the floor of the bin with the most hits in the histogram.
#   it requires two args.  first the bin width and second an array of the numbers to histogram.

sub histogramMaxBin {

    my ($bin_width, $list) = @_;

    # This calculates the frequencies for all available bins in the data set
    my %histogram;
    $histogram{ceil(($_ + 1)/$bin_width) -1}++ for @$list;

    my $max;
    my $min;

    # Calculate min and max
    while (my($key, $value) = each(%histogram)){
        $max = $key if !defined($min) || $key > $max;
        $min = $key if !defined($min) || $key < $min;
    }

    my $max_count = 0;
    my $max_count_bin = 0;

    for (my $i = $min; $i <= $max; $i++){
        my $frequency = $histogram{$i} || 0;
        if($frequency > $max_count){
            $max_count = $frequency;
            $max_count_bin = $i;
        }

    }
    return $max_count_bin*$bin_width;
}

#   getExpCov
#   it returns the expected coverage of short reads from an assembly by
#   performing a histogram on the short1_cov column in the stats.txt from
#   the assembly directory.

sub getExpCov {
    my $file = shift;
    my @nums;
    unless(open IN, $file){ croak "Unable to open $file for exp_cov determination.\n";}
    while(<IN>){
        my @tmp = split /\s+/, $_;
        if($_ =~ /^\d+/ && $tmp[1] > $min_len){
            push @nums, $tmp[5];
        }
    }
    return &histogramMaxBin(1, \@nums);
}
