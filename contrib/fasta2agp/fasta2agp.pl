#!/usr/bin/perl

### david.studholme@tsl.ac.uk

### Generates FastA and AGP from Velvet 'contigs.fa' file

### Use entirely at you own risk!!

use strict;
use warnings ;
use Bio::SeqIO ;

my $sequence_file = shift or die "Usage: $0 <sequence file>\n" ;


### Output file for contigs in Fasta format
my $fasta_outfile = "$sequence_file.contigs.fna";
open (FILE, ">$fasta_outfile") and
   warn "Will write contigs to file '$fasta_outfile'\n" or
   die "Failed to write to file '$fasta_outfile'\n";


my $inseq = Bio::SeqIO->new('-file' => "<$sequence_file",
               '-format' => 'fasta' ) ;

while (my $seq_obj = $inseq->next_seq ) {

   my $supercontig_id = $seq_obj->id ;
   my $supercontig_seq = $seq_obj->seq ;
   my $supercontig_desc = $seq_obj->description ;
   my $supercontig_length = length($supercontig_seq);


   my $i = 0;
   my $start_pos = 1;
     my %contig_sequences;
   foreach my $contig_sequence ( split /(N+)/i, $supercontig_seq ) {
   $i++;
     ### Define the AGP column contents
   my $object1 = $supercontig_id;
   my $object_beg2 = $start_pos;
   my $object_end3 = $start_pos + length($contig_sequence) - 1;
   my $part_number4 = $i;
   my $component_type5;
   my $component_id6a;
   my $gap_length6b;
   my $component_beg7a;
   my $gap_type7b;
   my $component_end8a;
   my $linkage8b;
   my $orientation9a;
   my $filler9b;
     if (  $contig_sequence =~ m/^N+$/ ) {
       ### This is poly-N padding
       $component_type5 = 'N';
       $gap_length6b = length($contig_sequence);
       $gap_type7b = 'fragment';
       $linkage8b = 'no';
       $filler9b = '';
         } elsif ( $contig_sequence =~ m/^[ACGT]+$/ ) {
       ### This is a contig
       $component_type5 = 'W';
       $component_id6a = "contig_$supercontig_id.$i";
       $component_beg7a = 1;
       $component_end8a = length($contig_sequence);
       $orientation9a = '+';
             ### Print FastA formatted contig
       print FILE ">$component_id6a\n$contig_sequence\n";


   } else {
       die "Illegal characters in sequence\n$contig_sequence\n";
   }
     $start_pos += length ($contig_sequence);
     if ($component_type5 eq 'N') {
       ### print AGP line for gap
       print "$object1\t$object_beg2\t$object_end3\t$part_number4\t$component_type5\t$gap_length6b\t$gap_type7b\t$linkage8b\t$filler9b\n";
   } else {
       ### print AGP line for contig
       print "$object1\t$object_beg2\t$object_end3\t$part_number4\t$component_type5\t$component_id6a\t$component_beg7a\t$component_end8a\t$orientation9a\n";
         }

   }
    }

