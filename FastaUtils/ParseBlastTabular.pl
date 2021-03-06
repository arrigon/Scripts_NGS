#!usr/bin/perl

#####################################
#### Parse Blast report and keep only hits with XXX length and YYY similarity percentage
####  Works from blast reports produced using the "-outfmt 6" option
####
####  perl ParseBlastTabular blastreport outfile length percent
#####################################
# Erik R Hanschen, modified Nils Arrigo
# Written: Jan 21, 2012

use warnings;
use strict;

# This script parses a Blast 2.2 output format #6 table to only include non-self 'good' hits.
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $length = $ARGV[2];
my $percent = $ARGV[3];

if (!defined $infile || !defined $outfile || !defined $length || !defined $percent) {
  die "Usage: perl $0 input_file output_file match_length percent_match\n";
  }

open (INFILE, $infile) or die "Couldn't open file $infile: $!\n";
open (OUTFILE, ">$outfile") or die "Couldn't open file $outfile: $!\n";

print "Working...\n";

while (my $line = <INFILE>) {
  chomp $line;
  my @array = split (/\t/, $line);
  if (($array[0] ne $array[1]) and ($array[2] gt $percent) and ($array[3] gt $length)) {
    my $orient = $array[9] - $array[8];
      if($orient >= 0){
      $orient = 1;
      } else {
      $orient = -1;
      }
    print OUTFILE "$array[0]\t$array[1]\t$orient\n";
    }
  }

print "Done!\n";
close INFILE;
close OUTFILE;
