#!/bin/perl

#####################################
#### Extracts randomly N sequences from Fasta file
####
####  perl SliceRandFasta.pl fastafile nseq
#####################################

use POSIX;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

my $file1 = $ARGV[0];
my $nseq = $ARGV[1];

# Open fasta
open(FILE, "$file1");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

my %fasta;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = $_;
  $fasta{$tmp[0]} = $seq; 
  }

# Save it into $splits chuncks
my @accs = keys(%fasta);
@accs = shuffle(@accs);
my @selection = @accs[0..$nseq-1];

print "Randomly extracting $nseq sequences from $file1\n";

open(OUT, ">rand$nseq\_$file1");
foreach $acc (@selection){
  $seq = $fasta{$acc};
  print(OUT ">$acc\n$seq\n");
  }  
close(OUT);
