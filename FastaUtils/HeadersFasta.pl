#!/bin/perl

#####################################
#### Get fasta headers and save them into an output file
####
####  perl HeadersFasta.pl fastafile
#####################################

my $file1 = $ARGV[0];

## import fasta file
open(FILE, "$file1");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

## load it into hash. Note that %fasta contains your sequences. Can be reused elsewhere.
my %fasta;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = $_;
  $fasta{$tmp[0]} = $seq; 
  }

## print keys to output file (txt file)
open(OUT, ">$file1.head");
print(OUT join("\n", keys(%fasta)), "\n");
close(OUT);

