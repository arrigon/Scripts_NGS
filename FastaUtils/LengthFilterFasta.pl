#!/bin/perl

# Usage: Keeps only sequences with length > minlength
#
# perl SliceFasta.pl fastafile nchunks
#
###

use POSIX;

my $file1 = $ARGV[0];
my $minlen = $ARGV[1];

# Open fasta
open(FILE, "$file1");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n"; #slurp, mode
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

my %fasta;
$kept = 0;
$removed = 0;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = $_;
  $len = split("", $seq);
  if($len >= $minlen){
    $fasta{$tmp[0]} = $seq;
    $kept++;
    } else {
    $removed++;
    }
  }
my $init = $kept + $removed;
my $pct = floor(100 * $kept / $init);
print "Length filtering $file1; keeping only sequences greater than $minlen bp\nKept $kept / $removed ($pct \%)\n";

open(OUT, ">minlen\_$minlen\_$file1");
foreach $acc (keys(%fasta)){
  $seq = $fasta{$acc};
  print(OUT ">$acc\n$seq\n");
  }  
close(OUT);

