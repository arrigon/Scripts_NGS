#!/bin/perl
#####################################
#### Complement / Reverse / reverseComplement a sequence
####
####  perl revcon fasta.in fasta.out operation
####  - operation = rev, comp, revcomp
#####################################

use File::Basename;

my $infile = $ARGV[0];
my $operation = $ARGV[1];

chomp($infile);
chomp($operation);

$bsn = $infile;
$bsn =~ s/\..*$//;


## import fasta file
open(FILE, "$infile");
open(OUT, ">$bsn.$operation.txt");

while(<FILE>){
  chomp();
  my $line = $_;
  my @fields = split(/\t/, $line, 2);
  my $specimen = $fields[0];
  my $seq = $fields[1];

  $_ = $seq;
  s/\r|\n//g;
  $seq = uc $_; #turn everything in uppercase

  if($operation =~ /comp/){
    $seq =~ tr/ATCG/TAGC/;
    }

  if($operation =~ /rev/){
    $seq = reverse($seq);
    }
  
  print(OUT "$specimen\t$seq\n");
  }
close(FILE);
close(OUT);
