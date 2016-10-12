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
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);


## load it into hash. Note that %fasta contains your sequences. Can be reused elsewhere.
my %fasta;
my @accs;
foreach my $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = uc $_; #turn everything in uppercase

  if($operation =~ /comp/){
    $seq =~ tr/ATCG/TAGC/;
    }

  if($operation =~ /rev/){
    $seq = reverse($seq);
    }
  $fasta{$tmp[0]} = $seq; 
  push(@accs, $tmp[0])
  }


## print sequences of interest into output
open(OUT, ">$bsn.$operation.fas");
foreach $head (@accs){
  chomp($head);
  my $seq = $fasta{$head};
  print(OUT ">$head\n$seq\n");
  }
close(OUT);
