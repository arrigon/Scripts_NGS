#!/bin/perl
#####################################
#### Filters out fasta file, based on coverage stats (Uclust output format)
#####################################

## Scripts params
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $mincov = $ARGV[2];
my $minfrac = $ARGV[3];
chomp($minfrac);


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
my @covs;
foreach my $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = uc $_; #turn everything in uppercase
  $fasta{$tmp[0]} = $seq;
  push(@accs, $tmp[0]);

  # get coverages
  my $head = $tmp[0];
  $head =~ /.*\_cov=(.*)/;
  my $cov = $1;
  push(@covs, $cov);
  }
my $best = $covs[0];
my $lwr = $best * $minfrac;


## print sequences of interest into output
open(OUT, ">$outfile");
my $cnt = 0;
foreach $head (@accs){
  chomp($head);
  $head =~ /.*\_cov=(.*)/;
  my $cov = $1;
  my $seq = $fasta{$head};
  
  if($cov >= $lwr && $cov >= $mincov){
    print(OUT ">$head\n$seq\n");
    }
  }
close(OUT);
