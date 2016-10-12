#!/bin/perl

#load packages
use File::Basename;
use POSIX;


# Parameters
my $infolder = $ARGV[0];
my $outfile = $ARGV[1];

my @files = `ls $infolder/*.fq`;
open(OUT, ">$outfile");

@files = sort {$a cmp $b} @files;

foreach my $infile (@files){
  chomp($infile);

  # Open fasta and load it into a hash
  open(FILE, "$infile");
  local $/ = undef; #slurp, mode
  my $input = <FILE>;
  local $/ = "\n";
  my @fields = split(/\+/, $input);
  close(FILE);

  my $seq = $fields[0];
  my $acc = basename($infile);
  $acc =~ s/\.bam\.cons\.fq//;
  $seq =~ s/\@.*\n/\>$acc\n/;

  print OUT "$seq";
  }

close(OUT);
