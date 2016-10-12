#!/bin/perl

####################################
# Removing specific base pairs from reads
#
### Params
# - infile = input file (takes fastq files; not fast.gz, in 4 lines format)
# - outfile = output file
# - pos = 1,2,3,4 positions to be removed, 
#	  if deleting several of them, use comma separator. Positions are in bp, not in intervals between bp.
#
### Usage: perl DropBpFastq.pl infile outfile 69,70,71
#
# Nils Arrigo. Uni Lausanne 2012
############################################################
use File::Basename;

### Get script arguments
my $infile = $ARGV[0];
my $bsnfile = basename($infile);

my $outfile = $ARGV[1];

my $pos = $ARGV[2];
chomp($pos);
die "Dropping bases in $infile: Stopped - please define positions to remove\n" if $pos == undef;
my @pos = split(",", $pos);
my @pos = map { $_ - 1 } @pos;


### load fastq file  into %store
open IN1, $infile;
my %store0;
my $cnt = 0;
my $nline = 0;
my $sequencer;
my @accsorder;


while(<IN1>){
  chomp();
  my $line = $_;

  if($nline == 0){ #Use the very first line of fastq to pick sequencer name
    $sequencer = $line;
    $sequencer =~ /(@[\w|\-,]+:\d+:\w+):.*/; # RegExp to catch sequencer name, here @VADER:36:C0MLBACXX @VADER:36:C0MLBACXX:8:2314:17186:5629 2:N:0:
    $sequencer = $1;
    }

  if($line =~ /^$sequencer/){ #we enter into a new accession
    $acc = $line; # get accession number
    $acc =~ /(^@.*)\s+.*/; #RegExp to clean accession name
    $acc = $1;
    $cnt = 0;
    push(@accsorder, $acc);
    }

  if($cnt < 4){ # store accession as long as there are < 4 lines (hence skip any weird 5th lines, if any)
    if($cnt == 1 | $cnt == 3){
      my @char = split("", $line);
      my $ref = \@char;
      $store0{$acc}->{$cnt} = $ref;
      } else {
      $store0{$acc}->{$cnt} = $line;
      }
    $cnt++;  
    }
  $nline++;
  }
close(IN1);


### Write fastq to output, but remove non-desired positions
open OUT, ">$outfile";

foreach my $acc (@accsorder){
  # retrieve seq and qual
  my $accline = $store0{$acc}{0};
  my @seq = @{$store0{$acc}{1}};
  my @qual = @{$store0{$acc}{3}};

  # cut out specific positions
  @seq[@pos] = ("") x ($#pos + 1);
  @qual[@pos] = ("") x ($#pos + 1);

  # convert back to string
  my $seq = join("", @seq);
  my $qual = join("", @qual);

  # print to output
  print OUT "$accline\n$seq\n+\n$qual\n";
  }
close(OUT);
