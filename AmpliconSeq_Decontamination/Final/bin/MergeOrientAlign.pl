#!/bin/perl
#####################################
#### Extracts first sequence of a fasta file
#### prefixes sequence header with file name
#####################################
use File::Basename;

my $infolder = $ARGV[0];
my $outfolder = $ARGV[1];

chomp($infolder);
chomp($outfolder);


my @allfiles = `ls -1 $infolder`;

# list existing markers
my %markers;
foreach my $infile (@allfiles){
  chomp($infile);
  my $bsn = basename($infile);

  $bsn =~ /(.*)\_(.*).fas/;  
  my $mark = $1;
  $markers{$mark} = 1;
  }

# append them
foreach $mark (keys(%markers)){

  # merging
  my $command = "cat tmp/postfilter/$mark* > tmp/merged/$mark.merged";
  print("$command\n");
  system($command);

  # orienting
  my $command = "perl bin/CheckOrientFasta.pl tmp/merged/$mark.merged tmp/merged/$mark.orient";
  print("$command\n");
  system($command); 

  # aligning
  my $command = "muscle -in tmp/merged/$mark.orient -out $outfolder/$mark.fas";
  print("$command\n");
  system($command); 
  }

