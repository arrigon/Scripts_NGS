#!/bin/perl
#####################################
#### Orients properly fasta files, using ublast outputs
#####################################
use File::Basename;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $ucfile = $ARGV[2];

chomp($infile);
chomp($outfile);
chomp($ucfile);

$bsn = basename($infile);
$bsn =~ s/\..*$//;


#### load ublast uc output
open IN, $ucfile;
my %listID;
while(<IN>){
  chomp();
  my @fields = split("\t", $_);
  my $marker = $fields[8];
  my $rc = $fields[4];
  $listID{$marker} = $rc;
  }


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
  $fasta{$tmp[0]} = $seq; 
  push(@accs, $tmp[0])
  }


## print sequences of interest into output
open(OUT, ">$outfile");
my $cnt = 0;
foreach $head (@accs){
  chomp($head);
  my $orient = $listID{$head};  
  my $seq = $fasta{$head};
  if($orient =~ /\-/){
    $seq =~ tr/ATCG/TAGC/;
    $seq = reverse($seq);
    }

  if($orient =~ /\-|\+/){
    print(OUT ">$head\n$seq\n");
    $cnt++;
    }
  }
close(OUT);
