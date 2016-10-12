#!/bin/perl

use File::Basename;

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $outinfo = $ARGV[2];


## Prepare IO
chomp($outinfo);
my $bsn = basename($infile);
$bsn =~ s/\>/\./g;


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
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $acc = trim($tmp[0]);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = trim($_);
  $fasta{$acc} = $seq; 
  push(@accs, $acc);
  }


## print sequences of interest into output
my $seqcnt = 0;


open(OUT, ">$outfile");
open(INF, ">$outinfo");
foreach $head (@accs){
  chomp($head);
  my $seq = $fasta{$head};
  if(defined($seq)){
    print(OUT ">$bsn.$seqcnt\n$seq\n");
    print(INF "$bsn.$seqcnt\t$head\n");
    $seqcnt++
    }
  }
close(OUT);
close(INF);


## remove whitespaces
sub trim($){
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

