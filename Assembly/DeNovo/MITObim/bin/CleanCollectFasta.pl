#!/bin/perl


use File::Basename;

my $infolder = $ARGV[0];
my $outfile = $ARGV[1];


my @files = `ls $infolder/*.fasta`;
open(OUT, ">$outfile");

@files = sort {$a cmp $b} @files;

foreach my $infile (@files){
  chomp($infile);

  my $bsn = basename($infile);
  $bsn =~ s/\.fasta$//;

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
  foreach $input (@fields){
    my @tmp = split(/\n/, $input, 2);
    $_ = $tmp[1];
    s/\r|\n//g;
    $seq = uc $_; #turn everything in uppercase
    $fasta{$tmp[0]} = $seq; 
    }

  ## print sequences of interest into output
  foreach $head (keys(%fasta)){
    chomp($head);
    my $seq = $fasta{$head};
    $seq =~ s/[x|X]/N/g;
    print(OUT ">$bsn\n$seq\n");
    }
  }

close(OUT);
