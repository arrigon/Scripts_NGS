#!/bin/perl

#####################################
#### Specific to DMPS pipeline
#### Collapses sequences from a same specimen, but sotred in distinct chuncks into a single consensus
#### e.g. 
####	>seqA.chk1
####	ATTCGC---------
####	>seqA.chk2
####	----------AATTCG
####	becomes
####
####	>seqA
####	ATTCGC----AATTCG
####################################
#### usage perl ALN_ConsenseChuncks.pl infile outfile
####
#### Assumes that sequence name contains chunck identifiers (sequences names are ending with .chk1, )
####
#### Nils Arrigo, UNIL 2014
#####################################

#load packages
use File::Basename;
use POSIX;
use List::MoreUtils qw(uniq);


# Parameters
my $infile = $ARGV[0];
my $outfile = $ARGV[1];



### Open fasta and load it into a hash
open(FILE, "$infile");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

my %fasta;
my @inds; #keep track of accession input order
my $seqlen = 0;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  my $seq = $_;
  $seq = uc($seq);
  $seq =~ s/-/N/g;
  
  my $acc = $tmp[0];
  my $ind = $acc;
  $ind =~ s/\.chk\d+$//;
  
  $fasta{$ind}{$acc} = $seq; 
  push(@inds, $ind);
  $seqlen = length($seq);
  }

@inds = uniq(@inds);


  
### Merge sequences form same individual
my %fastaclean;
foreach $ind (keys(%fasta)){
  my %tmp = %{$fasta{$ind}};
  my @consensus = ("N") x $seqlen;

  foreach $chk (keys(%tmp)){
    my $seq = $fasta{$ind}{$chk};
    my @seq = split("", $seq);
    
    foreach $pos (0..$#seq){
      my $bp = $seq[$pos];
      if($bp !~ /[N|n|-]/){
	$consensus[$pos] = $bp;
	}
      }
      
    my $seqclean = join("", @consensus);
    $fastaclean{$ind} = $seqclean;
    }
  }



### save to outfile
open(OUT, ">$outfile");
foreach $ind (@inds){
  my $seq = $fastaclean{$ind};
  print OUT ">$ind\n$seq\n";
  }
