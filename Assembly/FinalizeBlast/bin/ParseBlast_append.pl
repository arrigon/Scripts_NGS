#!/bin/perl

### Blast parser + building pseudo-scaffolds of mtDNA (used by RunDCMegaBlast.pl)
### Deals with one specimen at a time
# Usage:
# perl ParseBlast_append.pl blastfile contigfile outfile append
#
# blastfile = blast hits (-outfmt 6, as obtained with Nucleotide-Nucleotide BLAST 2.2.29+)
# contigfile = contigs obtained via assembly (fasta)
# outfile = outfile (fasta)
# append = 0|1 Concatenate (append = 1) or store hits in distinct sequences (append = 0)
#
# Nils Arrigo, Unil 2013
####

use List::Util qw(min max);
use File::Basename;


my $blastfile = $ARGV[0];
my $contigfile = $ARGV[1];
my $outfile = $ARGV[2];
my $append = $ARGV[3];
chomp($append);

my $accession = basename($contigfile);


### Load blast tabular
open(BLAST, "$blastfile");
open(OUT, ">>$outfile");

my %goodhits;
my %track;

while(<BLAST>){
  chomp();
  my @fields = split("\t", $_);

  # Get data
  my $ctg = $fields[0];
  my $pct = $fields[2];
  my $len = $fields[3];
  my $start = $fields[8];
  my $stop = $fields[9];
  my $qstart = $fields[6];
  my $idx = min($start, $stop);

  # check whether we need to rev-comp  
  my $needrevcomp = 0;
  if($start > $stop){
    $needrevcomp = 1;
    }

  # Identify interesting hits
  if($pct > 50 and $len > 100){
    # check if we have already seen  that guy
    if($previouslen = $track{$ctg}){ #yes, need to check length
      if($previouslen < $len){ #if longer hit, keep this one
    	$goodhits{$idx}{"name"} = $ctg;
	$goodhits{$idx}{"len"} = $len;
	$goodhits{$idx}{"pct"} = $pct;
	$goodhits{$idx}{"revcomp"} = $needrevcomp;
	$goodhits{$idx}{"qstart"} = $qstart;
	$track{$ctg} = $len; # keep track of visited contigs
	}
      } else { # no, go ahead
      $goodhits{$idx}{"name"} = $ctg;
      $goodhits{$idx}{"len"} = $len;
      $goodhits{$idx}{"pct"} = $pct;
      $goodhits{$idx}{"revcomp"} = $needrevcomp;
      $goodhits{$idx}{"qstart"} = $qstart;
      $track{$ctg} = $len; # keep track of visited contigs
      }
    }
  }
# print "Loading $blastfile: done\n";

close(BLAST);


### Load contigs from fasta file
open(FILE, "$contigfile");
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
  @acc = split(/\s/, $tmp[0], 2);
  $acc = trim($acc[0]);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = trim($_);
  $fasta{$acc} = $seq; 
  }

# print "Loading $contigfile: done\n";


### Produce output
my @pos = keys(%goodhits);
my @sortedpos = sort { $a <=> $b } @pos;
my $cnt = 0;

if($append == 1){
  print OUT ">$accession\n";
  }
foreach my $idx (@sortedpos){
  my $acc = $goodhits{$idx}{"name"};
  my $revcomp = $goodhits{$idx}{"revcomp"};
  my $seq = $fasta{$acc};
  my $finalseq;

  if($cnt == 0){ # we are looking at the first contig, make sure it starts at pos 1 (because mtDNA is circular)
		 # if we start at pos > 1, then we chop that contig into two chuncks.
		 # the first chunck is kept here as starting sequence, the last chunck is stored back into %goodhits
		 # and will be appended to our sequence at the last iteration
    my $qstart = $goodhits{$idx}{"qstart"};
    if($qstart > 10){
      my @tmp = split("", $seq);
      my $maxpos = $#tmp;
      $seq = join("", @tmp[$qstart..$maxpos]);

#       print "$contigfile - Adjust start to: $qstart\n";
      $finalseq = join("", @tmp[0..($qstart - 1)]);
      if($revcomp == 1){
	$finalseq = reverse_complement($finalseq);
	}
      $goodhits{"finalchunck"}{"seq"} = "$finalseq";
      }
    $cnt++;
    }

  if($revcomp == 1){
    $seq = reverse_complement($seq);
    }

  if($append != 1){
    print OUT ">$accession.$acc.$idx\n$seq\n";
    } else {
    print OUT "$seq\n";
    }
  }

undef(%fasta);

# finish with final chunck
if(my $seq = $goodhits{"finalchunck"}{"seq"}){
  if($append != 1){
  print OUT ">$accession.finalchunck\n$seq";
    } else {
    print OUT "$seq";
    }
  }
print OUT "\n";
close(OUT);


# print "Saving to $outfile: done\n";




##### Routines
## remove whitespaces
sub trim($){
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
  }

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);
# 	my $revcomp = $dna;

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

