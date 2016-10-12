#!/bin/perl

### Blast parser + building pseudo-scaffolds of mtDNA (used by RunDCMegaBlast.pl)
### Deals with one specimen at a time
# Usage:
# perl ParseBlast_append.pl blastfile contigfile outfile minpct minlen append mincov bckbonelen
#
# blastfile = blast hits (-outfmt 6, as obtained with Nucleotide-Nucleotide BLAST 2.2.29+)
# contigfile = contigs obtained via assembly (fasta)
# outfile = outfile (fasta)
# minpct = minimum % similarity between contig and reference
# minlen = minimum length (bp) of region matching in contig on reference
# append = 0|1 Concatenate (append = 1) or store hits in distinct sequences (append = 0)
# mincov = 5, number of reads supporting the considered contig (SOAP outputs only, disabled otherwise)
# bckbonelen = length of template sequence (N's) that will be filled with actual data. Set to something longer than a mitonchondria (e.g. 20000)
#
# Nils Arrigo, Unil 2015
####
use List::Util qw(min max);
use File::Basename;

### Get script args
my $blastfile = $ARGV[0];
my $contigfile = $ARGV[1];
my $outfile = $ARGV[2];
my $minpct = $ARGV[3];
my $minlen = $ARGV[4];
my $append = $ARGV[5];
my $mincov = $ARGV[6];
my $bckbonelen = $ARGV[7]; #length of mitochondrial backbone
chomp($bckbonelen);


### hard coded params
my $format = "spades"; #or SOAP



my $accession = basename($contigfile);

### Load blast tabular
open(BLAST, "$blastfile");
open(OUT, ">>$outfile");
open(COVS, ">>$outfile.covs");

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
  my $realstart = min($start, $stop);
  my $idx = min($start, $stop);
  my $strand = $fields[12];

  # check whether we need to rev-comp  
  my $modif = "ok";
  if($start > $stop){
    $modif = "rev";
    }

  if($strand eq "minus"){
    $modif = "$modif.comp";
    }


  # Identify interesting hits
  if($pct > $minpct and $len > $minlen){
    # check if we have already seen  that guy
    if($previouslen = $track{$ctg}){ #yes, need to check length
      if($previouslen < $len){ #if longer hit, keep this one
    	$goodhits{$idx}{"name"} = $ctg;
	$goodhits{$idx}{"len"} = $len;
	$goodhits{$idx}{"pct"} = $pct;
	$goodhits{$idx}{"revcomp"} = $modif;
	$goodhits{$idx}{"qstart"} = $qstart;
	$goodhits{$idx}{"bstart"} = $realstart;
	$track{$ctg} = $len; # keep track of visited contigs
	}
      } else { # no, go ahead
      $goodhits{$idx}{"name"} = $ctg;
      $goodhits{$idx}{"len"} = $len;
      $goodhits{$idx}{"pct"} = $pct;
      $goodhits{$idx}{"revcomp"} = $modif;
      $goodhits{$idx}{"qstart"} = $qstart;
      $goodhits{$idx}{"bstart"} = $realstart;
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
my %covs;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  
  # take of accession (and coverages)
  if($format =~ /SOAP/){
    @acc = split(/\s/, $tmp[0], 2);
    $acc = trim($acc[0]);
    $cov = trim($acc[1]);
    }

  if($format =~ /spades/){
    $acc = $tmp[0];
    $acc =~ /.*\_cov\_(.*)\_ID.*/;
    $cov = $1;
    print "coverage = $cov\n";
    }

  
  # get sequence
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = trim($_);
  $fasta{$acc} = $seq;	
  $covs{$acc} = $cov;
  }

# print "Loading $contigfile: done\n";

### filter blast hits, in a way to remove contigs with too low coverage
foreach my $idx (keys(%goodhits)){
  $acc = $goodhits{$idx}{"name"};
  $cov= $covs{$acc};
  if($cov < $mincov){
    delete $goodhits{$idx};
    print "skip contig $acc, coverage too low: $cov\n";
    }
  }


### Produce output
my @pos = keys(%goodhits);
my @sortedpos = sort { $a <=> $b } @pos;
my $cnt = 0;

my $backbone = "N" x $bckbonelen;
my @backbone = split("", $backbone);

foreach my $idx (@sortedpos){
  my $acc = $goodhits{$idx}{"name"};
  my $revcomp = $goodhits{$idx}{"revcomp"};
  my $bstart = $goodhits{$idx}{"bstart"};
  my $qstart = $goodhits{$idx}{"qstart"};
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
      if($revcomp =~ /comp/){
	$finalseq = comp($finalseq);
	}
      if($revcomp =~ /rev/){
	$finalseq = reverse($finalseq);
	}
      $goodhits{"finalchunck"}{"seq"} = "$finalseq";
      $goodhits{"finalchunck"}{"bstart"} = $bckbonelen - length($finalseq);
      }
    $cnt++;
    } else {
    $finalseq = $seq;
    }

  if($revcomp =~ /comp/){
    $finalseq = comp($finalseq);
    }
  if($revcomp =~ /rev/){
    $finalseq = reverse($finalseq);
    }
  
  my @insert = split("", $finalseq);
  my $offset = $bstart - $qstart - 1;
  if($offset < 0){
    $offset = 1;
    }
  my $end = $offset + $#insert;
  @backbone[$offset..$end] = @insert;

  $cov = $covs{$acc};
  print COVS "$cov\t";
  }
undef(%fasta);


## finish with final chunck
if(my $finalseq = $goodhits{"finalchunck"}{"seq"}){
  my $bstart = $goodhits{"finalchunck"}{"bstart"};
  my @insert = split("", $finalseq);
  my $offset = $bstart;
  if($offset < 0){
    $offset = 1;
    }
  my $end = $offset + $#insert;
  @backbone[$offset..$end] = @insert;

  $cov = $covs{$acc};
  print COVS "$cov\t";
  }

my $finalseq = join("", @backbone);
print OUT ">$accession\n$finalseq\n";
close(OUT);
print COVS "$cov\n";
close(COVS);

# print "Saving to $outfile: done\n";




##### Routines
## remove whitespaces
sub trim($){
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
  }

sub comp {
  my $dna = shift;
  $dna =~ tr/ACGTacgt/TGCAtgca/;
  return $dna;
  }

