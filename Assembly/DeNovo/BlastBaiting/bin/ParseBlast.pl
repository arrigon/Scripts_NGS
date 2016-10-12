#!/bin/perl

use List::Util qw(min max);
use File::Basename;


my $blastfile = $ARGV[0];
my $contigfile = $ARGV[1];
my $outfile = $ARGV[2];

open(OUT, ">>$outfile");
print OUT ">$accession\n";

my $accession = basename($contigfile);


### Load blast tabular
open(BLAST, "$blastfile");

my %goodhits;
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
  if($pct > 85 and $len > 90){
    $goodhits{$idx}{"name"} = $ctg;
    $goodhits{$idx}{"len"} = $len;
    $goodhits{$idx}{"pct"} = $pct;
    $goodhits{$idx}{"revcomp"} = $needrevcomp;
    $goodhits{$idx}{"qstart"} = $qstart;
    }
  }
print "Loading $blastfile: done\n";

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

print "Loading $contigfile: done\n";


### Produce output
my @pos = keys(%goodhits);
my @sortedpos = sort { $a <=> $b } @pos;
my $cnt = 0;

foreach my $idx (@sortedpos){
  my $acc = $goodhits{$idx}{"name"};
  my $revcomp = $goodhits{$idx}{"revcomp"};
  my $seq = $fasta{$acc};
  my $finalseq;

  if($cnt == 0){ #we are looking at the first contig, make sure it starts at pos 1
    my $qstart = $goodhits{$idx}{"qstart"};
    my @tmp = split("", $seq);
    my $maxpos = $#tmp;
    $seq = join("", @tmp[$qstart..$maxpos]);

    $finalseq = join("", @tmp[0..($qstart - 1)]);

    if($revcomp == 1){
      $finalseq = reverse_complement($finalseq);
      $goodhits{"finalchunck"}{"seq"} = $finalseq;
      }
    $cnt++;
    }

  if($revcomp == 1){
    $seq = reverse_complement($seq);
    }

  print OUT "$seq";
  }

undef(%fasta);

# finish with final chunck
my $seq = $goodhits{"finalchunck"}{"seq"};
print OUT "$seq\n";
close(OUT);


print "Saving to $outfile: done\n";




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

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}

