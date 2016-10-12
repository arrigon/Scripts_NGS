#!/bin/perl

#####################################
#### Takes a fasta file and its corresponding blast output as input (typical tabular format, outfmt = 8)
#### 1. check blast hits to bin fasta sequences in non-overlapping regions (i.e. chuncks)
#### 2. produces one fasta file per region
#### 
####################################
#### usage perl Blast_DefineBins.pl fasta ref blastfile outfolder
####
#### Nils Arrigo, UNIL 2013
#####################################

#load packages
use File::Basename;
use POSIX;


# Parameters
my $infile = $ARGV[0];
my $reffile = $ARGV[1];
my $blastfile = $ARGV[2];
my $outfolder = $ARGV[3];



# Prepare workfolder
system("mkdir $outfolder");



# Open fasta and load it into a hash
open(FILE, "$infile");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

my %fasta;
my @accs; #keep track of accession input order
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = $_;
  $seq = uc($seq);
  $seq =~ s/-/N/g;
  $fasta{$tmp[0]} = $seq; 
  push(@accs, $tmp[0]);
  }


# Open fasta and load it into a hash
my $refname;
open(FILE, "$reffile");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

my %ref;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = $_;
  $seq = uc($seq);
  $seq =~ s/-/N/g;
  $ref{$tmp[0]} = $seq; 
  $refname = $tmp[0];
  }


  
# Open blast file and load it into a hash
open(FILE, $blastfile);

my %coords;
my %track;
my $currentacc = "";
my $cnt = 0;

while(<FILE>){
  chomp();
  my $line = $_;
  my @fields = split(/\t/, $line);
  my $acc = $fields[0];
  my $len = $fields[3];
  my $startraw = $fields[8];
  my $stopraw = $fields[9];
  
  my $start;
  my $stop;
      
  if($startraw >= $stopraw){
    $start = $stopraw;
    $stop = $startraw;
    } else {
    $start = $startraw;
    $stop = $stopraw;      
    }
      
  $coords{$cnt}{"acc"} = $acc;
  $coords{$cnt}{"len"} = $len;
  $coords{$cnt}{"start"} = $start;
  $coords{$cnt}{"stop"} = $stop;
  
  $track{$acc}{$cnt} = $len;
  $cnt++;
  }
close(FILE);
  
print "Loading files... done\n";



# parse %track and retrieve longest blast hit of each accession
my @tokeep;
my $tmp;
my $bestlen;
my $bestacc;
  
foreach $acc (@accs){
  $bestlen = 0;
  $tmp = -5;
  my @cnts = keys(%{$track{$acc}});
  foreach $cnt (@cnts){
    my $len = $track{$acc}{$cnt};
    if($len > $bestlen){
      $tmp = $cnt;
      $bestlen = $len;
      }
    }
  if($tmp >= 0){
    # print "best hit: $acc -> hitnr $tmp -> $bestlen bp\n";
    push(@tokeep, $tmp);
    }
  }


# focus on these best hits, store them into %besthits in a format allowing them to be sorted easily
my %bests;
foreach $cnt (@tokeep){
  $start = $coords{$cnt}{"start"};
  $bests{$start}{$cnt}{"acc"} = $coords{$cnt}{"acc"};
  $bests{$start}{$cnt}{"start"} = $coords{$cnt}{"start"};
  $bests{$start}{$cnt}{"stop"} = $coords{$cnt}{"stop"};
  }

my @allstarts = keys(%bests);
@allstarts = sort{$a <=> $b} @allstarts;


print "Keeping best blast hits... done\n";



# Perform binning
# Second, visit %bests, following startpos order 
# define non-overlapping stacks; and keep longest member of each stack.
my $stacknr = 0; 
my $rightlim = 0; #rightmost extent of ongoing stack
my $leftlim = 0;
my %stacks;


# Attribute hits within stacks
foreach $start (@allstarts){ #visit CDS hash following increasing order of start pos
  foreach $cnt (keys(%{$bests{$start}})){
    my %tmp = %{$bests{$start}{$cnt}};
    my $acc = $tmp{"acc"};
    my $start = $tmp{"start"};
    my $stop = $tmp{"stop"};

    if($start <= $rightlim){ #startpos < rightmost limit of stack, we are still in ongoing stack
      ## keep track of stop position and update rightlim if needed
      if($stop >= $rightlim){ #start pos still in same stack, but stop pos implies to update rightmost limit
	$rightlim = $stop; #update rightmost limit of ongoing stack
	}

      ## store data
      $stacks{$stacknr}{$cnt}{"acc"} = $acc;
      $stacks{$stacknr}{$cnt}{"start"} = $start;
      $stacks{$stacknr}{$cnt}{"stop"} = $stop;
      $stacks{$stacknr}{$cnt}{"seq"} = $fasta{$acc};
            
      } else { #start pos > rightmost limit, we enter in a new stack   
      ## open new stack    
      $stacknr++; #update stack number
      $track{$stacknr}{"start"} = $start;
      $rightlim = $stop; #update rightmost limit of ongoing stack

      ## store data
      $stacks{$stacknr}{$cnt}{"acc"} = $acc;
      $stacks{$stacknr}{$cnt}{"start"} = $start;
      $stacks{$stacknr}{$cnt}{"stop"} = $stop;
      $stacks{$stacknr}{$cnt}{"seq"} = $fasta{$acc};
      }
    }
  }
my $stackmax = $stacknr;



### make some stats about obtained stacks, and find their limits
my %track;
open(OUT, ">$outfolder/ChunckStats.txt");
open(OUT2, ">$outfolder/ChunckLimits.txt");
print OUT "Accession\tChunck\tRegion\tStart\tStop\tLength\n";
    
foreach $stacknr (1..$stackmax){
  my $left = 1e500;
  my $right = 0;
  my $nacc = 0;
  foreach $cnt (keys(%{$stacks{$stacknr}})){
    $acc = $stacks{$stacknr}{$cnt}{"acc"};
    my $start = $stacks{$stacknr}{$cnt}{"start"};
    my $stop = $stacks{$stacknr}{$cnt}{"stop"};
          
    if($start < $left){
      $left = $start
      }
    
    if($stop > $right){
      $right = $stop
      }
    my $range = $stop - $start;
    
    my $chunck = $acc;
    $acc =~ s/\.chk.*$//;
    $nacc++;
    
    print OUT "$acc\t$chunck\t$stacknr\t$start\t$stop\t$range\n";
    }
    
  my $chunckrange = $right-$left;
  $nacc = $nacc+1; #we started at 0 -> need to correct counts.
  print "Chunck $stacknr:\tfrom $left to $right bp\t>> width $chunckrange bp\t>> contains $nacc\tsequences\n";
  print OUT2 "Chunck $stacknr:\tfrom $left to $right bp\t>> width $chunckrange bp\t>> contains $nacc\tsequences\n";
  $track{$stacknr}{"start"} = $left;
  $track{$stacknr}{"stop"} = $right;
  }
  
close(OUT);
close(OUT2);



foreach $stacknr (1..$stackmax){
  open(OUT, ">$outfolder/chunck$stacknr.fas");
  
  # ref seq
  my $start = $track{$stacknr}{"start"};
  my $stop = $track{$stacknr}{"stop"};
  
  my $refseq = $ref{$refname};
  my @tmp = split("", $refseq);
  $refseq = join("", @tmp[$start..$stop]);
  
  print OUT ">$refname\n$refseq\n";
  

  # actual refs
  foreach $cnt (keys(%{$stacks{$stacknr}})){
    $acc = $stacks{$stacknr}{$cnt}{"acc"};
    $seq = $stacks{$stacknr}{$cnt}{"seq"};    
    print OUT ">$acc\n$seq\n";
    }
    
  close(OUT);
  }


print "done\n";
