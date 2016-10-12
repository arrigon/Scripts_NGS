#!/bin/perl

##### Script for aligning DMPS datasets
### Parameters
# Usage: perl Pipeline_Aligner_std_queue.pl MaxUse
#
# MaxUse = max CPU usage
#
# Input: assumes that contigs to be aligned (ideally, a DMPS dataset) are in data.in/ (as a single fasta file)
# N.B. We need a reference to align against.
#
# Output: in data.out
#
# Performed steps:
# -- Step 1. Prepare input fasta; slincing it at gaps positions
# -- Step 2. Make blastn search between this fasta and ref sequence
# -- Step 3. Binning of blast hits into non-overlapping regions (usually corresponds to DMPS sites that were successfully amplified)
# -- Step 4. Align each region individually (we go parallel here)
# -- Step 5. Refine each chunck alignment
# -- Step 6. Collapse all chuncks back into a global alignment
# -- Step >= 7. Remove remaining assembly dusts, two passes; Remove empty and degenerated sites; produce output
#########################################################
use File::Basename;
use threads;
use threads::shared;



### Get script arguments
my $maxuse = $ARGV[0];
my $scriptname = "Aligner_Pipeline";



### Prepare workfolders
my $command = "rm -rf tmp/ data.out/";
print("\n$scriptname running:\n$command\n\n");
system($command);


my $command = "mkdir data.out tmp tmp/fasgapsliced tmp/blast tmp/alnI tmp/alnII tmp/alnIII";
print("\n$scriptname running:\n$command\n\n");
system($command);



### Prepare reference for blast searches
my @list = `ls ref/*.fa*`;
my $ref = $list[0];
chomp($ref);
$ref = basename($ref);

my $command = "makeblastdb -dbtype nucl -in ref/$ref";
print("\n$scriptname running:\n$command\n\n");
system($command);



### Step 1. Prepare input fasta; slincing it at gaps positions
my @list = `ls data.in/*.fa*`;
my $infile = $list[0];
chomp($infile);
$infile = basename($infile);
my $bsn = $infile;
$bsn =~ s/\..*$//;

my $outfile = "data.out/$bsn.final.fas";

my $command = "perl bin/Fas_SliceGapsMaxLen.pl data.in/$infile tmp/fasgapsliced/$infile.sliced.fas 30 8000";
print("\n$scriptname running:\n$command\n\n");
system($command);



### Step 2. Make blastn search between this fasta and ref sequence
my $command = "blastn -task blastn -db ref/$ref -query tmp/fasgapsliced/$infile.sliced.fas -outfmt 6 -num_threads 1 -evalue 0.01 -perc_identity 10 -out tmp/blast/$infile.sliced.blast";
print("\n$scriptname running:\n$command\n\n");
system($command);



### Step 3. Bin blast hits into non-overlapping regions (usually corresponds to DMPS sites that were successfully amplified)
my $command = "perl bin/Blast_DefineBins.pl tmp/fasgapsliced/$infile.sliced.fas ref/$ref tmp/blast/$infile.sliced.blast tmp/chuncks";
print("\n$scriptname running:\n$command\n\n");
system($command);



### Step 4. Align each chunck individually (we go parallel here)
# Get list of chuncks to align
my @chuncks = `ls tmp/chuncks/*.fas`;

### Loop over each chunck
my @threads;

foreach my $i (@chuncks){
  chomp($i); 
  
  # run params
  my $infile = "$i";

  my $bsn = basename($i);
  $bsn =~ s/\.fa.*$//;
  
  my $outfile = "tmp/alnI/$bsn.aln";
  my $workfolder = "tmp/alnI";


  # Get CPU load
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + 1;

  if($load < $maxuse) {
    my $t = threads->new(\&runaln, $infile, $outfile, $workfolder);
    push(@threads, $t);
    sleep(1);

    } else {
    sleep(2);
    redo;
    }
  }

### Survey ongoing threads
foreach (@threads) {
  my $infile = $_-> join;
  print "done with $infile\n";
  }

  
  
### Step 5. Refine each chunck alignment
my @chuncks = `ls tmp/alnI/*.aln`;

foreach my $i (@chuncks){
  chomp($i); 
  
  my $bsn = basename($i);
  $bsn =~ s/\.aln$//;
  
  my $command = "perl bin/ALN_ConsenseChuncks.pl tmp/alnI/$bsn.aln tmp/alnII/$bsn.aln";
  print("\n$scriptname running:\n$command\n\n");
  system($command);
  }

  

### Step 6. Collapse all chuncks back into a global alignment
my $command = "R CMD BATCH '--args infolder=\"./tmp/alnII/\" outfile=\"tmp/alnIII/$bsn.fas\"' bin/ALN_MergeAlignments2FAS.r bin/R.out";
print("\n$scriptname running:\n$command\n\n");
system($command);



## Step 8. Remove remaining assembly dusts, two passes
my $command = "perl bin/ALN_DustMasking_queue.pl tmp/alnIII/$bsn.fas tmp/alnIII/$bsn.dustsI.fas 5 20 $maxuse";
print("\n$scriptname running:\n$command\n\n");
system($command);

my $command = "perl bin/ALN_DustMasking_queue.pl tmp/alnIII/$bsn.dustsI.fas tmp/alnIII/$bsn.dustsII.fas 2 5 $maxuse";
print("\n$scriptname running:\n$command\n\n");
system($command);



## Step 9. Remove empty and degenerated sites
my $command = "R CMD BATCH '--args infile=\"tmp/alnIII/$bsn.dustsII.fas\" outfile=\"tmp/alnIII/$bsn.dustsIII.fas\"' bin/ALN_CleanAlignments.r bin/Clean.out";
print("\n$scriptname running:\n$command\n\n");
system($command);


## Step 10. Remove remaining assembly dusts
my $command = "perl bin/ALN_DustMasking_queue.pl tmp/alnIII/$bsn.dustsIII.fas $outfile 5 3 $maxuse";
print("\n$scriptname running:\n$command\n\n");
system($command);


### ROUTINES
sub runaln {
  my ($infile, $outfile, $workfolder) = ($_[0], $_[1], $_[2]);

  my $command = "perl bin/ALN_DivideAndConquer.pl $infile $outfile $workfolder";
  print("\n$scriptname running:\n$command\n\n");
  system($command);
  
  return($infile);
  }
  
