#!/bin/perl

### Blasts and scaffolds mtDNA against a ref. Deals with circular DNA issues.
### WARNING: will not check that contigs are not overlapping at ends
### Multithreaded
# Usage:
# perl RunDCMegaBlast.pl MaxUse
#
# MaxUse = Maximal load allowed on the computing server before starting new job.
#
### Performed steps:
# 1. Blasts all contig files provided in data.in against any ref provided in ref/ (multithreaded, each blastn takes 4 cores)
# 2. Parse blast hits to fetch regions in contigs that correspond to ref
# 3. Rev-Comp these regions if needed,
# 4. Concatenates matching regions into a "pseudo-scaffold", following synteny given by reference. 
#    Enforces that pseudo-scaffolds start at position 1 (dealing with circular DNA here)
#    Only regions that produce acceptable blast hits end-up in the output.
#
# An example to make it clearer:
# contig : ATTAGTGCTGATCGTAGATTAGGATGCTTAGAGTGC
# ref    : AAAAAAGCTGATCGTAGATCAGGATGCTTAAAAAAA
# Match  :       ************* *********
# Region :       <--------------------->
# Length = 23 bp
# %simil = 22/23 = 96%
#
# Only the region <---------------------> is passed to the output, based on length and % similarity filters (see hard-coded params below).
# 
# expects inputs:
#	- Reference sequences of each marker: to be stored in refs/*.fna
#	- Contig files (one fasta per specimen to analyse) in contigs/
#
# Results in data.out/AllMito.fas
#
# Nils Arrigo, Unil 2015
####
use threads;
use threads::shared;
use File::Basename;


#### Get script arguments
my $MaxUse = $ARGV[0];
my $scriptname = "RunDCBlast";

#### Hard-coded params
my $minpct = 75 # minimum % similarity between matching region and reference
my $minlen = 500 # minimum length (bp) of region matching between considered contig and reference
my $mincov = 1 #number of reads supporting the considered contig
my $bckbonelen = 20000 #length of template sequence (N's) that will be filled with actual data. Set to something longer than a mitonchondria (e.g. 20000)




##############################################################
# start fresh and prepare projects folder
my $command = "rm -rf tmp data.out refs/*.fa*.*";
print "### $scriptname : $command\n";
system("$command");

my $command = "mkdir -p data.out/ tmp/ tmp/blast tmp/fasta tmp/aln";
print "### $scriptname : $command\n";
system("$command");


# Prepare BLAST database 
my @refs = `ls refs/*.f*`;
my $ref = $refs[0];
chomp($ref);
system("makeblastdb -dbtype nucl -in $ref");


##############################################################
### Contigs: list contigs to analyse
my @ctgs = `ls contigs/*.fas`;

### Contigs: run blast searches
my @threads;
foreach $contig (@ctgs){
  $contig = basename($contig);
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + $CPU;
  
  if($load + 4 < $MaxUse) {
    chomp($contig);

    #### Automated loader
    my $t = threads->new(\&runblast, $contig, $ref);
    push(@threads, $t);
    sleep(20)
    } else {
    redo;
    sleep(20)
    }
  }
foreach (@threads) {
  my $genome = $_-> join;
  }


### Collect Assemblies, filter out bad contigs
$outfile = "AllMito.fas";

open(OUT, ">$outfile");
close(OUT);

foreach $contig (@ctgs){
  $contig = basename($contig);
  chomp($contig);
#   my $command = "perl bin/ParseBlast_append.pl tmp/blast/$i.contig.blast contigs/$i.contig tmp/fasta/$i.fasta";
#   my $command = "perl bin/ParseBlast_append2.pl tmp/blast/$contig.blast contigs/$contig data.out/AllMito.fas 0";
  my $command = "perl bin/PseudoAlign.pl tmp/blast/$contig.blast contigs/$contig data.out/AllMito.fas $minpct $minlen 1 $mincov $bckbonelen";
  print("$command\n");
  system($command);
  }


### Contigs: run alignments against reference
my @ctgs = `ls tmp/fasta/*.fasta`;
my @threads;
foreach $contig (@ctgs){
  $contig = basename($contig);
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + $CPU;
  
  if($load + 4 < $MaxUse) {
    chomp($contig);

    #### Automated loader
    my $t = threads->new(\&runmuscle, $contig, $ref);
    push(@threads, $t);
    sleep(3);

    } else {
    sleep(3);
    redo;
    }
  }
foreach (@threads) {
  my $genome = $_-> join;
  }


### Produce final consensus -> final sequence



###### ROUTINES
## Blast search
sub runblast {
  my ($read, $ref) = ($_[0], $_[1]);
  my $command = "blastn -task dc-megablast -db $ref -query contigs/$contig -outfmt \"6 qacc sacc pident length mismatch gaps qstart qend sstart send evalue bitscore sstrand\" -num_threads 1 -evalue 0.1 -out tmp/blast/$contig.blast";
#   my $command = "blastn -task dc-megablast -db $ref -query contigs/$contig -outfmt 6 -num_threads 4 -evalue 0.1 -out tmp/blast/$contig.blast";
  print "$command\n";
  system($command);
  }

## Muscle
sub runmuscle {
  my ($contig, $ref) = ($_[0], $_[1]);

  my $command = "muscle -profile -in1 $ref -in2 tmp/fasta/$contig -out tmp/aln/$contig.aln";
  print "$command\n";
  system($command);
  }






