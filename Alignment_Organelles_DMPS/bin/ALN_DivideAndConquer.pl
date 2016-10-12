#!/bin/perl

#####################################
#### We assume that first seq IS reference
####
####################################
#### usage perl ALN_DivideAndConquer.pl infile workfolder outfile
####
#### Nils Arrigo, UNIL 2014
#####################################

#load packages
use File::Basename;
use POSIX;



# Parameters
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $workfolder = $ARGV[2];


# Prepare workfolder
my $bsn = basename($infile);
$bsn =~ s/\.fa.*$//;
$workfolder = "$workfolder/DAC.$bsn";

system("mkdir $workfolder/ $workfolder/fas $workfolder/aln $workfolder/profiles");



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
  $fasta{$tmp[0]} = $seq; 
  push(@accs, $tmp[0]);
  }
my @originalaccs = @accs;
my $refname = $accs[0];
my $refseq = $fasta{$refname};



### Perform pairwise alnmts, (i.e. focal seq versus refseq)
print "Individual alignments against reference...\r";

my $cnt = 0;
foreach my $acc (@accs){
  open(OUT, ">$workfolder/fas/$cnt.fas");
  my $seq = $fasta{$acc};
  print OUT ">$refname\n$refseq\n>$acc\n$seq\n";
  close(OUT);
  
  system("muscle -quiet -in $workfolder/fas/$cnt.fas -out $workfolder/aln/$cnt.aln"); #muscle works better
#   system("./bin/clustalo --force --in=$workfolder/fas/$cnt.fas --out=$workfolder/aln/$cnt.aln");
  $cnt++;
  }
my $maxcnt = $cnt-1;
print "Individual alignments against reference... done\n";
 
  
  
### Reconcile all pairwise alnmts using profile alnmt
print "Merging individual alignments... \r";

# initiate merging with first chunck
system("cp $workfolder/aln/1.aln $workfolder/profiles/tmp.aln");

# perform profile alignments, iteratively to first chunck
foreach my $cnt (2..$maxcnt){
  system("./bin/clustalo --force --profile1=$workfolder/profiles/tmp.aln --profile2=$workfolder/aln/$cnt.aln --outfile=$workfolder/profiles/tmp.aln");
  $cnt++;
  }
print "Merging individual alignments... done\n";
  
  
  
### Perform final cleaning of alnmt  
# load fasta into hash
open(FILE, "$workfolder/profiles/tmp.aln");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

my %aln;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = $_;
  $aln{$tmp[0]} = $seq; 
  }


open(OUT, ">$outfile");
foreach $acc (@originalaccs){
  my $seq = $aln{$acc};
  print OUT ">$acc\n$seq\n";
  }
close(OUT);



# clean mess
system("rm -rf $workfolder/");
