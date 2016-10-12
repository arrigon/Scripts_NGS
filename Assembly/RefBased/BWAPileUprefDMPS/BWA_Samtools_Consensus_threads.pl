#!/bin/perl

#### Reconstructing an alignment, using ref-based read mapping.
# This script relies on a read aligner (BWA) and samtools (mpileup program) 
# to build a consensus sequence from read mappings.
# Works using parallel threads, to speed up data processing.
#
# WARNING: some key parameters are hard coded just after this help. Look for them and modify accordingly.
#
#
# Usage:
# perl BWA_Samtools_Consensus_threads.pl MaxUse
#
# MaxUse = Maximal load allowed on the computing server.
#
#
### Performed steps:
# 1. Align reads against reference (specimen by specimen, using BWA) 
#     -> you NEED specify the number of allowed mismatches to reference (see hard coded params below)
#     -> you NEED specify the data type (Paired-end vs Single-end, see hard coded params below)
# 2. Builds the consensus using samtools
# 3. Pushes back the alignment into data.out/
#     -> Checks coverage supporting each bp, you NEED specify the number minimum coverage (see hard coded params below)
#
#
# expects inputs:
#	- Reference sequence: to be stored in refs/*.fas #Put ONE reference at a time (e.g. Whole mitochondria / Chloroplast)
#	- Reads: to be stored in reads/ WARNING: names MUST finish with *.fq_1 / *.fq_2
#	  Paired-end data: (forward and reverse, *.fq_1 and *.fq_2, oriented as --> <--, if working with mate-pair libs, you'll need to rev-comp *.fq_2 before starting)
#	  N.B: One pair of fastq per specimen and per marker.
#	  File names must follow nomenclature: [SpecimenID].fq_1|2 e.g. PYD1.fq_1 and PYD1.fq_2
#	
#	  Single-end data: just one file per specimen.
#	  File names must follow nomenclature: [SpecimenID].fq_1 e.g. PYD1.fq_1
#
# Results are in folder data.out/
#
#
####
# Nils Arrigo, Unil September 2014
####

#### WARNING: Hard coded params
my $mismatches = 5; #number of allowed mismatches during read mapping 
		    #(absolute number, if reads are 100bp long, this would allow $mismatches/100)
		    #WARNING: setting this value is tricky, on the one hand, you want to avoid including many paralogs in a given assembly,
		    #especially when dealing with highly repeated regions (e.g. ITS, NUMTs in mitochondria). 
		    #On the other hand, if you work with a distantly related reference, 
		    #you'll be forced to authorize some divergence degree (otherwise no reads from the analyzed specimen will 
		    #map on your reference. Thus, be aware of the specifics of your study system.
my $mincov = 3; #minimum coverage to support a given base in consensus
my $paired = 1; #Are you using paired-end data ($paired = 1) or single-end ($paired = 0)?

#### Get script arguments
my $MaxUse = $ARGV[0];
my $scriptname = "RefGuidedConsensus_loader";
my $outfile = "RefBasedAssembly\_cov$mincov\_mism$mismatches.fasta";


##############################################################
#### Call needed libraries
use threads;
use threads::shared;
use File::Basename;

# ## start fresh and prepare projects folder
# my $command = "rm -rf tmp data.out";
# print "### $scriptname : $command\n";
# system("$command");

# my $command = "mkdir -p data.out/ tmp/ tmp/sam tmp/bam tmp/sai tmp/bcf tmp/vcf tmp/pileups";
# print "### $scriptname : $command\n";
# system("$command");
# 

## check contents of refs folder
my @refs = `ls ref/*.fa*`;
my $ref = @refs[0];
chomp($ref);

## index references
my $command = "bwa index $ref";
print "### $scriptname : $command\n";
system("$command");



##############################################################
### BAM files: Align reads against references
## Check contents of reads folder and iterate over fastq files
# my @reads = `ls reads/*.fq\_1`;
# 
# my @threads;
# my $startedthreads = 0;
# 
# foreach $read (@reads){
#   chomp($read);
# 
#   open PIPE, "uptime |";
#   my $line = <PIPE>;
#   close PIPE;
#   $line =~ s/\s//g;
#   my @lineArr =  split /:/, $line;
#   my $times = $lineArr[@lineArr-1];
#   my @timeArr = split /,/, $times;
#   my $load = $timeArr[0] + 4;
#   
#   if($load < $MaxUse) {
#     my $t = threads->new(\&runsamtools, $scriptname, $read, $ref, $mismatches, $paired);
#     push(@threads, $t);
#     sleep(20);
#     $startedthreads++;
# 
#     } else {
# 
#     sleep(20);
#     redo;
#     }
#   }
# my $cnt = 0;
# foreach (@threads) {
#   my $genome = $_-> join;
#   $cnt++;
#   print "Finished job $cnt / $startedthreads\n";
#   }
# 
# 
# # clean mess
# $command = "rm tmp/sai/* tmp/sam/*";
# print "### $scriptname : $command\n";
# system("$command");
# 


##############################################################
### Create consensus sequences
my @bams = `ls tmp/bam/*.bam`;

my @threads;
my $startedthreads = 0;
foreach my $bam (@bams){
  chomp($bam);
  my $bam = basename($bam);
  
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0];
  
  if($load < $MaxUse) {
    my $t = threads->new(\&runmpileup, $scriptname, $bam, $ref);
    push(@threads, $t);
    sleep(5);
    $startedthreads++;

    } else {

    sleep(5);
    redo;
    }
  }

my $cnt = 0;
foreach (@threads) {
  my $genome = $_-> join;
  $cnt++;
  print "Finished job $cnt / $startedthreads\n";
  }


### Final step: Collecting assemblies and checking coverage
my $command = "perl bin/CleanCollectFasta.pl tmp/pileups/$bam $mincov data.out/$outfile";
print("$command\n");
system($command);

# Collect assembly stats 
my $command = "R CMD BATCH \"--args infolder='tmp/pileups/' outfile='data.out/$outfile.pdf' \" bin/StatsCoverage.r";
print("\nRunning: $command\n\n");
system($command);


# # # clean mess
# $command = "rm -rf tmp/";
# print "### $scriptname : $command\n";
# system("$command");



sub runsamtools {
  my ($scriptname, $read, $ref, $mismatches, $paired) = ($_[0], $_[1], $_[2], $_[3], $_[4]);

  # run samtools
  if($paired == 1){
    $command = "perl bin/RunSamtools_PE.pl $read $ref $mismatches";
    print "### $scriptname : $command\n";
    system("$command");

    } else {
    $command = "perl bin/RunSamtools_SE.pl $read $ref $mismatches";
    print "### $scriptname : $command\n";
    system("$command");
    }
  }


sub runmpileup {
  my ($scriptname, $bam, $ref) = ($_[0], $_[1], $_[2]);

  my $bsn = $bam;
  $bsn =~ s/\.bam//;

  # run samtools mpileup WARNING: samtools mpileup needs -A "count anomalous read pairs" to be specified
  # this way, one can use ALL read pairs that were mapped on the ref. genome.
  # otherwise, samtools just discards any "weird" read pair, which ends-up discarding MOST of the data.
  my $command ="samtools mpileup -Auf $ref tmp/bam/$bam | bcftools view -cg - | vcfutils.pl vcf2fq | ./bin/seqtk seq -A - | sed -r \"s/\>.*/\>$bsn/g\" > tmp/pileups/$bsn.aln";
  print "### $scriptname : $command\n";
  system("$command");

  my $command = "samtools depth tmp/bam/$bam > tmp/pileups/$bsn.cov";
  print("$command\n");
  system($command);
  }

