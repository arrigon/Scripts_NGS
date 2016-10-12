#!/bin/perl

# Usage:
# perl CallSNPs_loader_queue2.pl Outfile MaxUse
#
# Outfile = name of global VCF file
# MaxUse = Maximal load allowed on the computing server.
#
### Performed steps:
# 1. Align reads against reference (specimen by specimen, and marker by marker)
# 2. Calls SNPs from these alignments, Do this marker by marker
# 3. Merge all SNPs into a global VCF matrix => first output
# 4. Imputes missing data and performs haplotype reconstruction (using reads) => second output
#
# expects inputs:
#	- Reference sequences of each marker: to be stored in refs/AllRefs.fa
#	- Reads (forward and reverse, *.fq_1 and *.fq_2) to be mapped against these refs: to be stored in reads/ 
#	  N.B: One pair of fastq per specimen and per marker.
#	  File names must follow nomenclature: Rad[NrOfLocus].[SpecimenID].fq_1; e.g. Rad1027.PYD1.fq_1
#
# Results are in folder data.out/
#
# TODO: implement subfolders having names based on input file, in a way to enable parallel runs
# TODO: NCPU should govern how many searches can be launched per input file
#
# Nils Arrigo, Uni of Arizona 2012
# use Cwd;
####
use threads;
use threads::shared;
use File::Basename;


#### Get script arguments
my $MaxUse = $ARGV[0];
my $scriptname = "BWA\_aligner";



##############################################################
## start fresh and prepare projects folder
my $command = "rm -rf tmp data.out";
print "### $scriptname : $command\n";
system("$command");

my $command = "mkdir -p data.out/ tmp/ tmp/sam tmp/bam tmp/sai";
print "### $scriptname : $command\n";
system("$command");


## check contents of refs folder
my @refs = `ls refs/*.fa*`;
my $ref = @refs[0];
chomp($ref);

## index references
my $command = "bwa index $ref";
print "### $scriptname : $command\n";
system("$command");



##############################################################
### BAM files: Align reads against references
# Check contents of reads folder and iterate over fastq files
my @reads = `ls reads/`;
foreach $read (@reads){
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + $CPU;
  
  if($load < $MaxUse) {
    chomp($read);

    #### Automated loader
    $command = "perl bin/RunSamtools.pl $read $ref &";
    system("$command");
    sleep(3);

    } else {
    sleep(3);
    redo;
    }
  }

# clean mess
$command = "rm tmp/sai/*";
print "### $scriptname : $command\n";
system("$command");