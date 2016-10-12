#!/bin/perl

# Usage:
# perl CreateLocusReference_loader.pl MaxUse
#
# MaxUse = Maximal load allowed on the computing server.
#
# expects input in data.in:
# 	- expects collated.fa, the contigs obtained for each Rad locus (one sequence per locus and per specimen)
# 	- expects ForwardReads/*.fq_1, the forward reads (those having the restriction site) of each locus / specimen. 
#		File names must follow nomenclature: Rad[NrOfLocus].[SpecimenID].fq_1; e.g. Rad1027.PYD1.fq_1
#		Typically obtained from PreparePairsVelvet.pl (from STACKs pipeline)
#
# expects folders
#   data.in/ contains SNPs from STACKs post-clean script
#   bin/ contains scripts toolbox
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

my $rootdir = getcwd;


#### Get script arguments
my $scriptname = "RunSAMTOOLS";
my $read = $ARGV[0];
my $ref = $ARGV[1];
chomp($read);
chomp($ref);


# get basenames
my $bsn = basename($read);
$bsn =~ s/\.fq\_1//g;
$bsn =~ s/\.fq\_2//g;


# produce sai files
my $command = "bwa aln -t 4 -n 0.2 -f tmp/sai/$read.sai $ref reads/$read 2> /dev/null";
# print "### $scriptname : $command\n";
system("$command");

if($read =~ /\.fq\_2$/){
  # produce sam files
  my $command = "bwa sampe -f tmp/sam/$bsn.sam $ref tmp/sai/$bsn.fq_1.sai tmp/sai/$bsn.fq_2.sai reads/$bsn.fq_1 reads/$bsn.fq_2 2> /dev/null";
#   print "### $scriptname : $command\n";
  system("$command");

  # sort sam files
  my $command = "samtools view -bhS tmp/sam/$bsn.sam | samtools sort - tmp/sam/$bsn 2> /dev/null";
#   print "### $scriptname : $command\n";
  system("$command");

  # clean mess
  $command = "mv tmp/sam/$bsn.bam tmp/bam";
#   print "### $scriptname : $command\n";
  system("$command");
  }