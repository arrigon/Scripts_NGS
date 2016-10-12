#!/bin/perl

# Usage (internal use only)
# perl RunSamTools read ref
#
# MaxUse = Maximal load allowed on the computing server.
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
$bsncomp = $bsn;
$bsn =~ s/\.fq\_1//g;
$bsn =~ s/\.fq\_2//g;

## index references
my $command = "bwa index $ref";
print "### $scriptname : $command\n";
system("$command");

# produce sai files
my $command = "bwa aln -n 10 -f tmp/sai/$bsncomp.sai $ref $read 2> /dev/null";
print "### $scriptname : $command\n";
system("$command");

if($read =~ /\.fq\_2$/){
  my $altread = $read;
  $altread =~ s/\.fq\_2/\.fq\_1/g;

  # produce sam files
  my $command = "bwa sampe -f tmp/sam/$bsn.sam $ref tmp/sai/$bsn.fq_1.sai tmp/sai/$bsn.fq_2.sai $altread $read 2> /dev/null";
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