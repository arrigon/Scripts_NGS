#!/bin/perl

### Performs a standard reads mapping against ref. for Single End libraries
# Usage: perl RunSamtools_SE.pl readfile ref
#
# Outputs directed to tmp/*
#
# Nils Arrigo, Unil 2014
# use Cwd;

use threads;
use threads::shared;
use File::Basename;

my $rootdir = getcwd;


#### Get script arguments
my $scriptname = "RunSAMTOOLS_SE";
my $read = $ARGV[0];
my $ref = $ARGV[1];
my $mism = $ARGV[2];
chomp($read);
chomp($ref);


# get basenames
my $bsn = basename($read);
$bsn =~ s/\.fq\_1//g;

# produce sai files
my $command = "bwa aln -t 4 -n $mism -f tmp/sai/$read.sai $ref reads/$read 2> /dev/null";
# print "### $scriptname : $command\n";
system("$command");


# produce sam files bwa samse ref.fa aln_sa.sai short_read.fq > aln-se.sam 

my $command = "bwa samse $ref tmp/sai/$bsn.fq_1.sai reads/$bsn.fq_1 > tmp/sam/$bsn.sam";
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
