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
my $scriptname = "RunfastPHASE";
my $workdir = $ARGV[0];
my $input = $ARGV[1];
chomp($workdir);
chomp($input);


chdir($workdir);

# run fastPHASE
my $command = "fastPHASE_Linux $input -KL 2 -KU 8 -KI 1 -T 50 2> /dev/null";
print "### $scriptname : $command\n";
system("$command");
