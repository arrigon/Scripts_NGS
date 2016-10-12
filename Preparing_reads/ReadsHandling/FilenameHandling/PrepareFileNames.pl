#!/bin/perl

##### Renaming fastq files, BEFORE starting the demultiplex / cleaning stuff
# 
# Input: $infolder = folder where raw fastq files are all stored
#
# Usage: perl PrepareFileNames.pl infolder
#
# Output: infolder/Original = ancient files
# 	  infolder/OK = renamed files; ready to be demultiplexed
#
# WARNING: this script is based on hard-coded RegExp, see comments below at line 48
########################################################
use File::Basename;

### Script params
my $infolder = $ARGV[0]; #folder that contains infiles


### Prepare outfolders
my $origfolder = "$infolder/Original/";
my $okfolder = "$infolder/OK/";
system("mkdir $origfolder");
system("mkdir $okfolder");
 
# list fastq files into infolder
my @list = `ls $infolder/*.fastq.gz`;

foreach my $file (@list){ #loop over all fastq files to handle
  chomp($file);
  $file = basename($file);

  # move to orig folder (using links to preserve disk space)
  my $command = "mv $infolder/$file $origfolder/$file";
  print("$command\n");
  system($command);
  }


### 2. Check file names consistency
# list fastq files into infolder
my @list = `ls $origfolder`;


foreach my $file (@list){ #loop over all fastq files to handle
  chomp($file);
  $file = basename($file);
  
  ### WARNING REGEXP: to be adapted here
  # Typically, the Illumina libraries come with varying naming schemes. These names contain three important informations:
  # - laneID
  # - read orientation
  # - fastq number
  ### We need to know these things before demultiplexing the files; notably because we have to keep read pairs synchronized.
  ### Most of the time, the file names look something like 
  # EPH1a_lfpT4Ix0oppa_L1_R1_001.fastq.gz, where 
  # "EPH1a" = laneID
  # "R1" = read orientation
  # "001" = fastq number
  # "lfpT4Ix0oppa" is some random junk; used to prevent files erasing the one another
  # "L1" is the index number (when using several indexes)
  ### Annoyingly, these items do not always come in the order; and we have to check manually what we are dealing with.
  ### This is where we need this script for renaming our fastq files; in a coherent way, so that the demultiplexing pipeline
  ### knows how to deal with it. 
    
  ### To this end, we use a regular expression that looks for the needed items.
  ### You must pick the RegExp that suits your needs here
  ### The syntax goes as follows (check http://www.troubleshooters.com/codecorn/littperl/perlreg.htm, paragraph Doing String Selections (Parsing))
  # $file =~ /(.*)_.*_.*_(R.*)_(.*)\.fastq\.gz/; #deals with EPH1a_lfpT4Ix0oppa_L1_R2_001.fastq.gz
  # EPH1a_lfpT4Ix0oppa_L1_R2_001.fastq.gz
  # (.*) _     .*     _.*_(R.*)_(.*).fastq.gz
  # every item matching with () can then be retrieved using $1, $2, $3
  
  ###
  # $file =~ /(.*)_.*_.*_(R.*)_(.*)\.fastq\.gz/; #MasRAD01_84S8HDfA0iWj_L7_R1_002.fastq.gz, we get ridd of YYY and L5
  $file =~ /(.*)_.*_.*_(R.*)_(.*)\.fastq\.gz/; #EPH1a_lfpT4Ix0oppa_L1_R2_001.fastq.gz
					       #LycCa01_89Ngb5BsvQFP_L1_R1_002.fastq.gz


  ### Extracting infos from raw file name, using RegExp mapping (make sure that $1, $2 and $3 are mapping to the right thing)
  my $tech = $1; # XXX, is laneID (to extract from random Illumina filenames), must be shared among all files of that lane
  my $readnr = $2; # R1 or R2
  my $specnr = $3; # can be either specimen or fastq chunck number (i.e. lane chuncks)

 
  ### Remove any weird character from that laneID
  $tech =~ s/[\.|_|-]//g;


  ### Generate the new file name: XXX.R1.001.fastq.gz (this is the "checked" namestyle)
  my $newname = "$tech.$readnr.$specnr.fastq.gz"; #name looks like laneID.R1.001.fastq.gz


  ### Copy-paste-rename the files to a new folder, use links to save disk space
  my $command = "cp -l $origfolder/$file $okfolder/$newname";
  print("$command\n");
  system($command);
  }

