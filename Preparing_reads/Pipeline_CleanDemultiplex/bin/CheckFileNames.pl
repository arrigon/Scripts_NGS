#!/bin/perl

##### Renaming fastq files, BEFORE starting the demultiplex / cleaning stuff
# 
# Input: $infolder = folder where raw fastq files are all stored
#        $outfolder = folder where renamed fastq files are saved
#
# Usage: perl CheckFileNames.pl infolder outfolder
#
# WARNING: this script is based on hard-coded RegExp, see comments below
########################################################


### Script params
my $infolder = $ARGV[0];
my $outfolder = $ARGV[1];


### prepare outfolder
$command = "mkdir $outfolder";
print("$command\n");
system($command);

 
### list fastq files into infolder
my @list = `ls $infolder`;

foreach my $file (@list){ #loop over all fastq files to handle
  chomp($file);

  ### Extract infos from current fastq file name
  # this is where you MUST adapt the following RegExp according to your needs
  # For instance, the typical CIG file names are XXX_YYY_R1_001.fastq.gz
  # where XXX is the lane ID, YYY is a random character string, R1 (or R2) is the read orientation and 001 is the specimen / lane chunck ID
  # All we need to keep is XXX, R1 and 001
  # To do so (and from the mentionned exemple), we can use the following RegExp : /(.*)_.*_(R.*)_(.*)\.fastq\.gz/ 
  # and extract XXX with $1
  # and extract R1 with $2
  # and extract 001 with $3
  ## Yet, not all files will look the same way and you need to adapt your regexp accordingly

  ### WARNING REGEXP: to be adapted here
  # $file =~ /(.*-)(.*)_(R.*)\.fastq\.gz/; #TGAC
  $file =~ /(.*)_.*_.*_(R.*)_(.*)\.fastq\.gz/; #Tom's lanes XXX_L5_YYY_R1_001.fastq.gz, we get ridd of L5 and YYY


  ### Extracting infos from raw file name, using RegExp mapping (make sure that $1, $2 and $3 are mapping to the right thing)
  my $tech = $1; # XXX, is laneID (to extract from random Illumina filenames), must be shared among all files of that lane
  my $readnr = $2; # R1 or R2
  my $specnr = $3; # can be either specimen or fastq chunck number (i.e. lane chuncks)

 
  ### Remove any weird character from that laneID
  $tech =~ s/[\.|_|-]//g;


  ### Generate the new file name: XXX.R1.001.fastq.gz (this is the "checked" namestyle)
  my $newname = "$tech.$readnr.$specnr.fastq.gz"; #name looks like laneID.R1.001.fastq.gz


  ### Copy-paste-rename the files to a new folder
  $command = "cp $infolder/$file $outfolder/$newname";
  print("$command\n");
  system($command);
  }
