#!/bin/perl

##### Fakes a second set of fastq files (SE -> PE)
# 
# Input: $infolder = folder where raw fastq files are all stored
#
# Usage: perl FakePEdata.pl infolder
#
########################################################


### Script params
my $infolder = $ARGV[0];
my $target = "\_R1\_";
my $dest = "\_R2\_";

 
### list fastq files into infolder
my @list = `ls $infolder`;

foreach my $file (@list){ #loop over all fastq files to handle
  chomp($file);

  # rename file
  my $newfile = $file;
  $newfile =~ s/$target/$dest/;

  # link to new copy (using links to preserve disk space)
  $command = "cp -l $infolder$file $infolder$newfile";
  print("$command\n");
  system($command);
  }
