#!/bin/perl

##### Quality checks
# 
# Input: $outfolder = folder where cleaned fastq are all stored
#
# The script will parse recursively all subdirectories of $outfolder
#  and will look for fastq files (using CollectSample.pl)
#
# Assumes that naming scheme of fastq files is as follow: SampleName.R1/R2.fastq (as produced by CollectSample.pl)
# 
# Outputs produced in data.out/$outfolder/final and data.out/$outfolder/logs
#
# Usage: perl Test.pl outfolder
########################################################

### Get script arguments
my $outfolder = $ARGV[0];
# $outfolder = data.out/$outfolder



### Collect files at the specimen level
my $command = "perl bin/CollectSample.pl $outfolder/tmp/qualclean $outfolder/final";
print("\nRunning: $command\n\n");
system($command);



### STEP 5 Make quality stats (run001 only)
system("cat $outfolder/tmp/qualclean/001/* > $outfolder/logs/allfiles.fastq");
system("./bin/FastQC/fastqc $outfolder/logs/allfiles.fastq");
system("rm $outfolder/logs/allfiles.fastq");



### STEP 6 Make quality stats (specimen level)
my @list = `ls $outfolder/final`;
my %corresp;

foreach my $i (@list){
  chomp($i);
  $i =~ /(.*)\.(.*)\.fastq/;
  $corresp{$1}{$2} = $i;
  }

foreach my $i (keys(%corresp)){
  my $file1 = $corresp{$i}{"R1"};
  my $file2 = $corresp{$i}{"R2"};
  system("cat $outfolder/final/$file1 $outfolder/final/$file2 > $outfolder/logs/complete.$i.fastq");
  system("./bin/FastQC/fastqc $outfolder/logs/complete.$i.fastq");
  system("rm $outfolder/logs/complete.$i.fastq");
  }

print "Done\n";



### STEP 7 Collect these stats 
my $command = "R CMD BATCH \"--args infolder='$outfolder/logs' outfolder='$outfolder' \" bin/QualStats.r";
print("\nRunning: $command\n\n");
system($command);

