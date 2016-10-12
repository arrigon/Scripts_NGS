#!/bin/perl

####################################
# Collecting fastq files of all samples, via recursive searching throughout all subdirectories of oufolder
# Assumes PE reads, keeps them in sync 
#
### Params
# - Takes fastq files (not fast.gz)
# - Scans all fastq files in infolder (which contains N subfolders)
#   Expects file names as "001/clean.ASPmef1.R1.001.fastq", i.e. runfolder/clean.SampleName.R1/R2.runfolder.fastq
#                          is used by RegExp /\/(\d*)\/clean\.(.*)\.(R.)\.(.*)\.fastq/
# - puts results into outfolder 
#
#
### Outputs: all files stored in following folder system
#		- final files: $outfolder/samples/$samplename
#
### Usage: perl CollectSample.pl infolder outfolder
#
# Nils Arrigo. Uni Lausanne 2012
############################################################

### Get script arguments
my $infolder = $ARGV[0];
chomp($infolder);

my $outfolder = $ARGV[1];
chomp($outfolder);



### Scan infolder to find all files
my @allfiles = `find $infolder -type f -name '*.fastq'`;



### Get specimen names from detected files
my %storage;
foreach my $path (@allfiles){
  chomp($path);
  $path =~ /\/(\d*)\/clean\.(.*)\.(R.)\.(.*)\.fastq/;
  my $foldernr = $1;
  my $samplename = $2;
  my $readnr = $3;
  
#   print "$foldernr\t$readnr\t$samplename\n";
  $storage{$samplename}{$foldernr}{$readnr} = $path;
  }



### Parse these guys, sequentially to produce a single R1 and R2 file per sample, in sync
my $cnt = 0;
my @allsamples = keys(%storage);
my $samplecount = $#allsamples + 1;

foreach my $samplename (@allsamples){
  my %files = %{$storage{$samplename}};
  
  open(R1File, ">$outfolder/$samplename.fq_1");
  open(R2File, ">$outfolder/$samplename.fq_2");

  foreach my $foldernr (keys(%files)){
    my $R1 = $storage{$samplename}{$foldernr}{"R1"};
    my $R2 = $storage{$samplename}{$foldernr}{"R2"};
    
    system("cat $R1 >> $outfolder/$samplename.fq_1");
    system("cat $R2 >> $outfolder/$samplename.fq_2");
    }
  
  $cnt++;
  
  print "Collected $cnt / $samplecount: specimen $samplename\n";

  close(R1File);
  close(R2File);
  }










