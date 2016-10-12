#!/bin/perl

##### Barcode trimming, Restriction site filtering, Uniform length trimming
# 
# Input: $infolder = folder where input fastq (*.fq_1 and *.fq_2) are all stored
# ASSUMES THAT *.fq_1 CONTAINS BARCODE AND RESTRICTION SITE
#
# One must provide:
# - BarcodeLength = length of barcodes, will be trimmed from sequence start
# - Motif = Restriction site (e.g. AATTC), will work using perfect matches, put "keepall" to bypass that filter
# - MaxLength = STACKs needs reads of uniform length, specify it here. Use MaxLength < 0 to keep full length
# - Keep barcode (0 = No, 1 = yes)
# 
# The script will apply the cleaning procedure (bin/RestrictionSiteFilterFastqII.pl) on all pairs of fastq files
# The resulting fastq will be kept in sync.
# 
# Outputs produced in $outfolder
#
# Usage: perl FinalClean.pl infolder outfolder Motif BarcodeLength MaxLength KeepBarcode
########################################################

### Get script arguments
my $infolder = $ARGV[0];
my $outfolder = $ARGV[1];
my $motif = $ARGV[2];
my $bcdelgth = $ARGV[3];
my $limlen = $ARGV[4];
my $keepbc = $ARGV[5];
my $trimR2 = $ARGV[6];
my $stitch = $ARGV[7];
chomp($stitch);



### Final Cleaning Step
# Get files in $outmultiplex
my @list = `ls $infolder`;

# reconstruct pairs of files
my %corresp;
my %chuncks;
foreach my $i (@list){
  chomp($i);
  if($i =~ /(.*)\.fq\_(.*)/){
    $corresp{$1}{$2} = $i;
    } 
  }

# loop over them, pair by pair
foreach my $i (keys(%corresp)){
  
  # get file names
  my $file1 = $corresp{$i}{"1"};
  my $file2 = $corresp{$i}{"2"};

  # repair, clean and sync fastq files (demultiplex introduces weird things)
  my$command = "perl bin/RestrictionSiteFilterFastqVI.pl $infolder/$file1 $infolder/$file2 $motif $bcdelgth $limlen $outfolder $keepbc $trimR2 $stitch";
  print "Running: $command\n\n";
  system($command);  
  } 
