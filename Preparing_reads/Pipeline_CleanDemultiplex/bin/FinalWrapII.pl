#!/bin/perl

##### Final cleaning, trimming of barcodes, Uniformization of R1 reads (R2 keep stay the way they are) and Quality checks
# 
# Input: $workfolder = folder where cleaned fastq are all stored
#        $barcodes = path to barcodes (we need to know their length)
#        $MaxLen = Length of reads R1
#	 $Keepbc = keep barcode (0 = no, 1 = yes)
# The script will parse recursively all subdirectories of $outfolder
#  and will look for fastq files (using CollectSample.pl)
#
# Assumes that naming scheme of fastq files is as follow: SampleName.R1/R2.fastq (as produced by CollectSample.pl)
# 
# Outputs produced in data.out/$outfolder/final and data.out/$outfolder/logs
#
# Usage: perl FinalWrapII.pl WorkingFolder BarcodesPath MaxReadLength keepbc
########################################################

### Get script arguments
my $workfolder = $ARGV[0];
my $barcodes = $ARGV[1];
my $MaxLen = $ARGV[2];
my $keepbc = $ARGV[3];
my $motif = $ARGV[4];
my $trimR2 = $ARGV[5];
my $stitch = $ARGV[6];
chomp($stitch);

# make working dirs
system("mkdir $workfolder/tmp/merged $workfolder/stacks");

if($barcodes ne 0){ #make sure to use "ne", because $barcode comes as a character string.
  ### STEP 1 Collect files at the specimen level, merge them into $outfolder/final
  my $command = "perl bin/CollectSample.pl $workfolder/tmp/qualclean $workfolder/tmp/merged";
  print("\nRunning: $command\n\n");
  system($command);

  ## clean mess
  my $command = "rm -rf $workfolder/tmp/QualClean/*";
  print("\nRunning: $command\n\n");
  system($command);



  ### STEP 2 Make final cleaning: clip barcode, filter out reads starting with restriction site and make uniform read length 
  # forward reads only, reverse kept at varying lengths)

  # get barcode length
  open(IN, "$barcodes");
  my @line = <IN>;
  my @fields = split("\t", $line[1]);
  my $bclength = length($fields[1]) - 1;
  close(IN);

  # clean reads, keep only those starting with restriction site and clip to given length
  # NB, if using several restriction enzymes, you need to specify several sites
  # the correct syntax for EcoR1 (AATTC) and MSE (ATT) is: EcoR1Mse1
  # the correct syntax for EcoR1 is: EcoR1
  # the correct syntax for Sbf1 is Sbf1
  # my $motif = "Sbf1";
  my $command = "perl bin/FinalClean.pl $workfolder/tmp/merged $workfolder/stacks $motif $bclength $MaxLen $keepbc $trimR2 $stitch";
  print("\nRunning: $command\n\n");
  system($command);



  } else { #no barcodes, cleaning only, we bypass all that RADs-related mess, perform final renaming (fq_1 and fq_2) and move results to stacks folder
  my @allfiles = `find $workfolder/tmp/qualclean -type f -name '*.fastq'`;

  foreach my $path (@allfiles){
    chomp($path);
    $path =~ /\/\d*\/clean\.(.*)\.R(\d*)\.(.*)\.fastq/;
    my $tech = $1;
    my $readnr = $2;
    my $specnr = $3;
    
    my $command = "mv $path $workfolder/stacks/$tech\_$specnr.fq\_$readnr";
    print("\nRunning: $command\n\n");
    system($command);    
    }
  }


### STEP 3 Make quality stats (complete set of reads)
print("\nRunning: FASTQc $workfolder\n\n");
system("cat $workfolder/stacks/* | grep -v unmatched.R > $workfolder/logs/allfiles.fastq");
system("./bin/FastQC/fastqc $workfolder/logs/allfiles.fastq");
system("rm $workfolder/logs/allfiles.fastq");



### STEP 4 Make quality stats (specimen level)
my @list = `ls $workfolder/stacks`;
my %corresp;

foreach my $i (@list){
  chomp($i);
  $i =~ /(.*)\.fq\_(.*)/;
  $corresp{$1}{$2} = $i;
  }

foreach my $i (keys(%corresp)){
  my $file1 = $corresp{$i}{"1"};
  my $file2 = $corresp{$i}{"2"};
  system("cat $workfolder/stacks/$file1 $workfolder/stacks/$file2 > $workfolder/logs/complete.$i.fastq");
  system("./bin/FastQC/fastqc $workfolder/logs/complete.$i.fastq");
  system("rm $workfolder/logs/complete.$i.fastq");
  }

print "Done\n";



### STEP 5 Collect these stats 
my $command = "R CMD BATCH \"--args infolder='$workfolder/logs' outfolder='$workfolder' \" bin/QualStats.r";
print("\nRunning: $command\n\n");
system($command);

