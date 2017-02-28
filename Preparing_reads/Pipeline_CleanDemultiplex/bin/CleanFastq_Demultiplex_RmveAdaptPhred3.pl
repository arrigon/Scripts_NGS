#!/bin/perl

####################################
# Demultiplexing ONE pair of FASTQ files (F W reading, must be in sync)
#
### Params
# - Takes fastq files (not fast.gz)
# - Assumes that barcodes are in param/barcodes.txt file, by default, keeps barcodes
# - Assumes that primers  are in param/Primers.fas
#
# - Uses eautils softs (fastq-multx), provided in bin folder.
#
### Outputs: all files stored in following folder system
#		- demultiplexed files: $outfolder/multiplex/$num
#		- final files: $outfolder/qualclean/$num
#
#
### Usage: perl CleanFastq_Demultiplex_RmveAdaptPhred.pl fileR1 fileR2 noFile outfolder keepbc
#
# Nils Arrigo. Uni Lausanne 2012
############################################################
use File::Basename;

### Get script arguments
my $file1 = $ARGV[0];
my $bsnfile1 = basename($file1);

my $file2 = $ARGV[1];
my $bsnfile2 = basename($file2);

my $num = $ARGV[2]; #file number

my $outfolder = $ARGV[3];
chomp($outfolder);

my $outlogsfolder = $ARGV[4];

my $barcodes = $ARGV[5];
chomp($barcodes);

my $mismatch = $ARGV[6];
chomp($mismatch);

my $minlen = $ARGV[7];
chomp($minlen);

my $phred = $ARGV[8];
chomp($phred);


# prepare outfolders
my $outfastq = "$outfolder/fastq/$num";
my $outmultiplex = "$outfolder/multiplex/$num";
my $outrepaired = "$outfolder/repaired/$num";
my $outqualclean = "$outfolder/qualclean/$num";


### Prepare subfolders
system("mkdir -p $outfolder $outlogsfolder $outlogsfolder/qualclean $outlogsfolder/multiplex $outfastq $outmultiplex $outrepaired $outqualclean");




### STEP 0 - Decompress files
my $command = "zcat -c $file1.gz > $outfastq/$bsnfile1";
print("\nRunning: $command\n\n");
system($command);

my $command = "zcat -c $file2.gz > $outfastq/$bsnfile2";
print("\nRunning: $command\n\n");
system($command); 




### STEP 1 - demultiplex files
if($barcodes eq 0){
  # simply move files to next folder, make sure they comply with the CheckFileName format
  my $command = "mv $outfastq/$bsnfile1 $outfastq/$bsnfile2 $outmultiplex/";
  print "Running: $command\n\n";
  system($command);
  
  } else {
  # Provide file1 and file2 simultaneously, to keep files in sync
  # We keep the barcode in the sequence, to make the debugging easier if needed
  my $command = "./bin/fastq-multx -b -x -B $barcodes $outfastq/$bsnfile1 $outfastq/$bsnfile2 -m $mismatch -o $outmultiplex/%.R1.$num.fastq -o $outmultiplex/%.R2.$num.fastq > $outlogsfolder/multiplex/demultiplex.$num.log";
  print "Running: $command\n\n";
  system($command);

  }




### STEP 2 - clip adapters and clean based on PHRED scores
# Get files from $outmultiplex, there is as much files as barcodes
# we first must list them and get pairs R1 vs R2
# then provide cleaning program with pairs of files, to keep everything in sync

# Get files in $outmultiplex
my @list = `ls $outmultiplex`;

# reconstruct pairs of files
my %corresp;
my %chuncks;
foreach my $i (@list){
  chomp($i);
  if($i =~ /(.*)\.(.*)\.(.*)\.fastq/){
    $corresp{$1}{$3}{$2} = $i;
    }
  }

# loop over them, pair by pair
foreach my $i (keys(%corresp)){
  my %tmp = %{$corresp{$i}};
  
  foreach my $j (keys(%tmp)){
    
    # get file names
    my $file1 = $corresp{$i}{$j}{"R1"};
    my $file2 = $corresp{$i}{$j}{"R2"};

    # repair, clean and sync fastq files (demultiplex introduces weird things)
    my$command = "perl bin/RepairSyncFastq2.pl $outmultiplex/$file1 $outmultiplex/$file2 $outrepaired/";
    print "Running: $command\n\n";
    system($command);
    
    # clean adapters and PHREDs
    # Take care here: make sure that the skewing trimming is not enabled
    # this would cut the start of the read since it is highly repeated and remove the barcode + restriction site.
    # While this is not a problem in genome assembly, this destroys our stacks
    my$command = "./bin/fastq-mcf -f -p 10 -F param/Contaminants.fas -q $phred -k 0 -l $minlen -o $outqualclean/clean.$file1 -o $outqualclean/clean.$file2 param/Primers.fas $outrepaired/$file1 $outrepaired/$file2 > $outlogsfolder/qualclean/Clean.$file1.log";
    print "Running: $command\n\n";
    system($command);
    }
  }



### STEP 4 clean mess
system("rm -rf $outfastq");
system("rm -rf $outrepaired $outmultiplex");
system("rm $outqualclean/clean.unmatched.*");