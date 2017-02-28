#!/bin/perl

##### Script pilot of Pipeline_CleanFastq.pl
### Performed steps:
# 1. Demultiplex (fastq-multx)
# 2. Clean reads according to PHRED scores (fastq-clean), clip out primers (see below)
# 3. Clip out barcodes (default, but see below), discard reads not starting with restriction site (but see below) 
#    and trim all R1 reads to the same length (R2 reads kept at varying lengths, but see below)
# 4. Collect and produce quality check stats 
#
### Parameters
# Usage: perl CleanFastq_build.pl inputfolder barcodes outfolder maxuse
# - inputfolder = repository of fastq files (keep them gzipped). 
#              WARNING: specimen names should NOT contain dots.
# - barcodes = path to file containing barcodes, by default should be stored in params/ folder, 
#              put 0 to skip demultiplexing -> will move directly to cleaning steps
#              WARNING: will not merge files if specimens were stored across several fastqs
#              WARNING: quite sensitive about file naming. 
#			If deciding to skip the demultiplexing, YOU MUST MAKE SURE THAT INPUT FILES 
# 			HAVE A NAME FORMAT THAT COMPLIES WITH CheckFileNames.pl script (i.e. TECHSTUFF.R1|2.SPECNR.fastq.gz")
# - outfolder = name of output folder. Results will be saved in data.out/outfolder
# - maxuse = Maximum number of CPUs allowed to run in parallel on the server
#
### General settings
# Barcodes: should be in param/xxx_barcodes.txt (two columns file, col 1 = specimen name, col2 = barcode, no headers)
#	    by default, the barcodes are clipped out from the cleaned reads (keepbc = 0)
#           this can be changed at line 41, by setting $keepbc = 1
#
# Primers : should be in param/Primers.fas
# Inputs: Works for paired-end data. For single end datasets, you need to fake a set of R2 reads (copy + rename files of R1 reads)
#         NB fastq.gz files
# Working scripts / external programs: must be in bin/ folder
#
### Dependencies:
# perl, R CRAN
# perl libraries: threads, File::Basename, 
#
### Usage
# perl CleanFastq_build.pl /data/Ber1 param/Ber1_barcodes.txt Ber1 10
# WARNING: several addidional params are left hard coded at the beginning of that script (see below)
# 
########################################################
use threads;
use threads::shared;



### WARNING Hard coded parameters WARNING: Juniper with Sbf1
my $bcmismatch = 2; #Maximum number of allowed mismatches in barcode
my $keepbc = 0; #keep ($keepbc = 1) or discard ($keepbc = 0) barcode
my $motif = "Sbf1"; # For stacks only: keep only reads starting with restriction enzyme. Sbf1, EcoR1 are available so far. If needing additional motifs, add them to bin/RestrictionSiteFilterFastqIII.pl
		    # N.B. only using $motif = keepall works in combination with barcode clipping options 
		    # i.e. $motif = "keepall" and $keepbc = 1 will leave the barcode in the reads
		    #      $motif = "keepall" and $keepbc = 0 will ALWAYS trim out the barcode
		    #      $motif = EcoR1 or others will always trim out the barcode
my $MaxLen = 130; # For stacks only: trim all reads to that length (usually 80), discard shorter reads. Put MaxLen = -1 to bypass that treatment
my $MinLen = 130; #discard reads shorter than that value (fastq-mcf step)
my $phred = 28;  #trim end of reads if phred scores goes below that value
my $namestyle = "checked"; # make sure we use the right RegExp to index raw fastq files (TODO). 
			  # Best option so far: rename your raw fastq files BEFORE starting, in order to comply with "checked" format
			  # To do so: use script bin/CheckFileNames.pl and read help in there to have it at work.
			  # format = checked (as produced by bin/CheckFileNames.pl): fastq files are as XXX.R1.001.fastq.gz, where XXX is laneID, R1 = readID (R1/R2), 001 = specimenID / fastqchunckID
			  # format = CIG : fastq files are named as XXX_YYY_R1_001.fastq.gz, where XXX is laneID, YYY is skipped, R1 = readID (R1/R2), 001 = specimen / chunckID
my $trimR2 = 0;	# Trim ($trimR2 = 1) second read to $MaxLen (used in double-digest + paired-end RADs) or leave as is ($trimR2 = 0)
		# Note that is inactive if MaxLen < -1.
my $stitch = 0; # Append ($stitch = 1) R1 (orientation standard) and R2 (revcomp) altogether in a single fastq file1  (used in double-digest + paired-end RADs)
		# if $stitch = 0, store R1 (orient. std) and R2 (orient. std) separately


### Get script arguments
my $inputfolder = $ARGV[0];
my $barcodes = $ARGV[1];
my $outfolder = $ARGV[2];
my $maxuse = $ARGV[3];



### Print logfile
system("mkdir data.out/$outfolder data.out/$outfolder/stacks");
open(LOG, ">data.out/$outfolder/CommandLog.txt");
print LOG "perl CleanFastq_Pipeline_queue3.pl $inputfolder $barcodes $outfolder $maxuse\nUsing:\nMotif = $motif\nMinLen = $MinLen\nMaxLen = $MaxLen\nBarcodeMismatches = $bcmismatch\nMinPhredQual = $phred\nkeepbc = $keepbc\ntrimR2 = $trimR2\nstichR1R2 = $stitch\n";
close(LOG);

### Get list of fastq files in data.in
### Reconstruct pairs of files (based on RegExp, specific to naming scheme of files)
my @list = `ls $inputfolder`;

my %corresp; # WARNING: this is where the namestyle trick happens. ADAPT REGEXP according to needs (so far, set up to comply with "checked" namesyles)
foreach my $i (@list){
  chomp($i);
  if($i =~ /(.*)\.(R.*)\.(.*)\.fastq\.gz/){ #format CheckFileNames.pl $1 = laneID, $2 = readID (R1/R2), $3 = specimen / chunckID
#   if($i =~ /(.*)_.*_(.*)_(.*)\.fastq.gz/){ #format CIG
#   if($i =~ /(.*_20x)_.*_(.*)_(.*)\.fastq.gz/){ #format CIG weirdo
    $corresp{$1}{$3}{$2} = $i;
    }
  }



### Loop over each file pair
my @threads;

foreach my $lane (keys(%corresp)){
  foreach my $i (keys(%{$corresp{$lane}})){
    my $file1 = $corresp{$lane}{$i}{"R1"};
    my $file2 = $corresp{$lane}{$i}{"R2"};

    # rename files
    $file1 =~ s/\.gz//;
    $file2 =~ s/\.gz//;


    # Get CPU load
    open PIPE, "uptime |";
    my $line = <PIPE>;
    close PIPE;
    $line =~ s/\s//g;
    my @lineArr =  split /:/, $line;
    my $times = $lineArr[@lineArr-1];
    my @timeArr = split /,/, $times;
    my $load = $timeArr[0] + 1;
    print "Current CPU load is $load \n\n\n";


    if($load < $maxuse) {
      print "OK --- start $file1 and $file2\n";

      ### STEP 1 - demultiplex, clean adapters and clip according to PHREDs
      my $file1 = "$inputfolder/$file1";
      my $file2 = "$inputfolder/$file2";
      my $num = $i;

      my $t = threads->new(\&rundataset, $file1, $file2, $num);
      push(@threads, $t);
      sleep(30);

      } else {
      print "WAIT... CPU load is maximised\n\n\n";
      sleep(60);
      redo;
      }
    }
  }


### Survey ongoing threads
foreach (@threads) {
  my $genome = $_-> join;
  print "done with $genome\n";
  }



### STEP 2 Collect all fastq and merge them at the sample level
# perform also quality checks and get stats.
# Trim all reads to 80bp (not counting barcode, but including restriction site), discard shorter reads.
my $command = "perl bin/FinalWrapII.pl data.out/$outfolder $barcodes $MaxLen $keepbc $motif $trimR2 $stitch";
print("\nRunning: $command\n\n");
system($command);



### ROUTINES
sub rundataset {
  my ($file1, $file2, $num) = ($_[0], $_[1], $_[2]);

  my $command = "perl bin/CleanFastq_Demultiplex_RmveAdaptPhred3.pl $file1 $file2 $num data.out/$outfolder/tmp data.out/$outfolder/logs $barcodes $bcmismatch $MinLen $phred";
  print("\nRunning: $command\n\n");
  system($command);
  }

