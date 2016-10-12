#!/bin/perl

##### Script pilot of Pipeline_CleanFastq.pl
# Assumes that barcodes / primers are in param/barcodes.txt and param/primers.fas files
# 
# Input fastq must me in $inputfolder. 
# 
# TAKES fastq.gz files
# Outputs produced in data.out/outfolder
#
# Launches cleaning process on as many as possible CPUs. Upper CPU usage limit must be given
#
# Usage: perl CleanFastq_build.pl inputfolder barcodes outfolder maxuse
########################################################
use threads;
use threads::shared;

my $maxuse = $ARGV[0];

### Loop over each file pair
my @threads;

foreach my $i (1..10){
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
    my $t = threads->new(\&rundataset);
    push(@threads, $t);
    sleep(1);

    } else {
    print "WAIT... CPU load is maximised\n\n\n";
    sleep(1);
    redo;
    }
  }



### Survey ongoing threads
foreach (@threads) {
  my $genome = $_-> join;
  print "done with $genome\n";
  }



### Wait at that point untill all threads are over
# sleep 1 until threads->exited;



### STEP 2 Collect all fastq and merge them at the sample level
# perform also quality checks and get stats.
print("\nAll threads are done...\n\n");




### ROUTINES
sub rundataset {
  print("Start thread\n");
  sleep(20)
  }
