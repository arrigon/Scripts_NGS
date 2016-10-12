#!/bin/perl

##### Part of Pipeline_CleanFastq.pl
### Assumes that barcodes are in param/barcodes.txt file
### Input fastq must me in data.in. Takes fastq files (not fast.gz)
### Outputs produced in data.out

### Uses eautils softs, must be already installed on system
##### Usage: perl CleanFastq.pl fileR1 fileR2 noFile barcode

$file1 = $ARGV[0];
$file2 = $ARGV[1];

$file1 =~ s/\.gz//;
$file2 =~ s/\.gz//;

## demultiplex files, just this is needed actually
my$command = "fastq-mcf -l 75 -o data.phred/clean.$file1.fastq -o data.phred/clean.$file2.fastq param/Primers.fas data.mult/$file1 data.mult/$file2"
print "Running:\n$command\n\n";
system($command);
