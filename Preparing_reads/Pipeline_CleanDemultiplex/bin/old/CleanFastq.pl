#!/bin/perl

##### Part of Pipeline_CleanFastq.pl
### Assumes that barcodes are in param/barcodes.txt file
### Input fastq must me in data.in. Takes fastq files (not fast.gz)
### Outputs produced in data.out

### Uses eautils softs, must be already installed on system

##### Usage: perl CleanFastq.pl fileR1 fileR2 noFile barcode

$file1 = $ARGV[0];
$file2 = $ARGV[1];
$num = $ARGV[2];

## copy file to tmp and gunzip it
my $command = "cp data.in/$file1 data.in/$file2 tmp/";
print "Running:\n$command\n\n";
system($command);

chdir("tmp/");
my $command = "gunzip $file1 $file2";
print "Running:\n$command\n\n";
system($command);
chdir("..");

## account for filename changes
$file1 =~ s/\.gz//;
$file2 =~ s/\.gz//;

## demultiplex files, just this is needed actually
my $command = "fastq-multx -B param/barcodes.txt tmp/$file1 tmp/$file2 -m 2 -o data.out/%.R1.$num.fastq -o data.out/%.R2.$num.fastq";
print "Running:\n$command\n\n";
system($command);

## clean tmp files
my $command = "rm tmp/$file1 tmp/$file2";
print "Running:\n$command\n\n";
system($command);
