#!/bin/perl

####################################
# Demultiplexing pair of FASTQ files
#
### Params
# - Takes fastq files (not fast.gz)
# - Assumes that barcodes are in param/barcodes.txt file
# - Outputs produced in outfolder
# - Uses eautils softs (fastq-multx), provided in bin folder.
#
### Usage: perl CleanFastq.pl fileR1 fileR2 noFile outfolder
#
# Nils Arrigo. Uni Lausanne 2012
############################################################


### Get script arguments
my $file1 = $ARGV[0];
my $file2 = $ARGV[1];
my $num = $ARGV[2]; #file number
my $outfolder = $ARGV[3];

# clean them a bit
my $file1 =~ s/\.gz//;
my $file2 =~ s/\.gz//;
chomp($outfolder);

# demultiplex files
my $command = "./bin/fastq-multx -B param/barcodes.txt $file1 $file2 -m 2 -o $outfolder/%.R1.$num.fastq.gz -o $outfolder/%.R2.$num.fastq.gz";
print "Running:\n$command\n\n";
system($command);