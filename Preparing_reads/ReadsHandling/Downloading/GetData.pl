#!/bin/perl

##### Downloads reads using a tab-delimited list of URL adresses (UHTS)
# 
# Input: $infile = tab-delimitted list of URLs to download
#
# Usage: perl CheckFileNames.pl infile
#
# Output: infile/ = folder, named as your infile, containing your reads
#
#
# Example of infile: just copy-paste the URLs from the UHTS
#
# http://uhts-lgtf.vital-it.ch/symlink/EphRd06_caXZT8KMb64V_L7_R1_001.fastq.gz  0.25 GB
# http://uhts-lgtf.vital-it.ch/symlink/EphRd06_YczaW7edQA1a_L7_R1_002.fastq.gz  0.25 GB
# http://uhts-lgtf.vital-it.ch/symlink/EphRd06_9s59jDLSr84o_L7_R1_003.fastq.gz  0.25 GB
# http://uhts-lgtf.vital-it.ch/symlink/EphRd06_OkZSNEBfdyz2_L7_R1_004.fastq.gz  0.25 GB
#
# NB only the first column of your infile is taken into account.
########################################################


my $infile = $ARGV[0];

open IN, "$infile";
$infile =~ s/\.txt$//;

system("mkdir $infile");
chdir($infile);
 
while(<IN>){
  chomp();
  my $tmp = $_;
  my @fields = split("\\s+", $tmp);
  my $url = $fields[0];
  chomp($url);
  system("wget $url");
  }
