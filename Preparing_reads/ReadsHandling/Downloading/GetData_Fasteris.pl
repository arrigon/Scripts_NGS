#!/bin/perl

my $infile = $ARGV[0];

open IN, "$infile";
$infile =~ s/\.txt$//;

system("mkdir $infile");
chdir($infile);
 
while(<IN>){
  chomp();
  my $url = $_;
  my @fields = split("\w+", $url);
  my $url = $fields[0];
  chomp($url);
  my $command = "wget --http-user=nadir.alvarez\@gmail.com --http-password=SrbV1K3L --no-check-certificate $url";
  print "$command\n";
  system($command);
    }

# "wget --http-user=nadir.alvarez\@gmail.com --http-password=SrbV1K3L --no-check-certificate https://data.fasteris.com/private/HBQ/HBQ-1-24/data/131029_SN365_A_L002_HBQ-3_R2.fastq.gz"