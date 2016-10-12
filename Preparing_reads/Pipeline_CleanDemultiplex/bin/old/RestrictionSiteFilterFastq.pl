#!/usr/bin/perl -w

######################################################################################
# This script filters out read pairs that are not containing the restriction site.
# Ensures that ouputs are still synchronized.
# Performs a perfect match of the restiction site Motif (EcoRI: AATTC) 
#
# By default, the script assumes that the restriction site is on the read provided in the file Fastq1.
# 
# 1. Assumes fastq stored following the 4 lines format
#
#  @HWI-ST885:67:D0GY1ACXX:2:1101:1176:1964
#  NGATAACCAGCAAAGTCGGGGCAACAGTAGAGGAGGTCATTCCTCACGGCCAAGGCTGCTACGAGCTGAGAANGAGCACGGGAGAGGGGAAGATTGCGAT
#  +
#  #1=DDFFEHHHHHJIEGHIJJJIJJJJGGHIIHGIJBFHIGGHIIGIJJJIHFEFBEFDEDEDDDDDDCDCC#,8?B@DDDD9B<BBDD9<@ACDDDD>A
#
#
# 2. Assumes that read headers are identical among fastq files
#
#
# 3. Warning: the script keeps both fastq files in RAM, maybe be resource consuming
#
#
# Usage: perl Fastq1 Fastq2 Motif Outfolder
#
######################################################################################
use File::Basename;

#### Get script inputs
my $input1 = $ARGV[0];
my $input2 = $ARGV[1];
my $motif = $ARGV[2];
my $limlen = $ARGV[3];
my $outfolder = $ARGV[4];

my $bsn1 = basename($input1);
my $bsn2 = basename($input2);
chomp($motif);
chomp($limlen);




#### Get sequencer name (used for marking transitions from one accession to the next)
#### Load fastqs into hashes
## First pass: parse fastq1 and put in hash
open IN1, $input1;
my %store0;
my $cnt = 0;
my $nline = 0;
my $sequencer;

while(<IN1>){
  chomp();
  my $line = $_;

  if($nline == 0){ #Use the very first line of fastq to pick sequencer name
    $sequencer = $line;
    $sequencer =~ /(@[\w|\-,]+:\d+:\w+):.*/; # RegExp to catch sequencer name, here @VADER:36:C0MLBACXX @VADER:36:C0MLBACXX:8:2314:17186:5629 2:N:0:
    $sequencer = $1;
    }

  if($line =~ /^$sequencer/){ #we enter into a new accession
    $acc = $line; # get accession number
    $acc =~ /(^@.*)\s+.*/; #RegExp to clean accession name
    $acc = $1;
    $cnt = 0;
    }

  if($cnt < 4){ # store accession as long as there are < 4 lines (hence skip any weird 5th lines, if any)
    push(@{ $store0{$acc} }, $line);
    $cnt++;  
    }
  $nline++;
  }
close(IN1);


## Second pass: parse hash of Fastq1 and keep only reads starting with Restriction Site.
# Perform also sequence length check and trim sequence to desired length
my %store1;
foreach my $read (keys(%store0)){
  my $seq = $store0{$read} -> [1];
   
    # Remove barcode from sequence
    $seq = substr($seq, 6);
    $store0{$read} -> [1] = $seq;

  # Check that resctriction is there
 # if($seq =~ /^$motif/){ #We use a RegExp; hence we perform a perfect match here

    if(length($seq) == $limlen){ #Keep only sequences with correct length
      my $qual = $store0{$read} -> [3];
      $qual = substr($qual, 0, $limlen);
      $store0{$read} -> [3] = $qual;

      $store1{$read} = $store0{$read};      
      }
  #  }
  }
@accs1 = keys(%store1); #get list of accessions.
undef(%store0); # free some RAM



## Parse fastq2 and put in hash
open IN2, $input2;
$cnt = 0;
$nline = 0;

while(<IN2>){
  chomp();
  my $line = $_;

  if($nline == 0){ #Use the very first line of fastq to pick sequencer name
    $sequencer = $line;
    $sequencer =~ /(@[\w|\-,]+:\d+:\w+):.*/; # RegExp to catch sequencer name, here @VADER:36:C0MLBACXX @VADER:36:C0MLBACXX:8:2314:17186:5629 2:N:0:
    $sequencer = $1;
    }

  if($line =~ /^$sequencer/){ #we enter into a new accession
    $acc = $line; # get accession number
    $acc =~ /(^@.*)\s+.*/; #RegExp to clean accession name
    $acc = $1;
    $cnt = 0;
    }

  if($cnt < 4){ # store accession as long as there are < 4 lines (hence skip any weird 5th lines, if any)
    push(@{ $store0{$acc} }, $line);
    $cnt++;  
    }
  $nline++;
  }
close(IN2);

## Second pass: perform sequence length check and trim sequence to desired length
my %store2;
foreach my $read (keys(%store0)){
  my $seq = $store0{$read} -> [1];

  # Trim sequence to needed length
  $seq = substr($seq, 6, $limlen);
  $store0{$read} -> [1] = $seq;

  if(length($seq) == $limlen){ #Keep only sequences with correct length
    my $qual = $store0{$read} -> [3];
    $qual = substr($qual, 6, $limlen);
    $store0{$read} -> [3] = $qual;

    $store2{$read} = $store0{$read};      
    }
  }
@accs2 = keys(%store2); #get list of accessions.
undef(%store0); # free some RAM



#### Find intersecting / orphan reads
my @pairs = my @orphs = ();
my %count = ();

foreach $e (@accs1, @accs2) { $count{$e}++ }
foreach $e (keys %count) {
    if ($count{$e} == 2) {
        push @pairs, $e;
    } else {
        push @orphs, $e;
    }
}

my $npairs = $#pairs;
my $norphs = $#orphs;
print "#### Checked $input1 vs $input2\nFound $npairs paired reads\n      $norphs orphan reads\n";


#### Produce outputs
if($npairs > 0){
  open(OUT1, ">$outfolder/$bsn1");
  open(OUT2, ">$outfolder/$bsn2");   
  
  foreach $acc (@pairs){
    $fq1 = $store1{$acc};
    print OUT1 join("\n", @{$fq1}), "\n";
    $fq2 = $store2{$acc};
    print OUT2 join("\n", @{$fq2}), "\n";
    }
  close(OUT1);
  close(OUT2);
  }

# if($norphs > 0){
#   open(OUT, ">orphs.$input1");
# 
#   foreach $acc (@orphs){
#     if($store1{$acc}){
#       $fq1 = $store1{$acc};
#       print OUT join("\n", @{$fq1}), "\n";
#       }
#     if($store2{$acc}){
#       $fq2 = $store2{$acc};
#       print OUT join("\n", @{$fq2}), "\n";
#       }  
#     }
#   close(OUT);
#   }

