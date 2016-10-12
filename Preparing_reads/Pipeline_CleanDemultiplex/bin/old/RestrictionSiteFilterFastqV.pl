#!/usr/bin/perl -w

######################################################################################
# This script filters out read pairs that are not containing the restriction site.
# Ensures that ouputs are still synchronized.
# 1. Clips the beginning of reads to get rid (or keep) of barcode sequence (parameter BarcodeLength + keepbc)
# 2. Performs a PERFECT match of the restiction site Motif (parameter Motif) 
#    NB if using several restriction enzymes, you need to specify several sites
#    the correct syntax for EcoR1 (AATTC) + MSE (TTAA) for motif parameter is "^AATTC|^TTAA"
# 3. Trims all reads at the same length (MaxLength, not including Barcode)
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
# 2. Assumes that read headers are identical among fastq files
#
# 3. Warning: the script keeps both fastq files in RAM, maybe be resource consuming
#
# Usage: perl RestrictionSiteFilterFastq.pl Fastq1 Fastq2 Motif BarcodeLength MaxLength Outfolder KeepBarcode
#
# Examples:
# Single digest with EcoR1 (AATTC)
# perl RestrictionSiteFilterFastq.pl FastqWithBarcodesAndRestrictionSite.fastq Fastq2.fastq Enzyme 5 85 Outfolder 0 0 0
#
# Double digest EcoR1 (AATTC) + Mse1 (TTAA), paired-end, trim R1 and R2 to limlen, and stitch R1R2 altogether in final file
# perl RestrictionSiteFilterFastq.pl FastqWithBarcodesAndRestrictionSite.fastq Fastq2.fastq EcoR1 5 85 Outfolder 0 1 1
#######################################################################################

use File::Basename;

#### Get script inputs
my $input1 = $ARGV[0]; # Read R1, contains barcode and restriction site
my $input2 = $ARGV[1]; # Read R2
my $motif = $ARGV[2];  # Enzyme: EcoR1, Sbf1, add whatever is needed, "keepall" = do not apply this filter
my $bcdelgth = $ARGV[3]; # Length of barcode
my $limlen = $ARGV[4]; # Maximal read length (barcode not counted), use limlen < 0 to bypass this cut.
my $outfolder = $ARGV[5]; # outfolder
my $keepbc = $ARGV[6]; # Keep (1) or remove (0) barcode
# my $trimR2 = $ARGV[7]; # Trim ($trimR2 = 1) second read to $limlen (used in double-digest + paired-end RADs) or leave as is ($trimR2 = 0)
# my $stitch = $ARGV[8]; # Append ($stitch = 1) R1 and R2 altogether (used in double-digest + paired-end RADs) or store R1 and R2 separately ($stitch = 0)


my $bsn1 = basename($input1);
my $bsn2 = basename($input2);
chomp($motif);
chomp($limlen);


#### Restriction sites to look for, add yours in here
if($motif =~ /EcoR1/){
  $motif = "AATTC";
  }

if($motif =~ /EcoR1Mse1/){
  $motif = "AATTC";
  }

if($motif =~ /Sbf1/){
  $motif = "TGCAGG";
  }


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
    $sequencer =~ /(@[\w|\-,]+:\d+:[\w|\-]+):.*/; # RegExp to catch sequencer name, here @VADER:36:C0MLBACXX @VADER:36:C0MLBACXX:8:2314:17186:5629 2:N:0:
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
my $limlenseq;
my $startpos;
foreach my $read (keys(%store0)){
  my $seq = $store0{$read} -> [1]; 
  my $qual = $store0{$read} -> [3];   

  my $seqshort;
  my $qualshort;

  my $startseq = substr($seq, 0, 25); #focus on start of sequence; used for looking for restriction sites

  ### Treatment 1: keep barcode or not
  if($keepbc == 0){ #discard barcode
    $seqshort = substr($seq, $bcdelgth); # using substr(seq, startpos) -> no stop pos defined, so goes up to end of string
    $qualshort = substr($qual, $bcdelgth); #do the same to qual infos 
    } else { #keep it
    $seqshort = $seq; # just pass on original sequence
    $qualshort = $qual; # do the same to qual infos 
    }

  ### Treatment 2: look for restriction site (if $motif is not "keepall")
  if($motif ne "keepall"){
    if($startseq =~ m/$motif/g){ # look for restriction motif
	$startpos = pos($startseq) - length($motif); #record where the restriction site starts
	$seqshort = substr($seq, $startpos); # trim to startpos (barcode will be discarded by doing so)
	$qualshort = substr($qual, $startpos); # do the same to qual infos 
	} else { # no restriction site found, stick to original sequence
	$seqshort = $seq; # just pass on original sequence
	$qualshort = $qual; # do the same to qual infos 
	}
      }

  ### Treatment 3: trim reads at uniform length
  if($limlen > 0){ # applies only if $limlen is defined
      $seqshort = substr($seqshort, 0, $limlen); # trim to limlen (does not account for barcode if it is still in the sequence)
      $qualshort = substr($qualshort, 0, $limlen); #do the same to qual infos 
      }


  ### Now let's inspect $seqshort
  # at that stage, $seqshort and $qualshort have been treated according to needs, we need to filter out non-desired reads
  
  # Filter 1: make sure that it starts with a restriction site (bypass filter if $motif = "keepall")
  if($seqshort =~ m/$motif/g or $motif eq "keepall"){ #This switch decides whether the focal sequence is passed on to %store1 (is ignored otherwise)
						     # if $motif = say AATTC, the regexp will look for restriction sites at beginning of sequence. If complies, the sequences is passed on. if not, it left aside.
						     # otherwise, if $motif = keepall, all sequences will be kept and passed into %store1
						     # this switch makes sure that either $motif =~ restr.site or $motif = "keepall" passes the focal sequence to the next stage
    
    # Filter 2: Check that read has desired length (barcode length is NOT included in that check)
    if($limlen > 0){ # limlen filtering active
      if(length($seqshort) == $limlen){ 
	$store1{$read} = $store0{$read}; #Yes, we pass on that sequence into %store1 and update its data accordingly
	$store1{$read} -> [1] = $seqshort; #just pass on $seqshort, it is already fine (i.e. BC removed and trimmed accordingly to limlen)
	$store1{$read} -> [3] = $qualshort; #idem with qual infos
	}
      } else { # limlen filtering inactive
	$store1{$read} = $store0{$read}; #Yes, we pass on that sequence into %store1 and update its data accordingly
	$store1{$read} -> [1] = $seqshort; #just pass on $seqshort, it is already fine (i.e. BC removed and trimmed accordingly to limlen)
	$store1{$read} -> [3] = $qualshort; #idem with qual infos
      } #end of Filter 2
    } #end of Filter 1
  } #end of looping over reads
@accs1 = keys(%store1); #get list of accessions.
undef(%store0); # free some RAM


#####
## Parse fastq2 and put in hash
open IN2, $input2;
$cnt = 0;
$nline = 0;

while(<IN2>){
  chomp();
  my $line = $_;

  if($nline == 0){ #Use the very first line of fastq to pick sequencer name
    $sequencer = $line;
    $sequencer =~ /(@[\w|\-,]+:\d+:[\w|\-]+):.*/; # RegExp to catch sequencer name, here @VADER:36:C0MLBACXX @VADER:36:C0MLBACXX:8:2314:17186:5629 2:N:0:
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


if($trimR2 == 1){
  # Second pass: perform sequence length check and trim sequence to desired length 
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

  } else {
  my %store2 = %store0;
  @accs2 = keys(%store2); #get list of accessions.
  undef(%store0); # free some RAM
  }


#### Resync both fastq files (find intersecting / orphan reads)
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

my $npairs = $#pairs + 1;
my $norphs = $#orphs + 1;
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

