#!/usr/bin/perl -w

######################################################################################
# This script synchronizes read pairs from distinct fastq files
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
# Usage: 
# mkdir outfolder
# perl RepairSyncFastq2.pl fastq.fq1 fastq.fq2 outfolder
#
######################################################################################
use File::Basename;

#### Get script inputs
my $input1 = $ARGV[0];
my $input2 = $ARGV[1];
my $outfolder = $ARGV[2];

my $bsn1 = basename($input1);
my $bsn2 = basename($input2);


#### List read names from fastq
$ref1 = IOfastq($input1, "acr");
$ref2 = IOfastq($input2, "acr");

my @accs1 = @{$ref1};
my @accs2 = @{$ref2};


#### Find intersecting / orphan reads
my %pairs;
my %orphs;
my %count;
my $npairs = 0;
my $norphs = 0;

foreach $e (@accs1, @accs2) { 
  chomp($e);
  $count{$e}++;
  }
  
foreach $e (keys %count) {
    if ($count{$e} == 2) {
        $pairs{$e} = 1;
        $npairs++;
    } else {
        $orphs{$e} = 1;
	$norphs++;
    }
}
print "#### Checked $input1 vs $input2\nFound $npairs paired reads\n      $norphs orphan reads\n";


#### Produce outputs
## save paired reads
my $pairsref = \%pairs;

my $outfile1 = "$outfolder/$bsn1";
my $outfile2 = "$outfolder/$bsn2";

IOfastq($input1, "efqp", $outfile1, $pairsref);
IOfastq($input2, "efqp", $outfile2, $pairsref);


## save orphaned reads (uncomment to use)
# my $orphsref = \%orphs;
# $outfile1 = "$outfolder/orph.$bsn1";
# $outfile2 = "$outfolder/orph.$bsn2";
# IOfastq($input1, "efqp", $outfile1, $orphsref);
# IOfastq($input2, "efqp", $outfile2, $orphsref);



### Routine
sub IOfastq {
  ################### Usage #######################
  #################################################
  # IOfastq(infile, job, outfile, index)
  #
  # infile = path to fastq infile
  #
  # job = [e][fa|fq|ac][p|r]
  # 	  [e] = optional: extract only subset of reads (to be provided in index)
  # 	  [fa|fq|ac] = returned info: fasta, fastq, accession
  # 	  [p|r] = output; either p to outfile, or r to ram, returns a reference to hash or array
  #
  # outfile = optional, path to outfile
  #
  # index = optional, list of reads to extract from fastq. MUST be provided as reference to indexing hash, 
  #        where the desired read names are loaded as hash's keys. e.g. $idx{"FCD0CBVABXX:1:1101:6477:2464#ATCACGAT/2"} = 1;
  # 
  # Nils Arrigo, Unil 2014
  ################## Examples #####################
  #################################################
  ## prepare indexing hash
  # my %idx;
  # $idx{"FCD0CBVABXX:1:1101:6477:2464#ATCACGAT/2"} = 1;
  # $idx{"FCD0CBVABXX:1:1101:6489:2288#ATCACGAT/2"} = 1;
  # $idx{"FCD0CBVABXX:1:1101:6712:2313#ATCACGAT/2"} = 1;
  # $listref = \%idx;
  #
  ## load into RAM, complete fastq
  # $ref = IOfastq($infile, "far");
  # $ref = IOfastq($infile, "fqr");
  # $ref = IOfastq($infile, "acr");
  #
  ## load into RAM, subset of reads
  # $ref = IOfastq($infile, "efar", "$infile.fas", $listref);
  # $ref = IOfastq($infile, "efqr", "$infile.fq", $listref);
  # $ref = IOfastq($infile, "eacr", "$infile.list", $listref);
  # 
  # # pipe to file, complete fastq
  # IOfastq($infile, "fap", $outfile); #fasta
  # IOfastq($infile, "fqp", $outfile); #fastq
  # IOfastq($infile, "acp", $outfile); #reads list
  # 
  # # pipe to file, subset of reads
  # IOfastq($infile, "efap", $outfile, $listref);
  # IOfastq($infile, "efqp", $outfile, $listref);
  # IOfastq($infile, "eacp", $outfile, $listref);
  #################################################
  
  my ($infile, $job, $outfile, $listref) = ($_[0], $_[1], $_[2], $_[3]);
  chomp($job);
  
  ## open infile
  open(IN, "$infile");
  
  ## init values
  my %list;
  my %fastq;
  my %fasta;
  my @index;
  my $acc;
  my $seq;
  my $qual;
  my $cnt = 0;
  my @accorder;
  
  ## prepare outfile, if needed
  if($job =~ /p$/){
    open OUT, ">$outfile";
    }
  

  ## prepare listref, if needed
  if($job =~ /^e/){
    %list = %{$listref};
    }
    
    
  ## parse fastq, line by line. Check consistency during parsing
  while(<IN>){
    # get focal line
    chomp();
    $line = $_;
    
    
    # If we look at accession
    if($cnt == 0){ 
      if($line =~ /^@/){ #make sure it starts with @ character
	$acc = $line;
	$acc =~ s/^@//;
	$acc =~ /([\w|\d|\-\:]*)\s*.*/;
	$acc = $1;
	
	push(@accorder, $acc);
	
	
	if($job =~ /^e/){ #if job has to do with Filtering reads
	  if(!defined($list{$acc})){
	    next;
	    }
	  }
	} else {
	next;
	}
      } #now will go to $cnt++, just before end of while

      
    # If we look at sequence
    if($cnt == 1){ 
      $seq = $line;
      } #now will go to $cnt++, just before end of while
      
      
    # If we look at separator
    if($cnt == 2){ 
      if($line !~ /^\+$/){ #must be + character
	$cnt = 0;
	next;
	}
      } #now will go to $cnt++, just before end of while
      

    # If we look at qual
    if($cnt == 3){
      $qual = $line;
      if(length($qual) != length($seq)){ #seq and qual must have the same length
	$cnt = 0;
	next;
	} else {
	
	# prepare outputs according to desired job
	if($job =~ /.*fq.*/){
	  if($job =~ /r$/){
	    $fastq{$acc}{"seq"} = $seq;
	    $fastq{$acc}{"qual"} = $qual;
	    }
	  if($job =~ /p$/){
	    if(defined($seq)){
	      print OUT "\@$acc\n$seq\n+\n$qual\n";
	      }
	    }
	  } 
	  
	if($job =~ /.*fa.*/){
	  if($job =~ /r$/){
	    $fasta{$acc} = $seq;
	    }
 	  if($job =~ /p$/){
	    print OUT ">$acc\n$seq\n";
	    }
	  } 
	
	if($job =~ /.*ac.*/){
	  if($job =~ /r$/){
	    push(@index, $acc);
	    }
	  if($job =~ /p$/){
	    print OUT "$acc\n";
	    }
	  }
	  
	# move to next accession
	$cnt = 0;
	next;
	}
      } #skips $cnt++, because of next
    
    $cnt++;
    }
    
  close(OUT);
  
  ## Prepare outputs 
  $fastq{"accorder"} = \@accorder;
  
  if($job =~ /fqr$/){
    return(\%fastq);
    } 
    
  if($job =~ /far$/){
    return(\%fasta);
    } 
  
  if($job =~ /acr$/){
    return(\@index);
    }
    
  if($job =~ /p$/){
    return(1);
    }
  }

  
  