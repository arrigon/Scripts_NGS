#!/usr/bin/perl -w

######################################################################################
# This script extracts reads from a fastq file, using read IDs specified from an external list
# The fastq files need to be in sync, no order check is performed, the file is visited only once and nothing is kept in RAM
#
### Arguments
# -fastq1
# -IDlist
# -outfolder
#
### Specifications
# 1. Assumes fastq stored following the 4 lines format
#
#  @HWI-ST885:67:D0GY1ACXX:2:1101:1176:1964
#  NGATAACCAGCAAAGTCGGGGCAACAGTAGAGGAGGTCATTCCTCACGGCCAAGGCTGCTACGAGCTGAGAANGAGCACGGGAGAGGGGAAGATTGCGAT
#  +
#  #1=DDFFEHHHHHJIEGHIJJJIJJJJGGHIIHGIJBFHIGGHIIGIJJJIHFEFBEFDEDEDDDDDDCDCC#,8?B@DDDD9B<BBDD9<@ACDDDD>A
#
# 2. Assumes that read headers are identical among fastq files
#
######################################################################################
use File::Basename;

#### Get script inputs
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $IDlist = $ARGV[2];


#### load list of IDs
open IN, $IDlist;
my %listID;
while(<IN>){
  chomp();
  my @fields = split("\t", $_);
  my $clstr = $fields[0];
  $clstr =~ s/^\>//;
  $clstr =~ s/ //;

  my $ID = $fields[1];
  $ID =~ s/^\>//;

  $listID{$ID} = $clstr;
  }
my $listref = \%listID;


#### load desired reads in RAM
my $ref = IOfastq($infile, "efqr", "$infile.fq", $listref);
my %fastq = %{$ref};


#### Parse these guys and save to correct outfile
foreach $acc (keys(%fastq)){
  # retrieve data
  my $seq = $fasta{$acc}{"seq"}
  my $qual = $fasta{$acc}{"qual"}

  # print to outfile
  $clstr = $listID{$acc};
  open(OUT, ">>$output.$clstr.fas");
  print OUT ">$acc\n$seq\n+\n$qual\n";
  close(OUT);
  }


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
  # $ref = IOfastq($infile, "efar", "outfile", $listref);
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
  my %list = %{$listref};
  my %fastq;
  my %fasta;
  my @index;
  my $acc;
  my $seq;
  my $qual;
  my $cnt = 0;

  
  ## prepare outfile, if needed
  if($job =~ /p$/){
    open OUT, ">$outfile";
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
	    print OUT "\@$acc\n$seq\n+\n$qual\n";
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
  close(IN);  
  close(OUT);
  
  ## Prepare outputs 
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