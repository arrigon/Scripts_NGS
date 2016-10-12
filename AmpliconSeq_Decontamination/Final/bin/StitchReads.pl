#!/usr/bin/perl -w

use File::Basename;

#### Get script inputs
my $input1 = $ARGV[0]; # Read R1
my $input2 = $ARGV[1]; # Read R2
my $trim = $ARGV[2]; # Trimming length of reads (applied to F and R)
my $outfile = $ARGV[3]; # outfile


#### Script params
my $bsn1 = basename($input1);
my $bsn2 = basename($input2);


#### load reads (trim at given length)
%store1 = %{IOfastq("$input1", "fqr", $trim)};
%store2 = %{IOfastq("$input2", "fqr", $trim)};


#### stitch together reads
open(OUT1, ">$outfile");
foreach $acc (keys(%store1)){
  my $accfull = $store1{$acc}{"accfull"};  
  my $seq1 = $store1{$acc}{"seq"};
  my $qual1 = $store1{$acc}{"qual"};

  my $seq2 = $store2{$acc}{"seq"};
  my $qual2 = $store2{$acc}{"qual"};

  $seq2 =~ tr/ATCG/TAGC/; #complement
  $seq2 = reverse($seq2); # and reverse
  my $line2 =  "$seq1$seq2";

  $qual2 = reverse($qual2); # and reverse
  my $line4 =  "$qual1$qual2";

  if(length($seq1) == $trim && length($seq2) == $trim){
    print OUT1 "\@$accfull\n$line2\n+\n$line4\n";   
    }
  }
close(OUT1);

print "done.\n";



### Routine
sub IOfastq {
  ################### Usage #######################
  #################################################
  # IOfastq(infile, job, trim, outfile, index)
  #
  # infile = path to fastq infile
  #
  # job = [e][fa|fq|ac][p|r]
  # 	  [e] = optional: extract only subset of reads (to be provided in index)
  # 	  [fa|fq|ac] = returned info: fasta, fastq, accession
  # 	  [p|r] = output; either p to outfile, or r to ram, returns a reference to hash or array
  #
  # trim = length of trimming
  #
  # outfile = optional, path to outfile
  #
  # index = optional, list of reads to extract from fastq. MUST be provided as reference to indexing hash, 
  #        where the desired read names are loaded as hash's keys. e.g. $idx{"FCD0CBVABXX:1:1101:6477:2464#ATCACGAT/2"} = 1;
  # 
  # Nils Arrigo, Unil 2014 WARNING: special version with read trimming
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
  
  my ($infile, $job, $trim, $outfile, $listref) = ($_[0], $_[1], $_[2], $_[3]);
  chomp($job);
  
  ## open infile
  open(IN, "$infile");
  
  ## init values
  my %list = %{$listref};
  my %fastq;
  my %fasta;
  my @index;
  my $acc;
  my $accfull;
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
	$accfull = $acc;
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
      if($trim > 0){
	$seq = substr($seq, 0, $trim);
	}
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
      if($trim > 0){
	$qual = substr($qual, 0, $trim);
	}
      if(length($qual) != length($seq)){ #seq and qual must have the same length
	$cnt = 0;
	next;
	} else {
	
	# prepare outputs according to desired job
	if($job =~ /.*fq.*/){
	  if($job =~ /r$/){
	    $fastq{$acc}{"accfull"} = $accfull;
	    $fastq{$acc}{"seq"} = $seq;
	    $fastq{$acc}{"qual"} = $qual;
	    }
	  if($job =~ /p$/){
	    print OUT "\@$accfull\n$seq\n+\n$qual\n";
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