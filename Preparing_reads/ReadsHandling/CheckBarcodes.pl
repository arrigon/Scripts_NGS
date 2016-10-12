#!/bin/perl

##### Generate some sanity check stats from a fastq file
# 
# Input: $infile
#
# Usage: perl CheckFileNames.pl infolder outfolder
#
# WARNING: this script is based on hard-coded RegExp, see comments below
########################################################

use File::Basename;

### Script params
my $infile = $ARGV[0];
my $lenbcde = 6;
my $lenrstr = 6;
my $bsn = basename($infile);
$bsn =~ s/\.gz$//;

# 
print "Sanity checks for $bsn\n-barcode = $lenbcde bp\n-rest.site = $lenrstr bp\n";

# prepare working folder
system("mkdir -p CheckBarcodes/");

# uncompress file
print "Uncompressing $bsn...\n";
my $command = "zcat -c $infile > CheckBarcodes/$bsn";
system($command);

my $command = "grep -c TGCAGG CheckBarcodes/$bsn > $bsn.restr.instances ";
system($command);


# load it into RAM
print "Loading $bsn in RAM...\n";
my $ref = IOfastq("CheckBarcodes/$bsn", "fqr");
my %fastq = %{$ref};

# Parse fastq
print "Parsing file...\n";
my %bcdes;
my %rstr;
foreach $acc (keys(%fastq)){
  my $seq = $fastq{$acc}{"seq"};
  my $barcode = substr($seq, 0, $lenbcde);
  my $restr = substr($seq, $lenbcde, $lenrstr);

  if($bcdes{$barcode}){
    $bcdes{$barcode} = $bcdes{$barcode} + 1;
    } else {
    $bcdes{$barcode} = 1;
    }

  if($rstr{$restr}){
    $rstr{$restr} = $rstr{$restr} + 1;
    } else {
    $rstr{$restr} = 1;
    }
  }

# Producing outputs
open(OUT, ">$bsn.barcodes.txt");
foreach my $bcde (sort { $bcdes{$b} <=> $bcdes{$a} } keys %bcdes) {
  my $occs = $bcdes{$bcde};
  print OUT "$bcde\t$occs\n";
  }
close OUT;

open(OUT, ">$bsn.restr.txt");
foreach my $restr (sort { $rstr{$b} <=> $rstr{$a} } keys %rstr) {
  my $occs = $rstr{$restr};
  print OUT "$restr\t$occs\n";
  }
close OUT;

# clean mess
# my $command = "rm -rf CheckBarcodes";  
# system($command);

print "Done. Outputs in $bsn.barcodes.txt and $bsn.restr.txt\n";



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
