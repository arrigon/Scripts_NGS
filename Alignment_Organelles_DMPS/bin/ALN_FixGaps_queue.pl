#!/bin/perl

#####################################
#### Recodes DNA gaps into DMPS alignments.
#### These alignments are typically crippled with N's (regions without sequence) AND gaps (real inserts), 
#### because most aligner give back gaps, we have to convert them into N's according to some criterion
####
#### Here, we use maxG, the longest tolerable insert length, inserts with length > maxG are converted to N's
####
#### e.g. with maxG = 6
####	>seqA
####	ATTCGC---AT------TGATGATTATGGATGATTGA------G------T----------TTGAGTCGTGATGAT
####
####	becomes
####
####	>seqA
####	ATTCGC---ATNNNNNNTGATGATTATGGATGATTGANNNNNNGNNNNNNTNNNNNNNNNNTTGAGTCGTGATGAT
####
####################################
#### usage perl ALN_FixGaps_queue.pl infile outfile maxG maxUse
####
#### -infile = fasta input (alignment)
#### -outfile = fasta output (alignment)
####
#### - maxG = longest tolerable gap length
####
#### - maxUse = Parallel run option. Launch jobs as long as current CPU load is smaller than $maxUse
####
#### NB. the script also cleans small idiosyncraties such as:
####     ~~~ -> --- (gaps manually inserted)
####     Also ensures that all sequences are starting and ending with either data (ATCG) or N
####
#### Nils Arrigo, UNIL 2014
#####################################

### Load packages
use threads;
use threads::shared;


### Script params
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $maxG = $ARGV[2];
my $maxUse = $ARGV[3];
chomp($maxUse);


### Load infile
# Open fasta and load it into a hash
open(FILE, "$infile");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

my %fasta;
my @accs; #keep track of accession input order
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = $_;
  $fasta{$tmp[0]} = $seq; 
  push(@accs, $tmp[0]);
  }

    
    
### Perform Gap checking, in parallel way
# prepare variables and pointers
my %fastaclean;
my $fastacleanref = \%fastaclean;


# launch parallel jobs
my @threads;
foreach my $acc (@accs){
  
  # Get CPU load
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + 1;

  if($load < $maxUse) {
    print "Converting gaps: $acc\n";
    my $seq = $fasta{$acc};
    my $t = threads->new(\&ConvertGaps, $acc, $seq, $maxG);
    push(@threads, $t);
    sleep(0.2);

    } else {
    sleep(0.2);
    redo;
    }
  }
 
# wait until last thread is done
foreach (@threads) {
  my $out = $_->join();
  my @fields = split("\t", $out);
  my $acc = $fields[0];
  $fastaclean{$acc} = $fields[1];
  print "Converting gaps: done with $acc\n";
  }

  

### Save to outfile
open(OUT, ">$outfile");
foreach $acc (@accs){
  my $seq = $fastaclean{$acc};
  print OUT ">$acc\n$seq\n";
  }
close(OUT);

# print "Gap checking in $infile is over\n";




##### Routines
sub ConvertGaps {  
  ### Initiate values
  my ($acc, $seq, $maxG) = ($_[0], $_[1], $_[2]);
     

  ### Step1. Make first cleaning:
  ## ~~~ -> ---
  $seq =~ tr/\~/-/;

  ## NNN -> ---
  $seq =~ tr/[N|n]/-/;

  ## --ATG -> NNATG (start of seq)
  if($seq =~ /^(-+)([^N|n])/){
    my $targ = $1;
    my $repl = $targ;
    my $RC = $2;
    $repl =~ tr/-/N/;
    $repl = "$repl$RC";
    $seq =~ s/^(-+)[^N|n]/$repl/;
    }

  ## ATG-- -> ATGNN (end of seq)
  if($seq =~ /([^N|n])(-+)$/){
    my $targ = $2;
    my $repl = $targ;
    my $LC = $1;
    $repl =~ tr/-/N/;
    $repl = "$LC$repl";
    $seq =~ s/([^N|n])(-+)$/$repl/;
    }
 
  
  ### Step2. Replace gaps with N's according to maxG
  # Initiate values
  my @stack;
  push(@stack, $seq); #working stack, where unfinished jobs go back
  my $tmp;
  my $clean;

  while($#stack >= 0){ #look for motifs as long as masking jobs need to be done
    $tmp = pop(@stack); #get sequence to check
  #   print "$tmp\n";

    if($tmp =~ /([^N|n|-]+)([N|n|-]{$maxG,})([^N|n|-]+)/g){ #check for presence of motif of type ATATNN--nnGTGA
	
#       print "hit\n";
      
      # We found a motif
      my $LN = $1; #extract and count chars on left hand of gap
      my $LNlen = length($LN);

      my $RN = $3; #extract and count chars on right hand of gap
      my $RNlen = length($RN);
      
      my $gap = $2; #extract gap and get its length as well
      my $gaplen = length($gap);
      
      
      # get actual position           
      my $start = pos($tmp) - $RNlen - $gaplen;        


      # replace gaps by N's
      my $gap = "=" x $gaplen; # prepare replacement string
      my $repl = "$LN$gap$RN";

      $tmp =~ s/([^N|n|-]+)([N|n|-]{$maxG,})([^N|n|-]+)/$repl/; #perform replacement
      
      # push back sequence into job stack
      push(@stack, $tmp);   
      
      
      } else { # the job stack is empty, stop and go have a beer.
      pop(@stack);
      }
      
    $clean = $tmp;
    }
  
  # finalize conversion
  my $seqclean = $clean;
  $seqclean =~ tr/=/N/;
  
  # return outputs
  my $out = "$acc\t$seqclean";
  
  return($out)
  }
 