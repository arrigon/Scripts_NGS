#!/bin/perl

#####################################
#### Removes DNA "dusts" within an alignment (as typically obtained during DMPS assemblies)
#### e.g. 
####	>seqA
####	ATTCGC---AT------TGATGATTATGGATGATTGA------G------T----------TTGAGTCGTGATGAT
####
####	becomes
####
####	>seqA
####	ATTCGC-----------TGATGATTATGGATGATTGA------------------------TTGAGTCGTGATGAT
####################################
#### usage perl ALN_DustMasking.pl infile outfile minN maxDust maxUse
####
#### -infile = fasta input (alignment)
#### -outfile = fasta output (alignment)
####
#### - minN and maxDust:
####
####  --------TGATGATTGTG--------
####  <-minN-><-MaxDust-><-minN->
####
#### -minN = Minimum number of N's (or n or -) flanking the left and right hand sides of a dust 
####	    (masking occurs if dust is flanked by AT LEAST minN N/n/- on each side -> allows N left being different of N right)
#### -maxDust = Maximum length of a dust (dusts smaller than maxDust are masked)
####
#### -maxUse = Parallel run option. Launch jobs as long as current CPU load is smaller than $maxUse
####
#### Nils Arrigo, UNIL 2014
#####################################

### Load packages
use threads;
use threads::shared;


### Script params
# my $seq = "ATGGCTGCNNNNNNRiRiNNNNNN--DoNaldDoNaldDoNald-NNNNFifiNNNNNNNGGGGNNNNNNNLouLouLouLouLouLouNNNDaisyNNNNNNNNGTTAAATNCCTTWGAATAAAATC";
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $minN = $ARGV[2];
my $maxDust = $ARGV[3];
my $maxUse = $ARGV[4];
chomp($maxUse);

### Prepare outfiles
# my $logfile = $outfile;
# $logfile = "$outfile.dustmask.txt";


# open(OUT, ">$logfile.maskdust.txt");
# print OUT "Accession\tStartPos\tLength\tMotif\tNleft\tNRight\n";


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


### Perform dust masking, in parallel way
# prepare variables and pointers
my $fastaref = \%fasta;
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
    print "Masking dusts: $acc\n";
    my $t = threads->new(\&MaskDustAcc, $acc, $fastaref);
    push(@threads, $t);
    sleep(1);

    } else {
    sleep(1);
    print "Waiting...\n";
    redo;
    }
  }
 
# wait until last thread is done
foreach (@threads) {
  my $out = $_->join();
  my @fields = split("\t", $out);
  my $acc = $fields[0];
  $fastaclean{$acc} = $fields[1];
  print "Masking dusts: done with $acc\n";
  }

  

### Save to outfile
open(OUT, ">$outfile");
foreach $acc (@accs){
  my $seq = $fastaclean{$acc};
  print OUT ">$acc\n$seq\n";
  }
close(OUT);

# print "Dust masking in $infile is over\n";


### Routines
sub MaskDustAcc {  
  my ($acc, $fastaref) = ($_[0], $_[1]);
    
  ### get variables back into function environment
  my %fasta = %{$fastaref};

  
  ### Initiate values
  my $seq = $fasta{$acc};
  my @stack;
  push(@stack, $seq); #working stack, where unfinished jobs go back
  my $tmp;
  my %clean;

  
  while($#stack >= 0){ #look for motifs as long as masking jobs need to be done
    $tmp = pop(@stack); #get sequence to check
  #   print "$tmp\n";
    
    if($tmp =~ /([N|n|-]{$minN,})([^N|n|-]{1,$maxDust})([N|n|-]{$minN,})/g){ #check for presence of motif of type NNNNNDUSTNNNNN
      
      # We found a motif
      my $LN = $1; #extract and count N's on left hand from motif
      my $LNlen = length($LN);

      my $RN = $3; #extract and count N's on right hand from motif
      my $RNlen = length($RN);
      
      my $dust = $2; #extract DNA dust and get its length as well
      my $dustlen = length($dust);
      
      
      # get actual position           
      my $start = pos($tmp) - $RNlen - $dustlen;        
    
    
      # mask motif by n's
      my $totlen = $LNlen + $dustlen + $RNlen; # get how mny n's are needed
      my $repl = "n" x $totlen; # prepare replacement string
      $tmp =~ s/([N|n|-]{$minN,})([^N|n|-]{1,$maxDust})([N|n|-]{$minN,})/$repl/; #perform replacement
      
      
      # push back sequence into job stack
      push(@stack, $tmp);   
      
      } else { # the job stack is empty, stop and go have a beer.
      # print "$acc: Dust masking is over\n";
      pop(@stack);
      }
      
    $clean{$acc} = $tmp;
    }
  my $seqclean = $clean{$acc};
  $out = "$acc\t$seqclean";
  return($out);
  }
