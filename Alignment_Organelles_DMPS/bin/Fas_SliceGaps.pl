#!/bin/perl

#####################################
#### Takes a fasta and slices each sequence in independant chuncks
#### e.g. 
####	>test
####	ATTCGCGCTTAGNNNNNNNNNGATGGTTA
####	becomes
####
####	>test.chk1
####    ATTCGCGCTTAG
####	>test.chk2
####	GATGGTTA
####
####	We split the sequences at gap positions.
####
####################################
#### usage perl Fas_SplitGaps.pl infile outfile limN
#### infile = input fasta
#### outfile = output fasta
#### limN = Max number of contiguous N's that are tolerated before splitting
####
####
#### Nils Arrigo, UNIL 2013
#####################################

#load packages
use File::Basename;
use POSIX;


# Parameters
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $limN = $ARGV[2];


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
  $seq = uc($seq);
  $seq =~ s/-/N/g;
  $fasta{$tmp[0]} = $seq; 
  push(@accs, $tmp[0]);
  }



# Parse that hash, accession by accession (following original order),
# look for sequence chuncks coordinates (start pos)
my %coords;
foreach my $acc (@accs){ #iterate accession
  my $cnt = 0;
  my $seq = $fasta{$acc}; #get sequence
  
  my @seq = split("", $seq);
  my $len = $#seq;

  my $startflag = 0;
  foreach $i (0..$len){
    if($seq[$i] !~ /N/){
      $Ncntdown = 0;
      if($startflag == 0){
	$startflag = 1;
	$cnt++;
	$coords{$acc}{$cnt}{"start"} = $i;
	}
      if($startflag == 1 && $i == $len){
	$coords{$acc}{$cnt}{"stop"} = $i;
	}
      } else {
	$Ncntdown++;
	if($startflag == 1 && $Ncntdown > $limN){
	$startflag = 0;
	$coords{$acc}{$cnt}{"stop"} = $i - 1;        
	}
	if($startflag == 1 && $i == $len){
	$startflag = 0;
	$coords{$acc}{$cnt}{"stop"} = $i - 1;        
	}
      }
    }  
  }



# open output file (update its name with bp positions)
my $bsn = basename($infile);
$bsn =~ s/\..*$//;
open(OUT, ">$outfile");

foreach $acc (@accs){
  my $seq = $fasta{$acc}; #get sequence
  my @seq = split("", $seq);
  my @chuncks = keys(%{$coords{$acc}});
  
  @chuncks = sort{$a <=> $b} @chuncks;
  foreach $chk (@chuncks){
    my $start = $coords{$acc}{$chk}{"start"};
    my $stop = $coords{$acc}{$chk}{"stop"};
    my @seqchk = @seq[$start..$stop];
    my $seqchk = join("", @seqchk);
    $seqchk =~ s/N//g;
    if(length($seqchk) > 10){
      print OUT ">$acc.chk$chk\n$seqchk\n"
      }
    }
  }



close(OUT); 
