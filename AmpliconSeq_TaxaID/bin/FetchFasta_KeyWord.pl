#!/bin/perl

#####################################
#### Get specific sequences out of fasta file, using headers from headersfile
####
#### Usage: perl FetchFasta_KeyWord.pl infile annotfile blastfile keyword outfile
#### Nils Arrigo, Unil 2015
#####################################


my $file1 = $ARGV[0]; #fasta file, as from GeneWise outputs
my $annot = $ARGV[1]; #annotations, as from blast2lca
my $blast = $ARGV[2]; #blast hits
my $keyword = $ARGV[3];
my $mode = $ARGV[4]; #CompleteSeq or HitsOnly
my $outfile = $ARGV[5];

chomp($file1);
chomp($ID);

## import fasta file
open(FILE, "$file1");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

## load it into hash. Note that %fasta contains your sequences. Can be reused elsewhere.
my %fasta;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $acc = $tmp[0];
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = trim($_);
  $fasta{$acc} = $seq; 
  }

## import targets
open(HD, "$annot");
@ids = ();
my %idshs;
while (<HD>){
  chomp();
  if($_ =~ /$keyword/){
    my @tmp = split(/\s+/, $_, 2);
    my $acc = $tmp[0];
    push (@ids, $acc);
    $idshs{$acc} = $tmp[1];
    }
  }
close HD;


## import blast hits
open(BH, "$blast");
my $hitcnt = 0;
my %cont2coords;
while (<BH>){
  chomp();
  my @tmp = split(/\s+/, $_);
  my $acc = $tmp[0];
  my $start = $tmp[6];
  my $stop = $tmp[7];
  
  $cont2coords{$acc}{$hitcnt}{"start"} = $start;
  $cont2coords{$acc}{$hitcnt}{"stop"} = $stop;
  $hitcnt++;
  }
close BH;


## print sequences of interest into output
open(OUT, ">$outfile");
foreach $acc (@ids){
  chomp($acc);
  my $seq = $fasta{$acc};
  
  # chop out hit regions
  my @hits = keys(%{$cont2coords{$acc}});  
  my $hitcnt = 0;
  my $tmpseq;
  
  foreach my $hit (@hits){
  
    if($mode =~ /HitsOnly/){
      
      # get coords
      my $start = $cont2coords{$acc}{$hit}{"start"} - 1;
      my $stop = $cont2coords{$acc}{$hit}{"stop"} - 1;
      	  
      # prepare sequence
      my @seq = split("", $seq);
      my @tmpseq = @seq[$start..$stop];
      $tmpseq = join("", @tmpseq);
      } 
      
    if($mode =~ /CompleteSeq/){      
      $tmpseq = $seq;
      }
      
    # print to output
    print(OUT ">$acc\.hit$hitcnt\n$tmpseq\n");
    $hitcnt++;
    }
  }
close(OUT);


## remove whitespaces
sub trim($){
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub shorten($){
	my $string = shift;
	$string =~ s/^\s+.*//;
	$string =~ s/\s+.*$//;
	return $string;
}