#!/bin/perl

#####################################
#### Get specific sequences out of fasta file, using headers from headersfile
####
####  perl FetchFasta.pl infile headersfile outfile
#####################################

my $file1 = $ARGV[0];
my $ID = $ARGV[1];
my $outfile = $ARGV[2];
my $doublons = 1; #don't care about that here, we want to keep file order intact.

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
  $acc = trim($tmp[0]);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = trim($_);
  $fasta{$acc} = $seq; 
  }

## import targets
open(HD, "$ID");
@ids = ();
my %idshs;
while (<HD>){
  chomp();
  push (@ids, $_);
  $idshs{$_} = 1;
  }
close HD;

if($doublons == 0) {
  @ids = keys(%idshs);
  } else {
  @ids = @ids;
  }

## print sequences of interest into output
open(OUT, ">$outfile");
foreach $head (@ids){
  chomp($head);
  my $seq = $fasta{$head};
  print(OUT ">$head\n$seq\n");
  }
close(OUT);


## remove whitespaces
sub trim($){
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}