#!/bin/perl

#####################################
#### Cross blast, annot, taxo and index tables to produce metadata
####
#### Usage: perl CrossParseAnnot.pl blast taxo annot idx infasta outfasta
#### Nils Arrigo, Unil 2015
#####################################


my $blast = $ARGV[0]; #blast file, as from GeneWise outputs
my $taxo = $ARGV[1]; #taxo file
my $annot = $ARGV[2]; #annot file
my $idx = $ARGV[3]; #annot file
my $file1 = $ARGV[4]; #fasta where to get headers
my $outfile = $ARGV[5]; #outfile
chomp($outfile);

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
my @accs;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $acc = $tmp[0];
  $acc =~ s/\.hit.*$//g;
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = trim($_);
  $fasta{$acc} = $seq; 
  push(@accs, $acc);
  }

  
### Reconstruct annotations  
## import Blast hits
open(BH, "$blast");
my %cont2gid;
my %gid2cont;
while (<BH>){
  chomp();
  my @tmp = split(/\s+/, $_, 2);
  my $acc = $tmp[0];
  
  my $gid = $tmp[1];
  $gid =~ /gi\|.*\|.*\|(.*)\|/;
  $gid = $1;
  
  $cont2gid{$acc}{"gid"} = $gid;
  $gid2cont{$gid} = $acc;
  }
close BH;


## import Taxo
open(TA, "$taxo");
while (<TA>){
  chomp();
  my @tmp = split(/\s+/, $_, 2);
  my $acc = $tmp[0]; 
  my $info = $tmp[1];
  
  $cont2gid{$acc}{"taxo"} = $info;
  }
close TA;


## import Annots
open(AN, "$annot");
while (<AN>){
  chomp();
  my @tmp = split(/\s+/, $_, 2);
  my $gid = $tmp[0];
  
  my $info = $tmp[1];
  $gid2cont{$gid} = $info;
  }
close AN;


## import Idx
open(ID, "$idx");
while (<ID>){
  chomp();
  my @tmp = split(/\s+/, $_, 2);
  my $acc = $tmp[0];
  my $head = $tmp[1];
  $cont2gid{$acc}{"head"} = $head;
  }
close ID;


## print sequences of interest into output
open(OUT, ">$outfile");
foreach $acc (@accs){
  chomp($acc);
  my $taxo = $cont2gid{$acc}{"taxo"};
  my $gid = $cont2gid{$acc}{"gid"};
  my $head = $cont2gid{$acc}{"head"};
  my $annot = $gid2cont{$gid};
  
  print(OUT "$acc\t$head\t$gid\t$annot\t$taxo\n");
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