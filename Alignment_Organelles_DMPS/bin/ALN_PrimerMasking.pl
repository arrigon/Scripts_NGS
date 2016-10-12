#!/bin/perl

#####################################
#### Mask primer positions (e.g. hits to primers), defined from an external position file
#
#### Position file : typical blastn tabular output
# Primer	Start	Stop
# TM-N200	27	49
# TI-J34	76	101
# C2-J-3696	488	510
# TW-J1301	568	586
#
## WARNING: this script assumes to work on a sequence alignment
#
## Outputs
# Produces outfile and outfile.hiddenpos.fas, where the first is ready to use (hiddenchar = N) 
# while the other outfile shows how the job was done (hiddenchar = !)
#
#### Usage: perl ALN_PrimerMasking.pl infile outfile posfile
#####################################

### load packages
use POSIX;
use List::Util qw[min max];


### Script params
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my $posfile = $ARGV[2];
my $hidechar = "!"; #adjust hidding character according to your needs
chomp($posfile);


## import fasta file, from infile
open(FILE, "$infile");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

## load it into hash. Note that %fasta contains your sequences. Can be reused elsewhere.
my %fasta;
my @accs; #keep track of accession input order
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $seq = $_;
  my @seq = split("", $seq);
  $fasta{$tmp[0]} = \@seq;
  
  push(@accs, $tmp[0]);
  }


## import positisions to hide, from posfile
# and modify directly into fasta store
my %pos;
open(POS, "$posfile");
while(<POS>){
  chomp();
  my $line = $_;
  my @fields = split("\t", $line);

  my @startstop = @fields[1..2]; #adjust here to extract infos from other position files
  
  foreach $head (@accs){ #loop over all accessions
    my $start = min(@startstop);
    my $stop = max(@startstop);

    # replace target regions with hidding character
    my @bp = @{$fasta{$head}};
    foreach $hide ($start..$stop){
      $bp[$hide] = "$hidechar";
      }
    
    # store back to %fasta
    my $bpref = \@bp;
    $fasta{$head} = $bpref;
    }
  }
close(POS);


## Print sequences to outfile
open(OUT, ">$outfile");
$outfile =~ s/\..*$//;
open(LOG, ">$outfile.hiddenpos.fas");
foreach $head (@accs){
  chomp($head);
  my @bp = @{$fasta{$head}};
  my $seq = join("", @bp);
  print(LOG ">$head\n$seq\n");

  $seq =~ s/$hidechar/N/g;
  print(OUT ">$head\n$seq\n");
  }
close(LOG);
close(OUT);


## remove whitespaces
sub trim($){
  my $string = shift;
  $string =~ s/^\s+//;
  $string =~ s/\s+$//;
  return $string;
  }