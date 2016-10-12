#!/bin/perl

#### Running mitoBIM, using a reference + SE libraries
# Nils Arrigo, Unil 2014
####


use threads;
use threads::shared;
use File::Basename;


#### Get script arguments
my $list = $ARGV[0];
my $MaxUse = $ARGV[1];
my $scriptname = "mitoBIM loader";
my $iter = 500;
my $outfile = "mitoBIMAssembly\_iter$iter.fasta";
my $paired = 1;
my $kbait = 31;


#############################################
## start fresh and prepare projects folder
my $command = "rm -rf data.out tmp/";
print "### $scriptname : $command\n";
system("$command");

my $command = "mkdir -p data.out tmp/ tmp/reads tmp/assemblies";
print "### $scriptname : $command\n";
system("$command");

## parse whishlist
my %store;
open(IN, $list);
while(<IN>){
  chomp();
  my @fields = split("\t", $_);
  my $read = $fields[0];
  my $ref = $fields[1];
  $store{$read} = $ref;
  }
  
## parse reads
my @reads = keys(%store);


##### Iterate over all datasets
my %threads;
my $startedthreads = 0;

foreach $read (@reads){
# $read = $reads[0];
  chomp($read);

  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $load = $timeArr[0] + $CPU;
  
  if($load < $MaxUse) {
    ## prepare working folder
    my $bsn = basename($read);
    $read = $bsn;
    $bsn =~ s/\.f.*$//;

    my $command ="mkdir tmp/$bsn";
    print "### $scriptname : $command\n";
    system("$command");

    ## get reference
    my $ref = $store{$read};
    
    
    if($paired == 0){
      my $command ="cp reads/$read tmp/reads/$bsn.fastq";
      print "### $scriptname : $command\n";
      system("$command");
      } else {
      my $command ="cat reads/$bsn.fq\_1 reads/$bsn.fq\_2 > tmp/reads/$bsn.fastq";
      print "### $scriptname : $command\n";
      system("$command");
      }


    my $command ="cp bin/MITObim_1.7_custom.pl tmp/$bsn/";
    print "### $scriptname : $command\n";
    system("$command");

    ## move there and run mitoBIM
    chdir("tmp/$bsn/");
    my $t = threads->new(\&runmitoBIM, $scriptname, $read, $ref, $paired);
    $threads{$bsn} = $t;
    sleep(5);
    $startedthreads++;

    ## come back to root dir
    chdir("../../");

    } else {

    sleep(5);
    redo;
    }
  }

my $cnt = 0;
foreach my $bsn (keys(%threads)) {
  my $t = $threads{$bsn};
  my $genome = $t -> join;
  $cnt++;
  print "Finished $bsn (job $cnt / $startedthreads)\n";

  ## parse outputs
  my @folders = `ls tmp/$bsn/`;
  my @nrs;
  foreach my $fld (@folders){
    chomp($fld);
    if($fld =~ /iteration(\d*)/){
      my $nr = $1;
      push(@nrs, $nr);
      }
    }
  my $last = max(@nrs);
#   print join("\n", @folders);

  my $command ="cp tmp/$bsn/iteration$last/$bsn-Seed\_assembly/$bsn-Seed\_d\_results/*$bsn.unpadded.fasta tmp/assemblies/$bsn.fasta";
  print "### $scriptname : $command\n";
  system("$command");

  my $command ="cp tmp/$bsn/$bsn.log tmp/assemblies/";
  print "### $scriptname : $command\n";
  system("$command");

  ## clean mess and go back to root folder
#   my $command ="rm -rf tmp/$bsn/ tmp/reads/$bsn*";
#   print "### $scriptname : $command\n";
#   system("$command");  
  }


### Collect fasta files
my $command ="perl bin/CleanCollectFasta.pl tmp/assemblies data.out/$outfile";
print "### $scriptname : $command\n";
system("$command");



### Routines
sub runmitoBIM {
  ## get arguments
  my ($scriptname, $read, $ref, $paired) = ($_[0], $_[1], $_[2], $_[3]);
  my $bsn = basename($read);
  $bsn =~ s/\.f.*$//;

  chomp($ref);
  my $ref = basename($ref);

  ## run MitoBIM
  if($paired == 0){
    my $command ="./MITObim_1.7_custom.pl -end $iter -sample $bsn -ref Seed -readpool ../reads/$bsn.fastq --clean --trimreads --quick ../../refs/$ref --kbait $kbait > $bsn.log";
    print "### $scriptname : $command\n";
    system("$command");
    } else {
    ## run MitoBIM
    my $command ="./MITObim_1.7_custom.pl -end $iter -sample $bsn -ref Seed -readpool ../reads/$bsn.fastq --clean --pair --trimreads --quick ../../refs/$ref --kbait $kbait > $bsn.log";
    print "### $scriptname : $command\n";
    system("$command");
    }
  }

sub max ($$) { $_[$_[0] < $_[1]] }


