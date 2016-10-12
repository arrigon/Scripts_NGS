#!/bin/perl

##### Assembly of PE data, from PCR amplicons.
#
#
### Usage and parameters perl Run_usearch.pl maxuse
#
# - maxuse = Maximum number of CPUs allowed to run in parallel on the server
#
# we assume that reads are stored in reads/
#
# for the moment, results are landing in tmp/postfiltering (fasta files and blast searches)
# 
###
# Working scripts / external programs: must be in bin/ folder
#
### Nils Arrigo, Ephemeroptera project 2014.
########################################################

#### libraries
use threads;
use threads::shared;
use File::Basename;
use Sys::MemInfo qw(totalmem freeswap freemem);
use POSIX;

#### Hard coded params
## Usearch command
my $usearch = "usearch7.0.1090_i86linux32";

## Clustering params: similarity threshold
my $perc_identity = 1; # 3% is Uclust default, I'd use something more stringent.
my $trim_len = 280; # For figure cases where F and R reads are not overlapping, we'll append them
		    # yet, this requires to trim F and R reads at the same length
		    # this value depends on your sequencing platform.

## Blast params: just to have an idea of what you got
# my $db = "nt.09"; #GenBank Database (use nt.09 for debug and nt for full search)
my $db = "nt"; #GenBank Database (use nt.09 for debug and nt for full search)
my $evalue = 1e-4;
my $blast_perc_identity = 30;

## Filtering params
my $mincov = 3; # minimum coverage that must be guaranteed
my $minfrac = .33; # coverage must be at least $minfrac * coverage of most abundant sequence


#################
#### Actual script
# Script arguments
my $maxuse = $ARGV[0];
my $readfolder = "reads/";
my $reffolder = "refs/";

my $scriptname = "Run\_Usearch\t";


# prepare working folders
system("rm -rf tmp data.out");
system("mkdir -p tmp tmp/fq tmp/fas tmp/cls tmp/ublast tmp/postfilter tmp/merged data.out");


# index all reads per specimens
my @allfiles = `ls -1 $readfolder`;
my %index;
foreach my $file (@allfiles){
  chomp($file);
  my $specimen = $file;
  $specimen =~ s/\.fq\_(.$)//;
  my $readnr = $1;
  $index{$specimen}{$readnr} = $file;
  }


### Start pipeline
my @threads;

my $startedthreads = 0;
# foreach my $spec (keys(%index)){

## for debug, running only one specimen
my @dddd = keys(%index);
# @dddd = @dddd[0];
foreach my $spec (@dddd){

  my $file1 = $index{$spec}{1};
  my $file2 = $index{$spec}{2};

  # inspect ram before firing up
  my $freeRAM = ceil((&freemem) / 1073741824);
  
  # inspect CPU
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $loadCPU = $timeArr[0] + $CPU;


  print "CPU load = $loadCPU\tFree memory = $freeRAM Gb\n";

  if(($loadCPU + 2) < $maxuse and $freeRAM > 0) {
    #### Automated loader
    my $t = threads->new(\&runUsearch, $readfolder, $file1, $file2);
    push(@threads, $t);
    $startedthreads++;
    sleep(1);

    } else {
    sleep(1);
    redo;
    }
  }


my $cnt = 0;
foreach (@threads) {
  my $genome = $_-> join;
  $cnt++;
  print "Finished job $cnt / $startedthreads\n";
  }

print "Done.\n";



### Routines
sub runUsearch {
  my ($readfolder, $file1, $file2) = ($_[0], $_[1], $_[2], $_[3]);
    
  my $bsn = $file1;
  $bsn =~ s/\.fq\_1//;

  ### Step 1. Preprocess reads
  # merge forward and reverse reads
  my $command = "$usearch -fastq_mergepairs $readfolder/$file1 -reverse $readfolder/$file2 -fastq_truncqual 3 -fastq_maxdiffs 0 -fastqout tmp/fq/$bsn.merged.fq";
  print("$command\n");
  system($command);

  my $test = `wc -l tmp/fq/$bsn.merged.fq`;
  if($test == 0){
    my $command = "perl bin/StitchReads.pl $readfolder/$file1 $readfolder/$file2 $trim_len tmp/fq/$bsn.merged.fq";
    print("$command\n");
    system($command);
    }

  # convert to fasta
  my $command = "perl bin/Fastq2Fasta.pl tmp/fq/$bsn.merged.fq tmp/fas/$bsn.merged.fas";
  print("$command\n");
  system($command);

  # check read orientations
  my $command = "head -n 2 tmp/fas/$bsn.merged.fas > tmp/ublast/$bsn.ref";
  print("$command\n");
  system($command);

  my $command = "$usearch -makeudb_ublast tmp/ublast/$bsn.ref -output tmp/ublast/$ref.udb";
  print("$command\n");
  system($command);

  my $command = "$usearch -ublast tmp/fas/$bsn.merged.fas -db tmp/ublast/$ref.udb -evalue 1e-5 -strand both -uc tmp/ublast/$bsn.ublast";
  print("$command\n");
  system($command);
  
  my $command = "perl bin/OrientFasta.pl tmp/fas/$bsn.merged.fas tmp/fas/$bsn.orient.fas tmp/ublast/$bsn.ublast";
  print("$command\n");
  system($command);
  
  # dereplicate full-length mode #seems to work better than prefix mode
  my $command = "$usearch -derep_fulllength tmp/fas/$bsn.orient.fas -output tmp/fas/$bsn.uniques.fas -sizeout -minseqlength 100";
  print("$command\n");
  system($command);

  my $command = "$usearch -derep_prefix tmp/fas/$bsn.uniques.fas -output tmp/fas/$bsn.uniques.fas -sizeout -minseqlength 100";
  print("$command\n");
  system($command);

  # Sort by cluster size (WARNING, this not exactly the coverage that is used here, leave to minsize = 0)
  my $command = "$usearch -sortbysize tmp/fas/$bsn.uniques.fas -output tmp/fas/$bsn.sorted.fas -minsize 0";
  print("$command\n");
  system($command);

  ### Step 2: cluster sequences
  my $command = "$usearch -cluster_otus tmp/fas/$bsn.sorted.fas -otu_radius_pct $perc_identity -otus tmp/cls/$bsn.markers";
  print("$command\n");
  system($command);

  ### Step 3: map back reads to markers and get coverages
  my $id = (100 - $perc_identity) / 100;
  my $command = "$usearch -usearch_global tmp/fas/$bsn.merged.fas -db tmp/cls/$bsn.markers -strand both -id $id  -uc tmp/cls/$bsn.clusters";
  print("$command\n");
  system($command);

  ### Step 3: cluster sequences
  my $command = "perl bin/FetchCovs.pl tmp/cls/$bsn.markers tmp/cls/$bsn.clusters tmp/fas/$bsn.clustered.fas";
  print("$command\n");
  system($command);

  ### Step 4: filter out best candidates
  my $command = "perl bin/FilterFasCov.pl tmp/fas/$bsn.clustered.fas tmp/postfilter/$bsn.fas $mincov $minfrac";
  print("$command\n");
  system($command);

# #   ### Step 5: produce blast report
#   my $command = "blastn -query tmp/postfilter/$bsn.fas -db $db -task blastn -dust yes -num_threads 1 -evalue $evalue -perc_identity $blast_perc_identity -out tmp/postfilter/$bsn.blast -max_target_seqs 10 -outfmt \"7 qacc sseqid length pident evalue sscinames\" ";
#   print("$command\n");
#   system($command);

  #clean mess
#   my $command = "rm tmp/fq/$bsn.merged.fq tmp/fas/$bsn.merged.fas tmp/fas/$bsn.uniques.fas tmp/fas/$bsn.sorted.fas tmp/cls/$bsn.markers";
#   print("$command\n");
#   system($command);
  }


