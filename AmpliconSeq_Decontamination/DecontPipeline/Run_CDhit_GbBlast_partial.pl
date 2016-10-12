#!/bin/perl

##### Decontamination of mutliplex PCR NGS libs. For Paired-End data
#
### Usage and parameters perl Run_CDhit_GbBlast.pl maxuse
#
# - maxuse = Maximum number of CPUs allowed to run in parallel on the server
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
## CDhit params: 
my $perc_identity_cdhit = .99; #within specimen sequence similarity [0-1]

## Blast params: first filter to eject hits. 30% min similarity, e-value 1e-4.
# Take care not using values that are too picky, you'd miss your sequences of interest. 
# Be comforted by the idea that only the best hit will be retrieved.
# my $db = "nt.09"; #GenBank Database (use nt.09 for debug and nt for full search)
my $db = "nt"; #GenBank Database (use nt.09 for debug and nt for full search)
my $perc_identity = 30; 
my $evalue = 1e-4;
my $resume = 0; #resume from current point (1 = yes, 0 = no, start fresh)?


#################
#### Actual script
# Script arguments
my $maxuse = $ARGV[0];

# hard coded: reads are assumed to be in reads/
my $readfolder = "reads/";
$readfolder =~ s/\/$//;

my $scriptname = "Run\_CDhit\_GbBlast";


# prepare working folders


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

my @todolist;
if($resume == 1){
  # check which specimens are already done
  my @done = `ls -1 data.out/*.fq\_1`;
  foreach my $file (@done){
    chomp($file);
    my $specimen = basename($file);
    $specimen =~ s/\.fq\_(.$)//;
    $specimen =~ s/^.*_//;

    if(defined($index{$specimen})){
      print "$specimen was done, skipping\n";
      delete $index{$specimen};
      }
    }
  }


### Start pipeline
my @threads;
my $startedthreads = 0;

## for debug, running only one specimen
my @todolist = keys(%index);
foreach my $spec (@todolist){

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

  if(($loadCPU + 1) < $maxuse and $freeRAM > 0) {
    #### Automated loader
    my $t = threads->new(\&runCDhit, $readfolder, $file1, $file2);
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
sub runCDhit {
  my ($readfolder, $file1, $file2) = ($_[0], $_[1], $_[2], $_[3]);

  ### step 5: extract these reads from original fastq WARNING: keywords and focal family are hard coded in bin/QualStats.r
  my $command = "R CMD BATCH \"--args intaxo='tmp/taxo/$file1.taxo' inclst='tmp/cdhit/clstr/$file1.fas.bak.clstr' inannot='tmp/taxo/$file1.annot' inhits='tmp/cdhit/lists/$file1.gb' inblast='tmp/cdhit/lists/$file1.blast' outfolder='tmp/stats' \" bin/CrossData.r";
  print("\nRunning: $command\n\n");
  system($command);

  my $command = "perl bin/FetchFastq.pl $readfolder/$file1 data.out/$file1 tmp/stats/$file1.cleanreads";
  print("$command\n");
  system($command);

  my $command = "perl bin/FetchFastq.pl $readfolder/$file2 data.out/$file2 tmp/stats/$file1.cleanreads";
  print("$command\n");
  system($command);
  }


