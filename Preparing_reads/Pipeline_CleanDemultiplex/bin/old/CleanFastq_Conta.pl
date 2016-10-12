#!/bin/perl

##### Part of Pipeline_CleanFastq.pl
### Assumes that barcodes are in param/barcodes.txt file
### Input fastq must me in data.in. Takes fastq files (not fast.gz)
### Outputs produced in data.out

### Uses eautils softs, must be already installed on system
##### Usage: perl CleanFastq.pl fileR1 fileR2 noFile barcode

my $CPU = $ARGV[0];

######### STEP 2 - clean adapters and clip according to phred
## takes data from data.in and puts them in data.mult
system("mkdir -p data.out");
system("mkdir -p bwatmp");

## prepare BWA database
system("./bin/deconseq/bwa64 index -p conta -a bwtsw param/contaminants.fas > conta.txt 2>&1 &");
system("mv conta* bwatmp");

# Get list of fastq files in data.in
my @list = `ls data.phred/`;

# reconstruct pairs of files
my %corresp;
my %chuncks;
foreach my $i (@list){
  chomp($i);
  if($i =~ /(.*)\.(.*)\.(.*)\.fastq/){
    $corresp{$1}{$3}{$2} = $i;
    }
  }

foreach my $i (keys(%corresp)){
  my %tmp = %{$corresp{$i}};
  
  foreach my $j (keys(%tmp)){
    my $file1 = $corresp{$i}{$j}{"R1"};
    my $file2 = $corresp{$i}{$j}{"R2"};
    
    ### Get CPU load
    open PIPE, "uptime |";
    my $line = <PIPE>;
    close PIPE;
    $line =~ s/\s//g;
    my @lineArr =  split /:/, $line;
    my $times = $lineArr[@lineArr-1];
    my @timeArr = split /,/, $times;
    my $load = $timeArr[0] + 1;
    print "Current CPU load is $load \n\n\n";

    if($load < $CPU) {
      print "OK --- start $file1 and $file2\n";

#       my $command = "fastq-mcf -l 75 -o data.phred/clean.$file1 -o data.phred/clean.$file2 param/Primers.fas data.mult/$file1 data.mult/$file2 &";
      my $command = "fastq-mcf -l 75 -o data.phred/clean.$file1 param/Primers.fas data.mult/$file1 &";

      print("\nRunning: $command\n\n");
      system($command);
      sleep(60);

      } else {

      print "WAIT... CPU load is maximised\n\n\n";
      sleep(60);
      redo;
      }
    }
  }