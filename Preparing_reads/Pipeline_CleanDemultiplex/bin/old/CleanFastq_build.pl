#!/bin/perl

##### Script pilot of Pipeline_CleanFastq.pl
### Assumes that barcodes / primers are in param/barcodes.txt and param/primers.fas files
### Input fastq must me in data.in
### Outputs produced in data.out

### Launches cleaning process on as many as possible CPUs. Upper CPU usage limit must be given

##### Usage: perl CleanFastq_build.pl maxCPUS

$CPU = $ARGV[0];


######### STEP 1 - demultiplex
## takes data from data.in and puts them in data.mult
system("mkdir -p data.mult");

# Get list of fastq files in data.in
my @list = `ls data.in/`;

# reconstruct pairs of files
my %corresp;
my %chuncks;
foreach my $i (@list){
  chomp($i);
  if($i =~ /(.*)_.*_.*_(.*)_(.*)\.fastq.gz/){
    $corresp{$3}{$2} = $i;
    }
  $chuncks{$3} = 1;
  $lane = $1;
  }

foreach my $i (keys(%chuncks)){
  my $file1 = $corresp{$i}{"R1"};
  my $file2 = $corresp{$i}{"R2"};
  
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

    my $command = "perl CleanFastq_Multiplex.pl $file1 $file2 $i &";
    print("\nRunning: $command\n\n");
    system($command);
#     sleep(60);

    } else {

    print "WAIT... CPU load is maximised\n\n\n";
#     sleep(60);
    redo;
    }
  }



######### STEP 2 - clean adapters and clip according to phred
## takes data from data.in and puts them in data.mult
system("mkdir -p data.phred");

# Get list of fastq files in data.in
my @list = `ls data.mult/`;

# reconstruct pairs of files
my %corresp;
my %chuncks;
foreach my $i (@list){
  chomp($i);
  if($i =~ /(.*)\.(.*)\.(.*)\.fastq.gz/){
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

      my $command = "perl CleanFastq_AdaptPhred.pl $file1 $file2 &";
      print("\nRunning: $command\n\n");
#       system($command);
      sleep(60);

      } else {

      print "WAIT... CPU load is maximised\n\n\n";
#       sleep(60);
      redo;
      }
    }
  }
