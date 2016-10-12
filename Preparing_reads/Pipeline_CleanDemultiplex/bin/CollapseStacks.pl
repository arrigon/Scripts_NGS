#!/bin/perl
# Run next to stacks folder, aims at merging fastq files 

system("mkdir stacksOK/");
my @allfiles = `ls stacks/`;

my %index;
foreach $file (@allfiles){
  chomp($file);
  $oldname = $file;
  $file =~ /(.*\_.*\_.*)_(\d*)\.fq\_(.*)/;
  
  my $ID = $1;
  my $lane = $2;
  my $read = $3;
  
  $index{$ID}{$lane}{$read} = $oldname;
  }



foreach my $ID (keys(%index)){
  my %tmp = %{$index{$ID}};
  foreach my $lane (keys(%tmp)){
    my $file1 = $index{$ID}{$lane}{"1"};
    my $file2 = $index{$ID}{$lane}{"2"};
    my $command = "cat stacks/$file1 >> stacksOK/$ID.fq\_1";
    print "$command\n";
    system($command);

    my $command = "cat stacks/$file2 >> stacksOK/$ID.fq\_2";
    print "$command\n";
    system($command);
    }
  }
