#!/bin/perl
# Run next to stacks folder, aims at merging fastq files 

my $dir1 = $ARGV[0];
my $dir2 = $ARGV[1];

system("mkdir readDMPSMerged/");


# load dir1
my @allfiles = `ls $dir1`;
my %index;
foreach $file (@allfiles){
  chomp($file);
  $oldname = $file;
  $file =~ /(.*)\.fq\_(.*)/;
  
  my $ID = $1;
  my $read = $2;
    
  $index{$ID}{$read} = $oldname;
  }


# load dir2
my @allfiles = `ls $dir2/`;
my %index;
foreach $file (@allfiles){
  chomp($file);
  $oldname = $file;
  $file =~ /(.*)\.fq\_(.*)/;
  
  my $ID = $1;
  my $read = $2;
  
  $index{$ID}{$read} = $oldname;
  }



foreach my $ID (keys(%index)){
  my $file1 = $index{$ID}{"1"};
  my $file2 = $index{$ID}{"2"};
  my $command = "cat $dir1/$file1 >> readDMPSMerged/$ID.fq\_1";
  print "$command\n";
  system($command);

  my $command = "cat $dir1/$file2 >> readDMPSMerged/$ID.fq\_2";
  print "$command\n";
  system($command);

  my $command = "cat $dir2/$file1 >> readDMPSMerged/$ID.fq\_1";
  print "$command\n";
  system($command);

  my $command = "cat $dir2/$file2 >> readDMPSMerged/$ID.fq\_2";
  print "$command\n";
  system($command);
  }
