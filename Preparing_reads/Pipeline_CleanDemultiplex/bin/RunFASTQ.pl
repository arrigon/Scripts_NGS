#!/bin/perl

##### Toolbox tool, to use for running FASTQC on files of interest
## Assumes standard folders architecture
## i.e. data.out/Project/stacks and data.out/Project/logs
## Just give data.out/Project as $workfolder, outputs generated in data.out/Project/logs



my $workfolder = $ARGV[0];

### STEP 4 Make quality stats (specimen level)
my @list = `ls $workfolder/stacks`;
my %corresp;

foreach my $i (@list){
  chomp($i);
  $i =~ /(.*)\.fq\_(.*)/;
  $corresp{$1}{$2} = $i;
  }

foreach my $i (keys(%corresp)){
  my $file1 = $corresp{$i}{"1"};
  my $file2 = $corresp{$i}{"2"};
  system("cat $workfolder/stacks/$file1 $workfolder/stacks/$file2 > $workfolder/logs/complete.$i.fastq");
  system("./bin/FastQC/fastqc $workfolder/logs/complete.$i.fastq");
  system("rm $workfolder/logs/complete.$i.fastq");
  }

print "Done\n";


### STEP 5 Collect these stats 
my $command = "R CMD BATCH \"--args infolder='$workfolder/logs' outfolder='$workfolder' \" bin/QualStats.r";
print("\nRunning: $command\n\n");
system($command);