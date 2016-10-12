#!/bin/perl


#### Get script arguments
my $scriptname = "MiraAssembly";


####
my $command = "mkdir -p data.out/ tmp/";
print "### $scriptname : $command\n";
system("$command");


####
foreach my $i (1..4){
  my $command = "perl bin/PrepareManifestMito.pl $i tmp/manifest.$i.conf";
  print("$command\n");
  system($command);

  my $command = "mira tmp/manifest.$i.conf > tmp/assembly.$i.log";
  print("$command\n");
  system($command);
                    
  my $command = "cp Assembly$i\_assembly/Assembly$i\_d\_results/Assembly$i\_out\_Mito$i.unpadded.fasta data.out/assembly.$i.fasta";
  print("$command\n");
  system($command); 

  my $command = "cp Assembly$i\_assembly/Assembly$i\_d\_info/Assembly$i\_info\_assembly.txt data.out/assembly.$i.info";
  print("$command\n");
  system($command); 

  my $command = "rm -rf Assembly$i\_assembly/";
  print("$command\n");
  system($command);
  }





