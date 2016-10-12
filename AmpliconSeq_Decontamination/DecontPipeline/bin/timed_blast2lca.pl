#!/usr/bin/perl -w

######################################################################################
# This is quick wrapper to launch blast2lca with a timer... 
# this program hangs stochastically and we need to monitor how long it take, and relaunch it when needed
######################################################################################

#### Get script inputs
my $file1 = $ARGV[0];
my $timeout = $ARGV[1];
chomp($timeout);


$SIG{ALRM} = sub {
  # make sure that we kill this damn job!! strong memory leak otherwise
  my $cmd = "kill -9 \$(ps aux \| grep tmp/cdhit/lists/$file1.blast \| awk -F \' \' \'{print \$2}\')";
  system($cmd);
  print STDOUT "timedout";
  exit(1); };

# make sure we start the timer
alarm($timeout);

# actual job
my $command = "./bin/blast2lca -names refs/names.dmp -nodes refs/nodes.dmp -dict refs/gi_taxid.bin -levels=superkingdom:kingdom:phylum:order:family -order=true tmp/cdhit/lists/$file1.blast --savemem > tmp/taxo/$file1.taxo";
system($command);

# if we completed the task, print "completed" to STDOUT (will be catched by other threads to move on)
print STDOUT "completed";
