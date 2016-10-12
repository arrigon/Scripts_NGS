#!/bin/perl

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

open(IN, "$infile");

my %store;
my @list;

while(<IN>){
  chomp();
  my $line = $_;
  my @fields = split(" ", $line);
  my $gbhit = $fields[1];
  my $readhit = $fields[0];
  $gbhit = "$readhit $gbhit";

  if(!defined($store{$gbhit})){
    $store{$gbhit} = $line;
    push(@list, $gbhit);
    } else {
#     print "skip\n";
    }
  }
close(IN);

open(OUT, ">$outfile");
foreach my $gbhit (@list){
  my $line = $store{$gbhit};
  print(OUT "$line\n");
  }
close(OUT);

# print "done\n";