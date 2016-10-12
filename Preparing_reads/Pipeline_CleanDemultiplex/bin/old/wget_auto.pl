#!/bin/perl 

$file = $ARGV[0];

open(IN, "$file");
while(<IN>){
  chomp();
  $url = $_;
  system("wget $url");
  }