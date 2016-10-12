#!/bin/perl

use File::Basename;

### script arguments
my $infile = $ARGV[0];
my $outfile = $ARGV[1];
chomp($infile);
chomp($outfile);


### Load VCF into memory
open(IN, "$infile");
my @allpos;
my %store;

while(<IN>){
  $line = $_;
  chomp($line);

  if($line !~ /^##/){ #we are in the VCF lines
    @fields = split("\t", $line);
    
    if($line =~ /#CHROM/){ #we are on the header line
      @accs = @fields;

      } else { #we are in data lines
      my $locus = @fields[1];
      push(@allpos, $locus);

      # store infos about locus
      my @bp = @fields[3..4];

      # store genotypes
      foreach $i (9..$#fields){
	my $spec = basename($accs[$i]);
	$spec =~ s/Rad$RAD\.//g;
	$spec =~ s/\.bam//g;

	my $genot = $fields[$i];
	my @data = split(":", $genot);
  
	if($data[2] > 3){
	  my @tmp = split("/", $data[0]);
	  my $all1 = $bp[$tmp[0]]; 
	  my $all2 = $bp[$tmp[1]];
	  $store{$spec}{$locus}{"all1"} = $all1;
	  $store{$spec}{$locus}{"all2"} = $all2;
	  } else {
	  $store{$spec}{$locus}{"all1"} = "?";
	  $store{$spec}{$locus}{"all2"} = "?";
	  }
	} #end of specimen parsing
      }  #end of data parsing
    }  #end of if !~ /^##/
  } #end of data importing
print "Loading VCF...\n";

my @allspec = keys(%store);


### Convert to PHASE input file
open(OUT, ">$outfile");

my $nspec = $#allspec + 1;
my $npos = $#allpos + 1;

print OUT "$nspec\n$npos\n";
print OUT "P ", join(" ", @allpos), "\n";
print OUT join(" ", ("S" x $npos)), "\n";

foreach my $spec (keys(%store)){
  my @all1;
  my @all2;
  
  foreach my $locus (@allpos){
    push(@all1, $store{$spec}{$locus}{"all1"});
    push(@all2, $store{$spec}{$locus}{"all2"});
    }

  my $all1 = join(" ", @all1);
  my $all2 = join(" ", @all2);

  print(OUT "#$spec\n$all1\n$all2\n");
  }
close(OUT);

print "test\n";




