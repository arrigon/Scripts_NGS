#!/bin/perl

use File::Basename;

### script arguments
my $haplo = $ARGV[0];
my $vcf = $ARGV[1];
my $ref = $ARGV[2];
my $outfile = $ARGV[3];
chomp($haplo);
chomp($vcf);
chomp($ref);
chomp($outfile);

my $marker = basename($vcf);
$marker =~ s/\..*$//;


### load FastPHASE results
my %haplo;
my @allspecs;
my $genoblock = 0;
my $spec;
my $cntallele;
my $cntpos;

open HAPLO, "$haplo";
while(<HAPLO>){
  my $line = $_;

  if($line =~ /BEGIN GENOTYPES/){
    $genoblock++;
    }

  if($line =~ /END GENOTYPES/){
    $genoblock++;
    }

  if($genoblock == 1){
    if($line =~ /#/){
      $spec = $line;
      $spec =~ s/#//;
      chomp($spec);
      push(@allspecs, $spec);
      $cntallele = 0;
      } else {
      my @fields = split(" ", $line);
      $cntpos = 0;
      foreach my $SNP (@fields){
	$haplo{$spec}{$cntpos}{$cntallele} = $SNP;
	$cntpos++;
	}
      $cntallele++;
      }
    }
  }
close(FAST);


### Load VCF coordinates
my %vcf;
my $cntpos = 0;

open VCF, "$vcf";
while(<VCF>){
  $line = $_;
  chomp($line);

  if($line !~ /^##/){ #we are in the VCF lines
    @fields = split("\t", $line);
    
    if($line =~ /#CHROM/){ #we are on the header line
      @accs = @fields;

      } else { #we are in data lines
      my $locus = @fields[1];
      $vcf{$cntpos} = $locus - 1;
      $cntpos++;
      }
    }
  }
close(VCF);


## Load sequence
my %refseq;
my $acc;

open REF, "$ref";
while(<REF>){
  $line = $_;
  chomp($line);

  if($line =~ /^\>.*/){
    $acc = $line;
    chomp($acc);
    $acc =~ s/\>//;
    } else {
    my @seq = split("", $line);
    my $seq_ref = \@seq;
    $refseq{$acc} = $seq_ref;
    }
  }
close(REF);


## Create fasta
open OUT, ">$outfile";

foreach my $spec (@allspecs){
  my %tmp = %{$haplo{$spec}};

  my @all1 = @{$refseq{$marker}};
  my @all2 = @{$refseq{$marker}};

  foreach my $pos (keys(%tmp)){
    my $realpos = $vcf{$pos};

    @all1[$realpos] = $tmp{$pos}{0};
    @all2[$realpos] = $tmp{$pos}{1};
    }
  
  my $all1 = join("", @all1);
  my $all2 = join("", @all2);

  print OUT ">$spec\_a\n$all1\n>$spec\_b\n$all2\n";
  }
close($outfile);

