#!/bin/perl
####################################
# Merging VCF files, and ensure proper treatment of columns containing different specimens
#
# WARNING: assumes that rads are named using the following format (in CHROM field): "Rad128988"
#
### Params
# - infolder = input folder with VCF files to process
# - outfile = where to save outputs
#
### Usage: perl MergeVCF.pl infolder outfolder
#
# Nils Arrigo. Uni Lausanne 2012
############################################################
use File::Basename;

### Get script arguments
my $infolder = $ARGV[0];
my $outfile = $ARGV[1];

chomp($infolder);
chomp($outfile);


### Open outfile
open(OUT, ">$outfile");


### Parse infolder and get all vcf in there
my @list = `ls $infolder/*.vcf`;
print "Loading of VCF files...\r";

### load each VCF in %store
my %store;
my %allspecimens;
my $filecnt = 0;

foreach my $infile (@list){
  chomp($infile);

  # get rad name
  my $RAD = basename($infile);
  $RAD =~ s/\.vcf//g;
  $RAD =~ s/Rad//g;

  # parse file
  open(IN, "$infile");
  my @accs;

  while(<IN>){
    $line = $_;
    chomp($line);

    if($line =~ /^##/){ #we are in the comments lines
      if($filecnt == 0){ #print the captions in outfile, do so only during visiting first file
	print OUT "$line\n";
	}
      } else { #we are in the VCF lines
      @fields = split("\t", $line);
      
      if($line =~ /#CHROM/){ #we are on the header line
	@accs = @fields;

	} else { #we are in data lines
	my $locus = @fields[1];
	
	# store infos about locus
	my @info = @fields[0..8];
	$store{$RAD}{$locus}{"info"} = \@info;

	# store genotypes
	foreach $i (9..$#fields){
	  my $spec = basename($accs[$i]);
	  $spec =~ s/Rad$RAD\.//g;
	  $spec =~ s/\.bam//g;
	  my $genot = $fields[$i];
	  $store{$RAD}{$locus}{$spec} = $genot;
	  $allspecimens{$spec} = $RAD;
	  } #end of specimen parsing
	}  #end of data parsing
      }  #end of if !~ /^##/
    } #end of data importing
  $filecnt++;
  } #end of file browsing
print "Loading of VCF files... done\n";


### List all rads and all specimens
my @specimens = keys(%allspecimens);
my @specimens = sort {$a cmp $b} @specimens;
my @rads = keys(%store);
my @rads = sort {$a <=> $b} @rads;

print "Outputting to $outfile...\r";


### Parse the store and reconstruct a complete VCF file
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
print OUT join("\t", @specimens);
print OUT "\n";

foreach my $RAD (@rads){ #visit each marker
  my @locuspos = keys(%{$store{$RAD}});
  my @locuspos = sort {$a <=> $b} @locuspos;

  foreach $locus (@locuspos){ #in each marker, visit each position

    # get metadata about locus
    my @info = @{$store{$RAD}{$locus}{"info"}};
    print OUT join("\t", @info), "\t";    

    # get genotypes
    foreach my $spec (@specimens){ #in there, visit each possible specimen
    if(my $genot = $store{$RAD}{$locus}{$spec}){ #the marker is present in the focal specimen
      print OUT "$genot\t";
      } else { #that marker is not present in a given specimen
      print OUT "0/0:0,0,0:0:0:0\t";
#       print OUT "Missing\t";
      }
      }
    print OUT "\n";
    }
  }
close(OUT);
print "Outputting to $outfile... done\n";


