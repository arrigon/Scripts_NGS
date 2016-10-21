#!/bin/perl

use File::Find;
use File::Basename;

### Script params
my $folder = $ARGV[0];
my $feature = "CDS";
my $outfolder = "CDS";


### Produce sink folder
system("mkdir -p $outfolder");

### Listing files in target folder
# get list of GFFs
find( sub {push @files, $File::Find::name if (/\.gff$/)},$folder);

# my @files = `ls $folder`;

# move to target folder
# chdir($folder);


# keep count of files, genes and features that are visited / exported
my $filecnt = 0;
my $totgenecnt = 0;
my $featcnt = 0;
my $totfeatcnt = 0;

### Parsing GFF files
# loop over GFFs
print "##########\nParsing GFF files...\n##########\n\n";

foreach $gff (@files){
  # iterate through each GFF file
  
  if($gff =~ /[\.gff]$/ ){ #open only if has *.gff extension
    chomp($gff);
      
    #increase file counter
    $filecnt++; 

    # open input / output files
    print "#####\nParsing $gff...\n";
    open(IN, $gff);

    # keep track of where we are at in GFF file, set fastablock to 1 as soon as we enter into the sequence data
    my $fastablock = 0;
    my $skip = 0; #we have to skip the two first lines in the FASTA block, so needs a counter to be initiated here

    # numbering of genes found within file
    my $localgenecnt = 0;
    my $localfeatcnt = 0;

    # hash for storing annotation infos
    my %annotations;

    # Parse lines of gff file
    while(<IN>){
      # iterate through each line
      chomp();
      $line = $_;
      
      # check where we are at into the GFF file
      if($line =~ /##FASTA/){
	$fastablock = 1;
	}

      ### GFF annotation block
      ## if we are still in the annotation block (i.e. fastablock == 0), look for features and store them into annotation hash
      if($fastablock == 0){      
	
	# track gene starts into GFF block
	if($line =~ /maker\sgene/){
	  $totgenecnt++;
	  $localgenecnt++;
	  $featnum = 0; #reinitiate features numbering

	  # initiate start and stop arrays
	  my @starts = undef;
	  my @stops = undef;
	  }

	# track features into GFF block
	if($line =~ /maker\s$feature/){

	  # keep track of features numbering; important since CDS must be spliced in that order out of sequence
	  $totfeatcnt++;
	  $localfeatcnt++;

	  # split line and retrieve annotation data (scaff name, start-stop positions)
	  @array = split(/\t/, $line);
	  $scaff = $array[0];

	  # Store data into annotation hash, will be used when parsing sequence data (structure = hash of arrays)
	  push @{ $annotations{$scaff}{"gene_$localgenecnt"}{"starts"} }, $array[3];
	  push @{ $annotations{$scaff}{"gene_$localgenecnt"}{"stops"} }, $array[4];
	  }

	} #end of if(fastblock == 0), so at that point, all annotation data should be stored into %annotation


	### FASTA block
	## if we are into the sequence block (i.e. fastablock == 1)
	if($localgenecnt > 0 && $fastablock == 1){
	  $_ = $gff;
	  s/\.gff//;
	  my $outname = $_;
	  $outname = basename($outname);
	  
	  print "Print output in $outfolder/cds.$outname.fa\n";
	  open(OUT, ">$outfolder/cds.$outname.fa");

	  $skip++;

	  if($skip == 2){
	    # switch to slurp mode
	    local $/ = undef; #switch to slurp, mode
	    my $seq = <IN>;
	    chomp($seq);
	    $_ = $seq;
	    s/\r|\n//g;
	    $seq = $_;

	    # split seq into nucleotides
	    my @nucl = split(undef, $seq);
	    my $test = join("", @nucl[0..59]);
	    print "First 60 nucleotides: $test\n";

	    # loop over stored features with annotations hash
	    foreach $gene (keys( %{$annotations{$scaff}} )){
	      
	      print(OUT ">$scaff\_$gene\n");

	      # retrieve start-stop positions and print sequence to output
	      my @starts = @{ $annotations{$scaff}{$gene}{"starts"} };
	      my @stops = @{ $annotations{$scaff}{$gene}{"stops"} };
	      my $nfeat = @starts;

# 	      print "Extracting $nfeat $feature from $gene from $scaff\n";	      
	      for $i (0..$nfeat-1){
		my $sta = $starts[$i];
		my $sto = $stops[$i];
		my @feat = @nucl[$sta .. $sto];
		my $prnt = join("", @feat);
# 		print("$sta\t$sto\t$prnt\n");
		print(OUT "$prnt");
		}
	      print(OUT "\n");
	      }
	    }
	  
	  }
	} #end of while(<IN>)
      print "found $localgenecnt genes and $localfeatcnt $feature\n#####\n\n";

    } #end of if($gff =~ //)
  close(OUT);    
  } #end of gff file loop
print "##########\nVisited $filecnt GFF files, found $totgenecnt genes in total.\n##########\n";
