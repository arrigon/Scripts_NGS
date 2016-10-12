#!/bin/perl

my $infile = "allseq.fas";


## import fasta file
open(FILE, "$infile");
local $/ = undef; #slurp, mode
my $input = <FILE>;
local $/ = "\n";
my @fields = split(/\>/, $input);
shift(@fields);
close(FILE);

## load it into hash. Note that %fasta contains your sequences. Can be reused elsewhere.
my %fasta;
my @accs;
foreach $input (@fields){
  my @tmp = split(/\n/, $input, 2);
  $_ = $tmp[1];
  s/\r|\n//g;
  $acc = $tmp[0];
  if($acc =~ /(.*)\_ref/){
    $seq = uc $_; #turn everything in uppercase
    $fasta{$1}{"ref"} = $seq; 
    push(@accs, $1);
    } else {
    $seq = uc $_; #turn everything in uppercase
    $fasta{$acc}{"seq"} = $seq; 
    push(@accs, $acc);
    }
  }


### Prepare alnmt hash
my %alnmt;
$alnmt{"ref"}{"ref"} = $fasta{$accs[0]}{"ref"};
$alnmt{"ref"}{"seq"} = $fasta{$accs[0]}{"seq"};


### loop over all sequences of fasta, add progressively each seq to the existaing aligmnt
foreach $acc (keys(%fasta)){

  # load new sequence into alnmt
  $alnmt{$acc}{"ref"} = $fasta{$acc}{"ref"};
  $alnmt{$acc}{"seq"} = $fasta{$acc}{"seq"};

  # check if it needs to be aligned, or if it modifies the existing alnmt
  COMP: #return point, come back here each time you corrected the alnmt for the focal sequence

  # get reference seq and new seq
  my $ref = $alnmt{"ref"}{"ref"};
  my $cand = $alnmt{$acc}{"ref"};
  my $candseq = $alnmt{$acc}{"seq"};

  my @ref = split("", $ref);
  my @cand = split("", $cand);
  my @candseq = split("", $candseq);
  
  # compare them, position by position. WARNING: we assume that only gaps are differing
  foreach $pos (0..$#ref){
    
    # get pos
    my $candpos = $cand[$pos];
    my $refpos = $ref[$pos];

    # compare it
    if($candpos ne $refpos){ # if mismatches, check where is the gap

      if($candpos eq "-"){ #gap is in new sequence, must add it to ref and propagate into complete alignment
	foreach $prop (keys(%alnmt)){ # propagate gap into alnmt existing so far (including refseq)
	  if($prop ne $acc){
	    print "Gap on candseq $acc at pos $pos, updating $prop\n";

	    my @tmp = split("", $alnmt{$prop}{"ref"});
	    my $repl = $tmp[$pos-1];
	    $tmp[$pos-1] = "$repl-";
	    $alnmt{$prop}{"ref"} = join("", @tmp);

	    my @tmp = split("", $alnmt{$prop}{"seq"});
	    my $repl = $tmp[$pos-1];
	    $tmp[$pos-1] = "$repl-";
	    $alnmt{$prop}{"seq"} = join("", @tmp);
	    }
	  }

	goto COMP; #keep checking pos after pos
	}

      if($refpos eq "-"){ #gap is in reference, must only update new sequence
	print "Gap on refseq on pos $pos, updating candseq\n";
	
	# propagate the gap into alnmt
	my $repl = $cand[$pos-1];
	$cand[$pos-1] = "$repl-";
	$alnmt{$acc}{"ref"} = join("", @cand);

	$repl = $candseq[$pos-1];
	$candseq[$pos-1] = "$repl-";
	$alnmt{$acc}{"seq"} = join("", @candseq);

	goto COMP; #keep checking pos after pos
	}
      }
    }
  }


open(OUT, ">ALN.fas");
foreach $acc (keys(%alnmt)){
  my $seq = $alnmt{$acc}{"ref"};
  print OUT ">$acc\n$seq\n";
  my $seq = $alnmt{$acc}{"seq"};
  print OUT ">$acc\n$seq\n";
  }
close(OUT);