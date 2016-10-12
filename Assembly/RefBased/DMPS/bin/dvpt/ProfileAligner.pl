#!/bin/perl
### Works only of reference sequences are strictly identical. not always the case...


# fake data
my $reference = "ATTTCG-----CGCGTT-";

my $seq1 = "ATTTCGCTT-GCG-TT";
my $seq1b = "--TTCGCTT-------";
my $seq2 = "AT-TTCGCGCG-T-T";
my $seq2b = "---------CG-T-T";
my $seq3 = "AT-TTC--GCGCG-T-T";
my $seq3b = "-----------CG-T-T";

my %fasta;
$fasta{"seq1"}{"ref"} = $seq1;
$fasta{"seq1"}{"seq"} = $seq1b;
$fasta{"seq2"}{"ref"} = $seq2;
$fasta{"seq2"}{"seq"} = $seq2b;
$fasta{"seq3"}{"ref"} = $seq3;
$fasta{"seq3"}{"seq"} = $seq3b;




### Prepare alnmt hash
my %alnmt;
$alnmt{"ref"}{"ref"} = $reference;
$alnmt{"ref"}{"seq"} = $reference;


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


foreach $acc (keys(%alnmt)){
  my $seq = $alnmt{$acc}{"ref"};
  print "$acc\t$seq\n";
  my $seq = $alnmt{$acc}{"seq"};
  print "$acc\t$seq\n";
  }
