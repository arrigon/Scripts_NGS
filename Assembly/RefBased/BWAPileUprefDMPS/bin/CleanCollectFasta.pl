#!/bin/perl


use File::Basename;

my $infolder = $ARGV[0];
my $mincov = $ARGV[1];
my $outfile = $ARGV[2];


my @files = `ls $infolder/*.aln`;
open(OUT, ">$outfile");

@files = sort {$a cmp $b} @files;

foreach my $infile (@files){
  chomp($infile);

  my $bsn = basename($infile);
  $bsn =~ s/\.fasta$//;

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
  foreach $input (@fields){
    my @tmp = split(/\n/, $input, 2);
    $_ = $tmp[1];
    s/\r|\n//g;
    $seq = uc $_; #turn everything in uppercase
    $fasta{$tmp[0]} = $seq; 
    }


  ## load coverage stats
  my $covfile = $infile;
  $covfile =~ s/\.aln$/\.cov/;

  open(COV, "$covfile");
  my @cov;
  while(<COV>){
    chomp();
    my $line = $_;
    my @tmp = split(/\t/, $line);
    push(@cov, $tmp[2] - 1);
    }


  ## print sequences of interest into output
  foreach $head (keys(%fasta)){
    chomp($head);
    my $seq = $fasta{$head};
    $seq =~ s/[x|X]/N/g;
    my @seq = split("", $seq);

    # clean correct sequences
    my $cnt = 0;
    my @clean;
    foreach $bp (@seq){
      if($cov[$cnt] <= $mincov){
	$bp = "N";
	push(@clean, $bp);
	} else {
	push(@clean, $bp);
	}
      $cnt++;
      }

    $seq = join("", @clean);
    $seq =~ s/\*//g;

    print(OUT ">$bsn\n$seq\n");
    }
  }

close(OUT);
