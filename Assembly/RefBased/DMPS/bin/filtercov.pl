my $infile = "1_2_1.fq_1.fas.clstr";
my $outfile = "test";
my $mincov = 2;
my $maxcov = 50;

filtercov($infile, $outfile, $mincov, $maxcov);


sub filtercov {
  my ($infile, $outfile, $mincov, $maxcov) = ($_[0], $_[1], $_[2], $_[3]);

  # open files
  open(IN, "$infile");
  
  # initiate values
  my $clstr;
  my %store;
  my $cnt;
  
  
  # load file in RAM
  while(<IN>){
    chomp();
    $line = $_;
    if($line =~ />Cluster/){
      $clstr = $line;

      } else {
	$line =~ /.*\t(\d+)nt, >(.*)\.\.\..*/;
	my $len = $1;
	my $rd = $2;
      
	$store{$clstr}{$rd} = $len;
      }          
   }
   close(IN);
   
  # now visit each cluster, and save to @keepers
  foreach $clstr (keys(%store)){
    my @tmp = keys(%{$store{$clstr}});
    if($#tmp > $mincov){
      push(@keepers, $clstr);
      }
    }
   
  
  # save to output
  open(OUT, ">$outfile.list");
  open(LOG, ">$outfile.log");

  foreach $clstr (@keepers){
    my %tmp = %{$store{$clstr}};
    
    # list included reads guys in descending length order
    $cnt = 1;
    foreach (sort { ($tmp{$b} <=> $tmp{$a}) || ($b <=> $a) } keys %tmp){
      my $rd = $_;
      if($cnt <= $maxcov){
	# print "$clstr\t$cnt\t$tmp{$rd}\n";
	print(OUT "$rd\n");
	print(LOG "$clstr\t$rd\n");
	$cnt++;
	} else {
	next();
	}
      }
    }
   
   print "done\n";
   close(OUT);
   close(LOG);
  }
   