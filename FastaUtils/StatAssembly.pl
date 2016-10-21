#!/usr/bin/perl -w

######################################################################################
# This script will translate a GenBank (*.seq) file into a fasta file        #
#                                                                                    #
# Dependencies: bioperl								     #
#                                                                                    #
# Usage: perl GB2FAS.pl						                     #
#                                                                                    #
# This script can be run in the folder where several GenBank files are stored        #
#                                                                                    #
# copyright: Nils Arrigo	                                                     #
#            EEB, Barker Lab			                                     #
#            University of Arizona                                                   #
#                                                                                    #
# You may distribute this script as long as this copyright notice is not removed     #
#                                                                                    #
######################################################################################

use Bio::SeqIO;
use IO::String;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

# get list of GenBank files to convert
@dir=`ls *.fasta *.fas *.fa *.scafSeq *.contig`;

open(OUT, ">>AssemblyStats.txt") or die;
print OUT "Dataset\tNscaff\tTotLen\tMaxLen\tN50\n";

foreach $gbfile (@dir){
  # prepare outfile
  chomp($gbfile);
  print "Analyze $gbfile\n";

  # open GenBank file and proceed to conversion
  $seqin=Bio::SeqIO->new('-file'=>$gbfile, '-format'=>'fasta');

  my @stc = ();
  while($seq=$seqin->next_seq) {
    $len = $seq->length;
    push(@stc, $len);
    }
    
    # Get Nscaff
    $nscaf = $#stc;
  
    # Get max scaflength
    my $maxl = max(@stc);
    
    # Get total length
    my $totlen = sum(@stc);

    ## Get N50
    @stc = sort par_num @stc;
    my @cumsum = (); 
    push(@cumsum, $stc[0]);
    
    my $i = 0;
    my $test = 0;
    until($test > $totlen/2){
      $N50 = $stc[$i];
      $test = $test + $N50;
      $i++;
      }
    
    print OUT "$gbfile\t$nscaf\t$totlen\t$maxl\t$N50\n";
  }
close OUT;
exit(0);


sub par_num { return $a <=> $b };
