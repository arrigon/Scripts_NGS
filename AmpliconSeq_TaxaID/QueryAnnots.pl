#!/bin/perl
#####################################
#### Take as input a set of contigs, blasts them against GenBank and extracts sequences of interest
#### using a keyword against taxonomic identifications
####
#### Usage: perl QuerryAnnot.pl ContigFile KeyWord AnalysisMode
####        - ContigFile = fasta file containing the contigs
####        - KeyWord = your taxonomic group of interest (from species to order)
####        - Analysis mode = either CompleteSeq (the complete contig is retrieved) 
####                          or HitsOnly (only regions producing blast hits are extracted)
####
#### The outputs are produced in data.out/
#### There are some hard-coded params about blast sensitivities, look at lines 28-32.
####
#### Dependencies: blastn, blastdbcmd, File::Basename
####
#### Nils Arrigo, Unil 2015
#####################################


## usage: perl QueryAnnot.pl data.in/Contigs.fasta Pinales keepall|hits
my $file = $ARGV[0]; #contig file, must be in data.in/, formatted as from Geneious assemblies.
my $keyword = $ARGV[1]; #keyword to look after, and extract contigs that match with it.
my $mode = $ARGV[2]; #either CompleteSeq = extract complete sequence or HitsOnly = keep only region generating a blast hit
chomp($mode);


## Important params (hard coded here)
# blast searches
my $nthreads = 4; #nthreads for blast search
my $eval = 0.1; #min e-value for blast search
my $pctsim = 30; #min % simil for blast search


#### checks
if($mode ne "CompleteSeq" && $mode ne "HitsOnly"){
  print "STOP\: Incorrect analysis mode, make sure to use either CompleteSeq or HitsOnly mode!\n";
  exit;
  }
  
print "\n###### QueryAnnot Pipeline ######\nContigs file: $file\nBlast params:\n-threads: $nthreads\n-evalue: $eval\n-pctsimil: $pctsim\n\nAnalysis mode: $mode\n\n\n\n";


#### Settings
use File::Basename;
# paths to databases
my $db = "/data/GenBank/db/nt";
my $refpath = "/archive/Scripts_NGS/AmpliconSeq_Decontamination/DecontPipeline/refs/";


### Step 0. Prepare working folders
my $command = "mkdir tmp/ data.out/";
print("$command\n");
system($command); 

my $bsn = basename($file);
$bsn =~ s/\..*$//;


### Step 0. Rename contigs
my $command = "perl bin/NumberHeadsFasta.pl $file tmp/$bsn.fas tmp/$bsn.idx";
print("$command\n");
system($command);


unless( -e "tmp/$bsn.blast"){
  ## Step 1. Run blast search
  my $command = "blastn -query tmp/$bsn.fas -db $db -task blastn -dust yes -num_threads $nthreads -evalue $eval -perc_identity $pctsim -out tmp/$bsn.blast -outfmt 6 -max_target_seqs 1";
  print("$command\n");
  system($command);

  if($mode =~ "CompleteSeq"){
    my $command = "perl bin/DeleteBlastDoublons2.pl tmp/$bsn.blast tmp/$bsn.blast"; #clean blast outputs
    print("$command\n");
    system($command);
    }
    

  ### Step 2. Collect annotations
  my $command = "awk -F \'\\t\' \'{print \$2}\' tmp/$bsn.blast | awk -F \'|\' \'{print \$4}\' > tmp/$bsn.gb";
  print("$command\n");
  system($command); 

  open(IN, "tmp/$bsn.gb");
  my @accs;
  while(<IN>){
    chomp();
    push(@accs, $_);
    }
  my $entries = join(",", @accs);
  close(IN);

  my $command = "blastdbcmd -db $db -entry $entries -outfmt \"%a\t%t\" > tmp/$bsn.annot";
  print("$command\n");
  system($command);   

  my $command = "./bin/blast2lca -names $refpath/names.dmp -nodes $refpath/nodes.dmp -dict $refpath/gi_taxid.bin -levels=superkingdom:kingdom:phylum:order:family -order=true tmp/$bsn.blast --savemem > tmp/$bsn.taxo";
  print("$command\n");
  system($command); 
  }
    
    
## Step 3. Collect contigs matching for a keyword of interest, extract hits or keep complete sequences
my $command = "perl bin/FetchFasta\_KeyWord.pl tmp/$bsn.fas tmp/$bsn.taxo tmp/$bsn.blast $keyword $mode data.out/$bsn.$keyword.$mode.clean";
print("$command\n");
system($command); 


### Step 4. Put in perspective of annotations
my $command = "perl bin/CrossParseAnnots.pl tmp/$bsn.blast tmp/$bsn.taxo tmp/$bsn.annot tmp/$bsn.idx data.out/$bsn.$keyword.$mode.clean data.out/$bsn.$keyword.$mode.infos";
print("$command\n");
system($command); 


### Step 5. Clean mess
# my $command = "rm -rf tmp/";
# print("$command\n");
# system($command); 

