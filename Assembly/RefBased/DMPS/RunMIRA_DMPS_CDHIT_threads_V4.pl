#!/bin/perl

##### Assembly of DMPS datasets (or any other ref-based stuff), using MIRA, Ref-guided assembly from PAIRED-END reads
### Performed steps:
# 1. Remove redundant reads, to make assembly more efficient and faster (reads clustering using cdhit, with defaults)
# 2. Assemble reads against reference.
# 3. Rename assembled sequences with final specimen names (if needed)
# 4. Collect and save results in data.out/
#
### Mutlithreading: by default, each job uses 5 cores (cdHIT) and 
#	            the script launches jobs untill 20cores OR 0Gb of free RAM 
#                   (this is safe, this script underestimates available ram)
#
### Usage and parameters
# perl RunMIRA_DMPS_CDHIT_threads.pl readsfolder reference specID
# - readsfolder = repository of fastq files. 
#		  Two fastq per specimen, using paired-end data
#
# - Reference: should be in ref/XXX.fna = reference sequence, in fasta format 
#	       (name with *.fna to make MIRA accept it as is) 
#
# - specID = should be in params/SpecimenID.txt 
#                   Used to rename specimens before assembly, optional
#
###
# Working scripts / external programs: must be in bin/ folder
#
### Nils Arrigo, Pyrgus project 2013.
########################################################

#### libraries
use threads;
use threads::shared;
use File::Basename;
use Sys::MemInfo qw(totalmem freeswap freemem);
use POSIX;


#### Script arguments
my $readfolder = $ARGV[0];
my $refseq = $ARGV[1];
my $specID = $ARGV[2];
$readfolder =~ s/\/$//;
my $useCDHIT = 1; #YES, use cdhit! (put 0 otherwise, but stick to cdhit if possible)
my $mincov = 3; #for cdhit step 
my $maxcov = 25; #for cdhit step

my $scriptname = "MIRA\_DMPS";


# prepare working folders
system("rm -rf tmp data.out");
system("mkdir -p tmp tmp/cdhit tmp/cdhit/fas tmp/cdhit/clstr tmp/cdhit/lists tmp/fq tmp/mira tmp/mira/wd tmp/mira/out tmp/sai tmp/sam tmp/bam data.out");


# index all reads per specimens
my @allfiles = `ls -1 $readfolder`;
my %index;
foreach my $file (@allfiles){
  chomp($file);
  my $specimen = $file;
  $specimen =~ s/\.fq\_(.$)//;
  my $readnr = $1;
  $index{$specimen}{$readnr} = $file;
  }


# load real specimen names
open IN, "$specID";
my %name;
while(<IN>){
  chomp();
  my @fields = split("\t", $_);
  $name{$fields[0]} = $fields[1]
  }



### Start pipeline
my @threads;

my $startedthreads = 0;
foreach my $spec (keys(%index)){

### for debug, running only one specimen
# # # my @dddd = keys(%index);
# # # @dddd = @dddd[0];
# # # foreach my $spec (@dddd){

  my $realname;
  if($name{$spec}){
    $realname = $name{$spec};
    } else {
    $realname = $spec;
    }

  my $file1 = $index{$spec}{1};
  my $file2 = $index{$spec}{2};

  # inspect ram before firing up
  my $freeRAM = ceil((&freemem) / 1073741824);
  
  # inspect CPU
  open PIPE, "uptime |";
  my $line = <PIPE>;
  close PIPE;
  $line =~ s/\s//g;
  my @lineArr =  split /:/, $line;
  my $times = $lineArr[@lineArr-1];
  my @timeArr = split /,/, $times;
  my $loadCPU = $timeArr[0] + $CPU;


  print "CPU load = $loadCPU\tFree memory = $freeRAM Gb\n";

  if(($loadCPU + 2) < 30 and $freeRAM > 0) {
    #### Automated loader
    my $t = threads->new(\&runMIRA, $realname, $readfolder, $file1, $file2, $refseq);
    push(@threads, $t);
    $startedthreads++;
    sleep(25);

    } else {
    sleep(25);
    redo;
    }
  }


my $cnt = 0;
foreach (@threads) {
  my $genome = $_-> join;
  $cnt++;
  print "Finished job $cnt / $startedthreads\n";
  }


print "Collecting assemblies...\n";
my $command = "perl bin/CleanCollectFasta.pl tmp/mira/out $mincov data.out/DMPSassembly.fasta";
print("$command\n");
system($command);


my $command = "cp -r tmp/mira/out data.out/assemblystats";
print("$command\n");
system($command);


my $command = "rm -rf tmp/";
print("$command\n");
system($command);


print "Done.\n";



### Routines
sub runMIRA {
  my ($realname, $readfolder, $file1, $file2, $refseq) = ($_[0], $_[1], $_[2], $_[3], $_[4]);

  
  if($useCDHIT eq 1){

    ### Run CDhit to reduce redundancy and to improve recovery of rare reads (i.e. PCR amplicons that were less successful)
    # Avoids overkilling MIRA and improves assembly
    # step 1: convert to fasta
    my $command = "perl bin/Fastq2Fasta.pl $readfolder/$file1 tmp/cdhit/fas/$file1.fasta";
    print("$command\n");
    system($command);

    my $command = "perl bin/Fastq2Fasta.pl $readfolder/$file2 tmp/cdhit/fas/$file2.fasta";
    print("$command\n");
    system($command);



    # step 2: run CDHit
    my $command = "./bin/cd-hit-454 -i tmp/cdhit/fas/$file1.fasta -o tmp/cdhit/clstr/$file1.fas -T 5 -M 20000 -c 0.9";
    print("$command\n");
    system($command);

    my $command = "./bin/cd-hit-454 -i tmp/cdhit/fas/$file2.fasta -o tmp/cdhit/clstr/$file2.fas -T 5 -M 20000 -c 0.9";
    print("$command\n");
    system($command);



    # step 4: shave-off excess of reads in file 1 (keep clusters supported by at least mincov reads, keep no more than maxcov reads)
    filtercov("tmp/cdhit/clstr/$file1.fas.clstr", "tmp/cdhit/lists/$file1", $mincov, $maxcov);
    filtercov("tmp/cdhit/clstr/$file2.fas.clstr", "tmp/cdhit/lists/$file2", $mincov, $maxcov);

    my $command = "cat tmp/cdhit/lists/*.list > tmp/cdhit/lists/$file1.all";
    print("$command\n");
    system($command);



    ### Step 5. Collapse selected reads from R1 and R2 into a single list
    open(IN, "tmp/cdhit/lists/$file1.all");
    my %readindex;
    while(<IN>){
      chomp();
      my $read = $_;
      $read =~ /(.*)/;
      my $ID = $1;
      $readindex{$ID} = 1;
      }
    close(IN);

    open(OUT, ">tmp/cdhit/lists/$file1.final");
    foreach my $read (keys(%readindex)){
      print OUT "$read\n";
      }
    close(OUT);



    # step 6: extract these reads from original fastq
    my $command = "perl bin/FetchFastq.pl $readfolder/$file1 tmp/fq/$file1 tmp/cdhit/lists/$file1.final";
    print("$command\n");
    system($command);

    my $command = "perl bin/FetchFastq.pl $readfolder/$file2 tmp/fq/$file2 tmp/cdhit/lists/$file1.final";
    print("$command\n");
    system($command);


    # step 7: measure deamination in these bastards
    my $command = "perl bin/RunSamtools.pl tmp/fq/$file1 $refseq";
    print("$command\n");
    system($command);

    # step 7: measure deamination in these bastards
    my $command = "perl bin/RunSamtools.pl tmp/fq/$file2 $refseq";
    print("$command\n");
    system($command);



    } else { 

    ### Adjust coverage to 100x (per bp), not worth going more than that. 
    # 
    # because that overkills MIRA, without improving the assembly... 
    # Here, Camille could expect at best 8Kbp out of DMPS, 
    # we cut the files to 50,000 reads here (adjust as function of expected coverage)
    # she's still get... 50'000 pairs * 150bp = 937x
    my $command = "head -n 200000 $readfolder/$file1 > tmp/fq/$file1";
    print("$command\n");
    system($command);

    my $command = "head -n 200000 $readfolder/$file2 > tmp/fq/$file2";
    print("$command\n");
    system($command);
    }



  ### Next step Ref guided assembly with MIRA
  ## Params section, run in draft mode, because refseq is not necessaraly close to new seqs. 
  ## will also use best voting mode and force avoiding IUPAC (fnicpst = yes | no)
  my $generalparams = "
  project = Assembly$realname
  job = genome,mapping,draft
  parameters = -NW:cmrnl=no -NW:cac=no -CO:mr=yes\\
  SOLEXA_SETTINGS -CO:msr=no -CO:fnicpst=yes
  ";

  my $refparams = "
  readgroup = Reference
  as_reference
  data = $refseq
  strain = Pyrgus
  ";

  my $specparam = "readgroup = $realname
  data = fastq::tmp/fq/$file1 fastq::tmp/fq/$file2
  technology = solexa
  strain = $realname
  template_size = 50 300
  segment_placement = ---> <---";
  

  # Prepare manifest
  open OUT, ">tmp/mira/wd/manifest.$realname.conf";
  print OUT "$generalparams\n$refparams\n$specparam\n";
  close OUT;



  # Run MIRA
  my $command = "mira tmp/mira/wd/manifest.$realname.conf > tmp/mira/out/assembly.$realname.log";
  print("$command\n");
  system($command);

  
  # Prepare SAM and BAM
  my $command = "miraconvert -f maf -t ace Assembly$realname\_assembly/Assembly$realname\_d\_results/Assembly$realname\_out.maf tmp/mira/out/assembly.$realname.ace";
  print("$command\n");
  system($command);

  my $command = "miraconvert -f maf -t sam Assembly$realname\_assembly/Assembly$realname\_d\_results/Assembly$realname\_out.maf Assembly$realname\_assembly/Assembly$realname\_d\_results/assembly.$realname.sam";
  print("$command\n");
  system($command);

  my $command = "samtools view -bT Assembly$realname\_assembly/Assembly$realname\_d\_results/Assembly$realname\_out\_Pyrgus.padded.fasta Assembly$realname\_assembly/Assembly$realname\_d\_results/assembly.$realname.sam > tmp/mira/out/assembly.$realname.bam";
  print("$command\n");
  system($command);

  my $command = "samtools depth tmp/mira/out/assembly.$realname.bam > tmp/mira/out/assembly.$realname.cov";
  print("$command\n");
  system($command);


  # move outputs to data.out/
  my $command = "cp Assembly$realname\_assembly/Assembly$realname\_d\_results/Assembly$realname\_out\_$realname.unpadded.fasta tmp/mira/out/assembly.$realname.fasta";
  print("$command\n");
  system($command);

  my $command = "cp Assembly$realname\_assembly/Assembly$realname\_d\_results/Assembly$realname\_out\_$realname.padded.fasta tmp/mira/out/assembly.$realname.aln";
  print("$command\n");
  system($command);

  my $command = "cp Assembly$realname\_assembly/Assembly$realname\_d\_results/Assembly$realname\_out\_Pyrgus.padded.fasta tmp/mira/out/assembly.$realname.refaln";
  print("$command\n");
  system($command);

  # clean mess
  my $command = "rm -rf Assembly$realname\_assembly/";
  print("$command\n");
  system($command);

  # clean mess
#   my $command = "rm tmp/fq/$file1 tmp/fq/$file2 tmp/cdhit/fas/$file1.f* tmp/cdhit/fas/$file2.f* tmp/cdhit/clstr/$file1.f* tmp/cdhit/clstr/$file2.f* tmp/cdhit/lists/$file1.* tmp/cdhit/lists/$file2.*";
#   print("$command\n");
#   system($command);
  }


### routine filtering coverages
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
    if(($#tmp + 1) > $mincov){
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
   
  close(OUT);
  close(LOG);
  }



