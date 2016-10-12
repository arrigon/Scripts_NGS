#!/bin/perl

my $readfolder = $ARGV[0];
$readfolder =~ s/\/$//;

# my $outfile = $ARGV[1];
my $outfile = "DMPS.miraOK.conf";


##### Params section
my $generalparams = "
project = MappingPyrgusDMPS
job = genome,mapping,accurate
parameters = -NW:cmrnl=no \\
SOLEXA_SETTINGS -CO:msr=no
";

my $refparams = "
readgroup = Reference
is_reference
data = ref/Erynnis.fna
strain = Pyrgus
";

# index all reads / specimens
my @allfiles = `ls -1 readsDMPS/`;
my %index;
foreach my $file (@allfiles){
  chomp($file);
  my $specimen = $file;
  $specimen =~ s/\.fq\_(.$)//;
  my $readnr = $1;
  $index{$specimen}{$readnr} = $file;
  }


open IN, "bin/SampleIDs";
my %name;
while(<IN>){
  chomp();
  my @fields = split("\t", $_);
  $name{$fields[0]} = $fields[1]
  }


#### Start file production
open OUT, ">$outfile";

# print general params
print OUT "$generalparams\n$refparams\n\n";

# print specimen params
foreach my $spec (keys(%index)){
  my $realname;
  if($name{$spec}){
    $realname = $name{$spec};
    } else {
    $realname = $spec;
    }

  my $file1 = $index{$spec}{1};
  my $file2 = $index{$spec}{2};

  my $specparam = "readgroup = $spec
  data = fastq::$readfolder/$file1 fastq::$readfolder/$file2
  technology = solexa
  strain = $realname
  template_size = 50 300
  segment_placement = ---> <---";

  print OUT "$specparam\n\n";
  }

print "done";
