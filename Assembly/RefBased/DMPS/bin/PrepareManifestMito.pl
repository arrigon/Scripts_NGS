#!/bin/perl

my $specimen = $ARGV[0];
my $outfile = $ARGV[1];

$readfolder = "reads";

chomp($specimen);
chomp($outfile);

##### Params section
my $generalparams = "
project = Assembly$specimen
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

my $specparams = "readgroup = $specimen
data = fastq::$readfolder/131029SN365AL00*HBQ_$specimen.fq_1 fastq::$readfolder/131029SN365AL00*HBQ_$specimen.fq_2
technology = solexa
strain = Mito$specimen
template_size = 50 300
segment_placement = ---> <---";


#### Start file production
open OUT, ">$outfile";
print OUT "$generalparams\n$refparams\n\n$specparams";
close(OUT);


print "done";
