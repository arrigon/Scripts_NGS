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
data = ref/Melipona.fna
strain = contig1
";

my $specparams = "readgroup = $specimen
data = fastq::$readfolder/Meg$specimen.fq_1
technology = solexa
strain = Mito$specimen";


#### Start file production
open OUT, ">$outfile";
print OUT "$generalparams\n$refparams\n\n$specparams";
close(OUT);


print "done";
