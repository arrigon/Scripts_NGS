### How to best use this pipeline suite:
1. Copy the complete pipeline to your home folder; using a link -> This will ensure that you always use the most up-to-date version.
cp -l /archive/Scripts_NGS/Pipeline_CleanDemultiplex/ /home/YourUserName/YourProjectName

2. Run it from your home.


### Usage
The cleaning pipeline is designed to work with:
- The naming schemes of reads are for the moment following those made by TGAC and the CIG
- The naming of files is based on the CIG practices 
(if needing to rename files: use PrepareFileNames.pl, in 
/archive/Scripts_NGS/Preparing_reads/ReadsHandling/FilenameHandling
- Paired ends; if you have only single ends: you have to fake some inputs, using script FakePEdata.pl, in 
/archive/Scripts_NGS/Preparing_reads/ReadsHandling/FilenameHandling


then make sure that $namestyle = "checked"
(e.g. XXX_YYY_R1_001.fastq.gz, we use a regexp to tell apart XXX, R1 and 001 (YYY is ignored) -> your file names MUST finish with R1_number.fastq.gz)

Note that you can use a file renaming script to enforce these names. This is the most annoying part.


USAGE:
perl CleanFastq_Pipeline_queue3.pl PathToRawReads param/PathToBarcodes NameOfOutfolder MaxCPUUsage


NOTE: you MUST edit a bunch of parameters in CleanFastq_Pipeline_queue3.pl

my $bcmismatch = 2; #Maximum number of allowed mismatches in barcode
my $keepbc = 0; #keep ($keepbc = 1) or discard ($keepbc = 0) barcode
my $motif = "keepall"; # For stacks only: keep only reads starting with restriction enzyme. Sbf1, EcoR1 are available so far. If needing additional motifs, add them to bin/RestrictionSiteFilterFastqIII.pl
		    # N.B. only using $motif = keepall works in combination with barcode clipping options 
		    # i.e. $motif = "keepall" and $keepbc = 1 will leave the barcode in the reads
		    #      $motif = "keepall" and $keepbc = 0 will ALWAYS trim out the barcode
		    #      $motif = EcoR1 or others will always trim out the barcode
my $MaxLen = -1; # For stacks only: trim all reads to that length (usually 80), discard shorter reads. Put MaxLen < 0 to bypass that treatment
my $MinLen = 25; #discard reads shorter than that value (fastq-mcf step)
my $phred = 28;  #trim end of reads if phred scores goes below that value
my $namestyle = "checked"; # make sure we use the right RegExp to index raw fastq files (rename your input fastq with CheckFileNames.pl, in /data/data). 
			  # Best option so far: rename your raw fastq files BEFORE starting, in order to comply with "checked" format
			  # To do so: use script bin/CheckFileNames.pl and read help in there to have it at work.
			  # format = checked (as produced by bin/CheckFileNames.pl): fastq files are as XXX.R1.001.fastq.gz, where XXX is laneID, R1 = readID (R1/R2), 001 = specimenID / fastqchunckID
			  # format = CIG : fastq files are named as XXX_YYY_R1_001.fastq.gz, where XXX is laneID, YYY is skipped, R1 = readID (R1/R2), 001 = specimen / chunckID
my $trimR2 = 0;	# Trim ($trimR2 = 1) second read to $MaxLen (used in double-digest + paired-end RADs) or leave as is ($trimR2 = 0)
		# Note that is inactive if MaxLen < -1.
my $stitch = 0; # Append ($stitch = 1) R1 (orientation standard) and R2 (revcomp) altogether in a single fastq file1  (used in double-digest + paired-end RADs)
		# if $stitch = 0, store R1 (orient. std) and R2 (orient. std) separately



Also, look at help and comments in CleanFastq_Pipeline_queue.pl
Good luck



#### Running the script on several libraries: look at BatchDemultiplex_Rad.sh and adapt it according to your needs.