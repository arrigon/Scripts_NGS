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

inputs: keep file name as simple as possible. Avoid long names with spaces.