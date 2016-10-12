rm(list = ls())


### Get script arguments (infolder and outfile)
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }
  
# infolder = paste(getwd(), infolder, sep = "/")
# outfile = paste(getwd(), outfile, sep = "/")

infolder
outfile
getwd()


### load packages
require(ape)

# either scan folder for fasta files
listfasta = dir(infolder, pattern = ".aln", full.names = T)

listord = gsub(".*chunck","", listfasta)
listord = as.numeric(gsub(".aln","", listord))
listfasta = listfasta[order(listord)]
listfasta 

# # or use user selection
# listfasta = c("FAS/1rpl32.seq.fas", "FAS/1rpl32.gaps.fas", "FAS/1trnC.seq.fas", "FAS/1trnC.gaps.fas")
# 
# # give outfile name
# outfile = "CompleteChloro_inclGaps"


# load datasets
ALNS = list()
allspecimens = NULL
for(fas in listfasta){
 data = read.dna(fas, format = "fasta", as.character = T)
 allspecimens = c(allspecimens, rownames(data))
 ALNS[[fas]] = data
 }


# make list of all specimens
allspecimens = levels(as.factor(allspecimens))

# match all alns against that unified list
ALNS.ordered = list()
check = NULL
for(i in 1:length(ALNS)){
  file = names(ALNS)[i]
  aln1 = ALNS[[i]]
  aln1.b = aln1[match(allspecimens, rownames(aln1)),]  
  ALNS.ordered[[file]] = aln1.b
  check = cbind(check, rownames(aln1.b))
  }

# quick check that things look OK
check


# concatenate all matrices
aln.tot = NULL
partitions = NULL
for(i in 1:length(ALNS.ordered)){
  aln1 = ALNS.ordered[[i]]
  aln.tot = cbind(aln.tot, aln1)
  partitions = c(partitions, ncol(aln1))
  }
rownames(aln.tot) = allspecimens


# recode missing data as N
aln.tot[is.na(aln.tot)] = "N"
aln.tot = toupper(aln.tot)


# save matrix, in fasta
cat(file = outfile, append = F)
for(i in 1:nrow(aln.tot)){
  acc = rownames(aln.tot)[i]
  seq = paste(aln.tot[i, ], collapse = "")
  cat(">", acc, "\n", seq, "\n", file = outfile, append = T, sep = "")
  }

# # save matrix, in phylip
# cat(nrow(aln.tot), ncol(aln.tot), "\n", file = paste(outfile, ".phy", sep = ""), append = F)
# for(i in 1:nrow(aln.tot)){
#   acc = rownames(aln.tot)[i]
#   seq = paste(aln.tot[i, ], collapse = "")
#   cat(acc, "\t", seq, "\n", file = paste(outfile, ".phy", sep = ""), append = T, sep = "")
#   }


