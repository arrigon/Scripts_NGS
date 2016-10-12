rm(list = ls())


# Get arguments from command line (should be uncommented in production version)
# args=(commandArgs(TRUE))
# for(i in 1:length(args)){
#   eval(parse(text=args[[i]]))
#   }



# DEBUG params (should be commented in production version)
ncores = 10
tsv = "../data.in/STI1.SNP"
cov = "../data.in/STI1.COV"
workfolder = "../savedanalyses"
outfolder = "../data.out"
bestscorers = 0 #parameter to pick-up best specimens. these are defined as having more than bestscorers * mean(NrStacks / specimen)
minTaxDepth = 0.1 #parameter to pick-up loci, minimum proportion of taxa in which locus is observed
maxAlleles = 3
maxlenSNP = 5
keepmono = 1 #0 = remove monomorphic stacks or 1 leave them in there


# for log purposes, display input parameters. Leave as it is.
ncores
tsv
cov
bestscorers
minTaxDepth
maxAlleles


# prepare outfile names
bsn = gsub(".tsv", "", basename(tsv))
outfile = paste(bsn, "_bscor", bestscorers, "_taxdepth", minTaxDepth, "_maxallele", maxAlleles, "_maxLen", maxlenSNP, "_mono", keepmono, sep = "")



######## GENERAL TOOLS
### Get number of alleles per specimen
nall_cell = function(val){
  nall = NA
  if(is.na(val) == F){
    if(val != ""){
      nall = length(as.numeric(strsplit(val, "/")[[1]]))
      }
    } else {
    nall = NA
    }
  nall 
  }
nall_row = function(vec) sapply(vec, nall_cell)

### Get average coverage per rad and per specimen
mcov_cell = function(val){
  nall = NA
  if(is.na(val) == F){
    if(val != ""){
      nall = mean(as.numeric(strsplit(val, "/")[[1]]))
      }
    } else {
    nall = NA
    }
  nall 
  }
mcov_row = function(vec) sapply(vec, mcov_cell)


######## ANALYSIS
# First, we check whether these data have not been already processed
path_to_workfiles = dir(workfolder, pattern = ".Rdata")
path_to_workfiles = path_to_workfiles[grep(basename(tsv), path_to_workfiles)] #make sure we pick the right file in there

if(length(path_to_workfiles) == 0){ #no file saved from a previous session is available, we have to compute everything.
  ## Get data
  # initiate multithreading
  library(doMC)
  library(foreach)
  library(multicore)
  registerDoMC(cores=ncores)

  # open files
  data = read.delim(tsv, header = T, row.names = 1)
  data.cov = read.delim(cov, header = T, row.names = 1)
  cls = colnames(data)
  colnames(data) = cls

  # get generic filenames
  tsv = basename(tsv)
  cov = basename(cov)

  # get out occurrences
  occs = data[, 12:ncol(data)]
  occs.cov = data.cov[, 12:ncol(data.cov)]
  occs = as.matrix(occs)
  occs.cov = as.matrix(occs.cov)

  # Compute number of alleles per sample
  # NALL = t(apply(occs, 1, nall_row)) #usual, single core method
  NALL = foreach(i=1:nrow(occs), .combine=rbind)%dopar%{ #multithreaded
    nall_row(occs[i, ])
    }
  rownames(NALL) = rownames(occs)

  # compute average coverage per rad / sample (average over alleles)
  MCOV = foreach(i=1:nrow(occs.cov), .combine=rbind)%dopar%{ #multithreaded
    mcov_row(occs.cov[i, ])
    }
  rownames(MCOV) = rownames(occs.cov)

  # make sure we save the produced file in a binary format. will save time for future tests
  save(list = c("data", "data.cov", "tsv", "cov", "NALL", "MCOV"), file = paste(workfolder, "/WorkingFiles_", tsv, ".Rdata", sep = ""))

  } else { #Bingo, we have a file already in here
  ## Get data
  load(paste(workfolder, path_to_workfiles, sep = "/"))

  # get out occurrences
  occs = data[, 12:ncol(data)]
  occs.cov = data.cov[, 12:ncol(data.cov)]
  occs = as.matrix(occs)
  occs.cov = as.matrix(occs.cov)

  }

##### Collect informations about the dataset
# Count how many stacks are present in each specimen
OUT = NALL
OUT[OUT > 1] = 1
OUT[is.na(OUT)] = 0
LocDepth = colSums(OUT)

# Keep only best scoring specimens (selection made relative to mean)
passedspecimens = LocDepth >= (bestscorers * mean(LocDepth))
OUT = OUT[, passedspecimens]
NALL = NALL[, passedspecimens]
MCOV = MCOV[, passedspecimens]


# Detect loci with more than 2 alleles / specimen (filtering step applied later)
TaxDepth = rowMeans(OUT)

# Detect loci with more than 2 alleles / specimen (filtering step applied later)
nbAlleles = apply(NALL, 1, max, na.rm = T)

# Get how many SNPs segregate per locus
lenSNP = data$Num.SNPs

# Check whether a given locus is polymorphic
if(keepmono == 1){
  mono = rep(0, nrow(occs))
  } else {
  mono = occs == "consensus"
  mono = rowSums(mono, na.rm = T)
  }


##### Clean dataset
# identify STACKs of interest
target = which(TaxDepth > minTaxDepth & nbAlleles < maxAlleles & lenSNP < maxlenSNP & mono == 0)

image(t(log10(MCOV[target,])))
missing = 1 - mean(OUT[target,])
missing

heatmap(t(1-OUT[ round(seq(1, nrow(OUT), length.out = 1000)),]), Colv = NA, scale = "none", labCol = NULL)
heatmap(t(1-OUT[ target,]), Colv = NA, scale = "none", labCol = NULL)

info = read.delim("../info/STI1.txt")
info = info[match(colnames(OUT), info$Tag),]
cols = as.factor(as.character(info[, 2]))
levels(cols) = rainbow(nlevels(cols))
cols = as.character(cols)

pdf(file = "MapSTI.pdf")
heatmap(t(1-OUT[ target,]), Colv = NA, scale = "none", labCol = NULL, RowSideColors = cols, cexRow = .4)
dev.off()



image(t(log10(MCOV[target, colSums(OUT[target,]) > 200])))
mean(is.na(MCOV[target, colSums(OUT[target,]) > 200]))

length(target)

boxplot(apply(MCOV[target,], 1, mean, na.rm=T), 
	apply(MCOV[target,], 2, mean, na.rm=T), 
	ylab = "SD Coverage",
	names = c("SD Intra Specimen", "SD Intra Locus"))



# filter them out of dataset
final = data.frame(CatalogID = rownames(data[target, ]), Consensus=data[target, ]$Consensus.Sequence, SNPs = data[target, ]$SNPs, occs[ target, passedspecimens])
final.cov = data.frame(CatalogID = rownames(data[target, ]), Consensus=data[target, ]$Consensus.Sequence, SNPs = data[target, ]$SNPs, occs.cov[ target, passedspecimens])
final[final == ""] = -9
final.cov[final.cov == ""] = -9

##### Save out results
write.table(final, file = paste(outfolder, "/", outfile, ".SNPs", sep = ""), quote = F, sep = "\t", row.names = F)
write.table(final.cov, file = paste(outfolder, "/", outfile, ".COVs", sep = ""), quote = F, sep = "\t", row.names = F)
