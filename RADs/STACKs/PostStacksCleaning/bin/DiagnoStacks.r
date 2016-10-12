rm(list = ls())


require(hexbin)

# Get arguments from command line (should be uncommented in production version)
# args=(commandArgs(TRUE))
# for(i in 1:length(args)){
#   eval(parse(text=args[[i]]))
#   }



# DEBUG params (should be commented in production version)
ncores = 3
tsv = "../data.in/Tom.SNP"
cov = "../data.in/Tom.COV"
info = "../info/tom.infos"
workfolder = "../savedanalyses"
outfolder = "../data.out"
bestscorers = 0 #parameter to pick-up best specimens. these are defined as having more than bestscorers * mean(NrStacks / specimen)
minTaxDepth = 0
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



######## GENERAL TOOLS
### diagplot function
diagplot = function(x, y, nspl = 30, aggfun = "median", smpl = 1000, ...){
  x.fact = cut(x, nspl)
  x.mw = tapply(x, x.fact, aggfun, na.rm = T)
  y.mw = tapply(y, x.fact, aggfun, na.rm = T)

  smp = sample(1:length(x), smpl, replace = F)
  plot(hexbin(x, y))
  lines(x.mw, y.mw, col = "black")
  }


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


### Compare replicates
comprep = function(i, mat){
  abs(mat[, i] - mat[, i + 1])
  }

### merging function
mergingC = function(A, B, linkA, linkB){
  indsA = linkA
  indsB = linkB
  common = intersect(indsA, indsB)
  matA = A[ match(common, indsA),]
  matB = B[ match(common, indsB),]
  cbind(matA, matB)
  }

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
passedspecimens = which(LocDepth > bestscorers * mean(LocDepth))
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
mono = occs == "consensus"
mono = rowSums(mono, na.rm = T)

# Check heterozyg status
HET = NALL
HET[HET < 2] = 0
HET[HET >= 2] = 1

# Get average coverage
NALL.avg = rowMeans(NALL, na.rm = T)
MCOV.avg = rowMeans(MCOV, na.rm = T)
HET.avg = rowMeans(HET, na.rm = T)

### Get idea of dropout rates, using replicates
# find replicates in there
repl = grep("[_r$|r$]", colnames(NALL))
repl = colnames(NALL)[repl]
orig = gsub("[_r$|r$]", "", repl)
comp = data.frame(couple = rep(1:length(repl), each = 2), Tag = sort(c(orig, repl)))
comp = data.frame(comp, usable = match(comp$Tag, colnames(NALL)))
comp.ok = NULL
for(i in 1:max(comp$couple)){
  comp.sub = comp[ comp$couple == i,]
  if(any(is.na(comp.sub$usable)) == F) comp.ok = rbind(comp.ok, comp.sub)
  }

# get data of replicates
NALL.rep = NALL[, comp.ok$usable]
MCOV.rep = MCOV[, comp.ok$usable]
HET.rep = HET[, comp.ok$usable]


# check allele dropout (via HET)
deltaNALL = sapply(seq(1, ncol(NALL.rep), by = 2), comprep, NALL.rep)
deltaNALL.avg = rowMeans(deltaNALL, na.rm = T)
plot(hexbin(log10(MCOV.avg), deltaNALL.avg, xbins = 200), style = "colorscale")
diagplot(log10(MCOV.avg), deltaNALL.avg, median, nspl = 20, smpl = length(MCOV.avg))

########
diagplot(log10(MCOV.avg), deltaNALL.avg, median, nspl = 20, smpl = length(MCOV.avg))
points(log10(MCOV.avg[target]), deltaNALL.avg[target], cex = .5, pch = 16, col = "red")

diagplot(log10(MCOV.avg), deltaHET.avg, median, nspl = 20, smpl = length(MCOV.avg))
points(log10(MCOV.avg[target]), deltaNALL.avg[target], cex = .5, pch = 16, col = "red")



deltaHET = sapply(seq(1, ncol(HET.rep), by = 2), comprep, HET.rep)
deltaHET.avg = rowMeans(deltaHET, na.rm = T)
plot(hexbin(log10(MCOV.avg), deltaHET.avg, xbins = 200), style = "colorscale")


deltaMCOV = sapply(seq(1, ncol(MCOV.rep), by = 2), comprep, MCOV.rep)
deltaMCOV.avg = rowMeans(deltaMCOV, na.rm = T)
plot(hexbin(log10(MCOV.avg), deltaMCOV.avg, xbins = 200), style = "colorscale")




##### Collect metadata
matinfo = read.delim(info, header = T, stringsAsFactors = F)

# target = which(MCOV.avg > 10 & TaxDepth > minTaxDepth & nbAlleles < maxAlleles & lenSNP < maxlenSNP & mono == keepmono)
target = which(TaxDepth > minTaxDepth & nbAlleles < maxAlleles & lenSNP < maxlenSNP & mono == keepmono)

diagplot(log10(MCOV.avg), log10(NALL.avg), nspl = 200)
target = which(nbAlleles < maxAlleles & mono == keepmono)
points(log10(MCOV.avg[target]), log10(NALL.avg[target]), cex = .5, col = "red", pch = 16)
target = which(TaxDepth > .8 & nbAlleles < maxAlleles &lenSNP < 5 & mono == keepmono)
points(log10(MCOV.avg[target]), log10(NALL.avg[target]), cex = .5, col = "blue", pch = 16)
target = which(TaxDepth > .8 & nbAlleles < maxAlleles &lenSNP < 5 & mono == keepmono & MCOV.avg > 10)
points(log10(MCOV.avg[target]), log10(NALL.avg[target]), cex = .5, col = "green", pch = 16)



diagplot(log10(MCOV.avg), deltaHet.avg, nspl = 50)

diagplot(log10(MCOV.avg), deltaNALL.avg, nspl = 50)

plot(hexbin(log10(as.vector(MCOV)), log10(as.vector(NALL)), xbins = 100))

diagplot(log10(MCOV.avg), deltaMCOV.avg, mean, nspl = 25, smpl = length(deltaNALL.avg))


diagplot(log10(MCOV.avg), log10(NALL.avg), nspl = 100, median)
points(log10(MCOV.avg[target]), log10(NALL.avg[target]), cex = .5, col = "red", pch = 16)


# plot(log10(MCOV.avg), log10(TaxDepth), cex = .5, col = "grey", bty = "n")
# points(log10(MCOV.avg[target]), log10(TaxDepth[target]), cex = .5, col = "red", pch = 16)
# abline(h = log10(seq(0, 1, length.out = 5)))







