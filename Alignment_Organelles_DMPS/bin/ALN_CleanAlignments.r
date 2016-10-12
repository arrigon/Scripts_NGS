rm(list = ls())


### Get script arguments (infolder and outfile)
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }
  

### DEBUG (comment when script is in use)
# infile = "../tmp/alnIII/AllMito.dustsIII.fas"
# outfile = "test.fas"

### Params
width = 1
cDEG = 0.99
minPres = .80





infile
outfile
getwd()


### load packages
require(ape)


# load datasets
data = read.dna(infile, format = "fasta", as.character = T)


# estimate proportion of missing and degenerated data
scanfun = function(x){
  x = toupper(x)
  x = table(x)

  MIS = sum(x[grep("[N|-]", names(x))])
  DEG = sum(x[grep("[^A|T|C|G|N|-]", names(x))])
  
  propMIS = MIS/sum(x)
  propDEG = DEG/sum(x)
  
  c(propMIS, propDEG)
  }

stats = apply(data, 2, scanfun)
stats = rbind(1:ncol(stats), stats)
  

# Sliding window to mesure proportion of degenerated sites
wMIS = rep(0, ncol(stats))
wDEG = rep(0, ncol(stats))
for(i in (1+width):(ncol(stats) - (width + 1))){
  focus = (i-width):(i+width)
  MIS = mean(stats[2, focus])
  DEG = mean(stats[3, focus])
  wMIS[focus] = MIS
  wDEG[focus] = DEG
  }


# Identify loci with degenerated sites
cutoff = quantile(wDEG, probs = cDEG)
sDEG = which(wDEG > cutoff)


# Identify sites with missing data
sMIS = which(stats[2 ,] >= minPres)


# Discard these sites from alignment
tormve = sort(unique(c(sMIS, sDEG)))
data.final = data[, -tormve]


# save to outfasta
cat(file = outfile, append = F)
for(i in 1:nrow(data.final)){
  acc = rownames(data.final)[i]
  seq = paste(data.final[i, ], collapse = "")
  seq = toupper(seq)
  cat(">", acc, "\n", seq, "\n", file = outfile, append = T, sep = "")
  }


