args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }

## WARNING: hard coded params:
keywords = c("cytochrome", "ribosomal")
family = "Ephemeroptera"


### Debug (comment these lines while in production)
# intaxo = "MG0128.fq_1.taxo"
# inclst = "MG0128.fq_1.fas.bak.clstr"
# inannot = "MG0128.fq_1.annot"
# inhits = "MG0128.fq_1.gb"
# inblast = "MG0128.fq_1.blast"
# outfolder = "tmp/"


## Prepare working envmt.
bsn = basename(inhits)
bsn = gsub("\\.gb$", "", bsn, perl = T)
outstats = paste(outfolder, "/", bsn, ".stats", sep = "")
outlist = paste(outfolder, "/", bsn, ".cleanreads", sep = "")


## load data
# taxo
taxo = read.delim(file = intaxo, header = F, sep = "\t", stringsAsFactors = F)

# read clusters and centroids
clust = read.delim(file = inclst, header = F, sep = "\t", stringsAsFactors = F)
clust[, 2] = gsub("^.*, >", "", clust[, 2])
centroids = grep("\\*", clust[, 2])
clust[, 2] = gsub("\\.\\.\\..*$", "", clust[, 2])
centroids = clust[ centroids,]


# GB annotations
hits = read.delim(file = inhits, header = F, sep = "\t", stringsAsFactors = F)
annot = read.delim(file = inannot, header = F, sep = "\t", stringsAsFactors = F)
getfirst = function(x){
  levs = unique(x)
  match(levs, x)
  }
annot.unique = annot[getfirst(annot[, 1]), ]
annot = annot.unique[ match(hits[, 1], annot.unique[, 1]),]


# Blast hits (check also if we need to revcomp the reads
blast = read.delim(file = inblast, header = F, sep = "\t", stringsAsFactors = F)
blast[, 2] = read.delim(textConnection(blast[, 2]), sep = "|", stringsAsFactors = F, header = F)[, 4]
revcomp = blast[, 9] > blast[, 10]
revcomp[revcomp == T] = "c"
revcomp[revcomp == F] = "ok"
blast = blast[, 1:2]
blast[, 3] = revcomp



## link taxo, hits and annot
taxo.OK = taxo[ match(blast[, 1], taxo[, 1]), ]
annot.OK = annot[ match(blast[, 2], annot[, 1]), ]
centroids.OK = centroids[ match(blast[, 1], centroids[, 2]), ]
metadata = data.frame(taxo.OK, 
		      annot.OK, 
		      centroids.OK, 
		      read.table(textConnection(taxo.OK[, 4]), sep = ";", stringsAsFactors = F),
		      revcomp = blast[, 3], stringsAsFactors = F)


## hunt for keywords (markers)
metadata$marker = rep("NA", nrow(metadata))
for(word in keywords){
  tmp = grep(word, metadata[, 6])
  if(length(tmp) > 0) metadata[tmp, ]$marker = word
  }
metadata = metadata[ !is.na(match(metadata$marker, keywords)),]  


## Count cluster coverages 
covs = table(clust[, 1])
metadata$cov = covs[match(metadata[, 7], names(covs))]


## Extract reads of interest
metadata.target = metadata[ metadata[, 12] == family,]

clust.info = metadata.target[match(clust[, 1], metadata.target[, 7]), ]
clust.info = data.frame(clust, clust.info)
clust.info = clust.info[ !is.na(clust.info[, 5]),]

reads = clust.info[, 2]
marker = clust.info$marker
rc = clust.info$revcomp
write.table(cbind(marker, reads, rc), file = outlist, sep = "\t", col.names = F, row.names = F, quote = F)


## Summary stats
crossval = tapply(metadata$cov, list(metadata[, 12], metadata$marker), sum, na.rm = T)
crossval[is.na(crossval)] = 0
barplot((t(crossval)), las = 2)
write.table(crossval, file = outstats, sep = "\t", col.names = T, quote = F, row.names = T)
