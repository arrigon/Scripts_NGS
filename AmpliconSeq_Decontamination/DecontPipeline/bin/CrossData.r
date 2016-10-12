args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }

## WARNING: hard coded params:
keywords = c("cytochrome", "ribosomal")
family = "Ephemeroptera"


### Debug (comment these lines while in production)
# intaxo = "MG0037.fq_1.taxo"
# inclst = "MG0037.fq_1.fas.bak.clstr"
# inannot = "MG0037.fq_1.annot"
# inhits = "MG0037.fq_1.gb"
# inblast = "MG0037.fq_1.blast"
# outfolder = "tmp/"


## Prepare working envmt.
bsn = basename(inhits)
bsn = gsub("\\.gb$", "", bsn, perl = T)
outstats = paste(outfolder, "/", bsn, ".stats", sep = "")
outlist = paste(outfolder, "/", bsn, ".cleanreads", sep = "")


## load data
# read clusters and centroids
clust = read.delim(file = inclst, header = F, sep = "\t", stringsAsFactors = F)
clust[, 2] = gsub("^.*, >", "", clust[, 2])
centroids = grep("\\*", clust[, 2])
clust[, 2] = gsub("\\.\\.\\..*$", "", clust[, 2])
centroids = clust[ centroids,]

# Blast hits (check also if we need to revcomp the reads
blast = read.delim(file = inblast, header = F, sep = "\t", stringsAsFactors = F)
blast[, 2] = read.delim(textConnection(blast[, 2]), sep = "|", stringsAsFactors = F, header = F)[, 4]
revcomp = blast[, 9] > blast[, 10]
revcomp[revcomp == T] = "rc"
revcomp[revcomp == F] = "ok"
blast = blast[, c(1:4, 11)]
blast[, 6] = revcomp

# taxo
taxo = read.delim(file = intaxo, header = F, sep = "\t", stringsAsFactors = F)

# GB annotations
annot = read.delim(file = inannot, header = F, sep = "\t", stringsAsFactors = F)
getfirst = function(x){
  levs = unique(x)
  match(levs, x)
  }
annot.unique = annot[getfirst(annot[, 1]), ]


## link taxo, hits and annot
taxo.OK = taxo
centroids.OK = centroids[ match(blast[, 1], centroids[, 2]), ]
annot.OK = annot.unique[ match(blast[, 2], annot.unique[, 1]),]
metadata = data.frame(taxo.OK, 
		      annot.OK, 
		      centroids.OK, 
		      read.table(textConnection(taxo.OK[, 4]), sep = ";", stringsAsFactors = F),
		      revcomp = blast[, 6],
		      eval = blast[, 5], stringsAsFactors = F)


## hunt for keywords (markers)
metadata$marker = rep("NA", nrow(metadata))
for(word in keywords){
  tmp = grep(word, metadata[, 6])
  if(length(tmp) > 0) metadata[tmp, ]$marker = word
  }
metadata = metadata[ !is.na(match(metadata$marker, keywords)),]  


## clean double hits (keep best hit, or the one matching against family of interest)
FINAL = NULL
for(i in unique(metadata[, 1])){
  tmp = metadata[ metadata[, 1] == i,]
  if(nrow(tmp) > 1){
    taxo = tmp[, 12]
    if(any(taxo == family)) {
      tmp = tmp[ taxo == family,]
      if(nrow(tmp) > 1){
	eval = tmp$eval
	tmp = tmp[which.min(eval),]
	}
      } else {
      eval = tmp$eval
      tmp = tmp[which.min(eval),]
      }
    }
  FINAL = rbind(FINAL, tmp)
  }
metadata = FINAL



## Count cluster coverages 
covs = table(clust[, 1])
metadata$cov = covs[match(as.character(metadata[, 7]), names(covs))]

## Summary stats
crossval = tapply(metadata$cov, list(metadata[, 12], metadata$marker), sum, na.rm = T)
crossval[is.na(crossval)] = 0
barplot((t(crossval)), las = 2)
write.table(crossval, file = outstats, sep = "\t", col.names = T, quote = F, row.names = T)



## Extract reads of interest
metadata.target = metadata[ metadata[, 12] == family,]

clust.info = metadata.target[match(clust[, 1], metadata.target[, 7]), ]
clust.info = data.frame(clust, clust.info)
clust.info = clust.info[ !is.na(clust.info[, 3]),]

reads = clust.info[, 2]
marker = clust.info$marker
rc = clust.info$revcomp
write.table(cbind(marker, reads, rc), file = outlist, sep = "\t", col.names = F, row.names = F, quote = F)


