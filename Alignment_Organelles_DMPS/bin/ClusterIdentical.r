require(ape)

infile = "Pyrgus_BoldGenbankWholeMitoDMPS.fix.aln"

alnmt = read.dna(infile, format = "fasta")
alnmt.d = dist.dna(alnmt, model = "N", pairwise.deletion = T)


### regroup identical sequences
clusts = cutree(hclust(alnmt.d, "ward"), h = 0)


### get their respective lengths
cntchar = function(x){
  length(x[x != "n" & x != "N"])
  }

test = as.character(alnmt)
lens = apply(test, 1, cntchar)
  

### Pickup one seq (the longest) per cluster
index = data.frame(clust = clusts, len = lens, acc = rownames(alnmt), stringsAsFactors = F)
index = index[order(index[, 1], index[, 2]),]

keepers = NULL
for(i in 1:max(clusts)){
  tmp = index[ index[,1] == i,]
  out = tmp[ which.max(tmp[, 2])[1], 3]
  keepers = c(keepers, out)
  }
index.keepers = data.frame(clust = 1:length(keepers), acc = keepers)

write.table(index, file = paste("clusters_", infile), quote = F, sep = "\t")
write.table(index.keepers, file = paste("clusters.kept_", infile), quote = F, sep = "\t")


alnmt.prune = alnmt[match(keepers, rownames(alnmt)),]
write.dna(alnmt.prune, file = paste("pruned_", infile), format = "fasta", nbcol = -1)







