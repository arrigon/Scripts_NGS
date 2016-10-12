alldata = read.delim("EPH1.csv", header = T, sep = "\t")

for(i in levels(alldata[, 1])){
  subset = alldata[ alldata[, 1] == i,]
  write.table(subset[, c(3, 2)], file = paste(paste("Barcodes_EPH1", i, sep = "_"), ".txt", sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)
  }
