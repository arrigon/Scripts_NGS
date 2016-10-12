alldata = read.delim("MAS_lib1.csv", header = T, sep = "\t")

for(i in 1:max(alldata[, 1])){
  subset = alldata[ alldata[, 1] == i,]
  for(prefix in c("RAD", "BIS")){
    subset$Tag = paste(prefix, subset[, 3], sep = "_")
    write.table(subset[, c(4, 2)], file = paste(paste("Barcodes_MAS", prefix, i, sep = "_"), ".txt", sep = ""), col.names = F, row.names = F, sep = "\t", quote = F)
    }
  }
