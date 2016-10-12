args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }

### Script IO
# infolder = ...
# outfolder = ...


### Find data into fastqc reports
# Get file names
command = paste("find", infolder, "-type f -name 'fastqc_data*'")
xx = system(command, intern = T)
xx = dirname(xx)
xx = gsub("./complete.", "", xx)
xx = gsub("_fastqc", "", xx)
xx= basename(xx)

# Get read numbers
command = paste("find", infolder, "-type f -name 'fastqc_data*' -exec grep \'Total Sequences\' {} \\;")
yy = system(command, intern = T)
yy = gsub("Total Sequences", "", yy)
yy = as.numeric(as.character(gsub("\t", "", yy)))

# Produce output
OUT = data.frame(folder = infolder, specimen = xx, nreads = yy)
OUT = OUT[ order(OUT[, 1]),]

write.table(OUT, file = paste(outfolder, "/ReadStats.txt", sep = ""), quote = F, col.names = T, row.names = T, sep = "\t")

