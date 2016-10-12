rm(list = ls())


### Get script arguments (infolder and outfile)
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }
  
infile
outfile


### load data
data = read.delim(infile, header = F)


### loop over all blast hits and get primer positions (using quantiles)
limsraw = aggregate(data[, 7:8], by = list(data[, 2]), quantile, probs = c(.25, .5, .75), simplify = T)


### keep 25% (startpos) - 75% bounds (stoppos), make sure that orientation is fine
limsOK = NULL
for(i in 1:nrow(limsraw)){
  test = limsraw[i, 2][2] <= limsraw[i, 2][2]
  if(test == T){
    out = round(c(limsraw[i, 2][, 1], limsraw[i, 3][, 3]))
    } else {
    out = round(c(limsraw[i, 2][, 3], limsraw[i, 3][, 1]))
    }
  limsOK = rbind(limsOK, out)
  }
limsOK = data.frame(limsraw[, 1], limsOK)  
  
  
### ensure that we are not keeping an outlier blast hit, primers should be no more than 40bp long
test = apply(limsOK[, -1], 1, diff)
limsOK = limsOK[ test < 40,]


### save to output
limsOK = limsOK[ order(limsOK[, 2]),]
write.table(limsOK, file = outfile, col.names = F, row.names = F, quote = F, sep = "\t")