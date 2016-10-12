# infolder = "../tmp/pileups"
# outfile = "../data.out/AssemblyStats.pdf"

### Get script arguments
args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
  }


### load coverage files
allcovs = dir(infolder, pattern = ".cov$", full.names = T)
ALL = NULL
for(file in allcovs){
  bsn = basename(file)
  bsn = gsub(".cov", "", bsn)
  data = read.delim(file = file, header = F)
  data[, 1] = rep(bsn, nrow(data))
  ALL = rbind(ALL, data)
  }
ALL[is.na(ALL[, 3]), 3] = 0


## Compute coverage stats
pdf(file = outfile)
  hist(log(ALL[, 3] + 1), col = "dark red", breaks = 20, axes = F, xlab = "Coverage (nreads / bp)", main = "Assembly coverage")
  axis(side = 2)
  xlims = ceiling(range(log(ALL[, 3] + 1)))
  ticks = seq(xlims[1], xlims[2], by = 1)
  labs = 10^ticks
  axis(side = 1, at = ticks, labels = labs)

  # View local coverage
  local = aggregate(ALL[, 3], by = list(ALL[, 2]), median)
  plot(local, type = "l", col = "dark red", xlab = "Assembly position (bp)", ylab = "median coverage (x)", bty = "n")
  wng = local[ local[, 2] <= 1,]
  if(nrow(wng) > 1) points(wng, pch = 16, col = "red", cex = 0.5)
  abline(h = median(ALL[, 3]), lty = 2)

  # view coverage of each specimen
  for(dset in unique(ALL[, 1])){
    set = ALL[ ALL[, 1] == dset, 2:3]
    plot(set, type = "l", col = "dark red", xlab = "Assembly position (bp)", ylab = "median coverage (x)", bty = "n", main = paste("Coverage for", dset))
    wng = set[ set[, 2] <= 1,]
    if(nrow(wng) > 1) points(wng, pch = 16, col = "red", cex = 0.5)
    abline(h = median(set[, 2]), lty = 2)
    }

  dev.off()
