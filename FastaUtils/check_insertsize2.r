data = read.delim("Megablast_800.out", header=F) 

coords = data.frame('rd'=data[,1], 'sca'=data[,2], 'pos' =(data[, 10] + data[, 9])/2)
coords = aggregate(coords$pos, by=list(coords$rd,coords$sca), max)
reads = levels(as.factor(gsub('/[1-2]','',as.character(coords[,1]))))

ones = data.frame(tag = paste(reads,'/1',sep=''), NA, NA)
twoes = data.frame(tag = paste(reads,'/2',sep=''), NA, NA)

for(i in 1:nrow(coords)){
  stack = coords[i, ]
  
  idx1 = match(as.character(unlist(stack[1])), ones[,1])
  idx2 = match(as.character(unlist(stack[1])), twoes[,1])
  if(is.na(idx1) == F){
    ones[idx1, 2] = as.character(unlist(stack[2]))
    ones[idx1, 3] = stack[3]
    } 
  if(is.na(idx2) == F){
    twoes[idx2, 2] = as.character(unlist(stack[2]))
    twoes[idx2, 3] = stack[3]
    } 
  }

focus = abs(ones[,3] - twoes[,3])
scaf1 = ones[is.na(focus)==F,2]
scaf2 = twoes[is.na(focus)==F,2]

final = focus[scaf1 == scaf2]
final = final[is.na(final)==F]
c(length(final), mean(final), sd(final))
