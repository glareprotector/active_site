chain_file = "../chains_to_use.txt"

data = read.csv('../CSA_2_2_12.dat')
data_lit = data[data$EVIDENCE.TYPE=='LIT',]

u = unique(data_lit[,c(1,4)])
us = cbind(as.character(u[,1]), as.character(u[,2]))

num_u = dim(us)[1]

num_to_sample = 200

to_write = us[sample(num_u, num_to_sample),]
ans = c()

for(i in 1:num_to_sample){
  if(to_write[i,2] != ""){
    ans = c(ans, paste(to_write[i,1], to_write[i,2], sep='_'))
  }
}

write.table(ans, chain_file, quote=F, row.name=F, col.name=F)