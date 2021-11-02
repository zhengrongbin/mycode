args = commandArgs(T)

library('sva')

f = args[1] # the path of expression matrix

print('++ readin file')
exp_mat = read.csv(f, sep = '\t', row.names = 1)

# read in the expression matrix
dat = log2(exp_mat+1)

print('++ estimate varibles')
# get the number of surrogate varibles
mod = matrix(1, nrow = ncol(dat), ncol = 1) # row = gene, columns = samples in dat
num.pc.estimates = num.sv(dat, mod,method="be")
print(num.pc.estimates)

print('++ pc loading')
# compute pc loads
pc.loadings = svd(scale(dat))$u

# pc correction
n.pc <- c(1:num.pc.estimates)
print(paste("removing", n.pc, "PCs", nrow(t(dat))))
## use residuals from top n.pc principal components

pc.correct = function(df){
  n.pc <- c(1:num.pc.estimates)
  dat.adjusted <- lm(df ~ pc.loadings[,n.pc])$residuals # regress out estimated number of PCs' loads
  return(dat.adjusted)
}

print('++ correct exp')
dat.adjusted = apply(dat, 2, pc.correct)#mcmapply(pc.correct, dat, pc.loadings, num.pc.estimates)

write.table(dat.adjusted, paste0(f, '_pca_adjusted.xls'), sep = '\t', quote = F)