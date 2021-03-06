# clear all
rm(list=ls(all=TRUE))

# Source
source("calcCorrelationParallel.R")

# Settings
# number of columns per chunk
cols.per.subset = 40

# sample data
n.rows = 1e6
n.cols = 1e2
data = matrix(rnorm(n.rows*n.cols), ncol = n.cols)

# calculate correlation matrix
corMatrix = calcCorrelationParBT(data, cols.per.subset
                                  , cor.use = "pairwise.complete.obs"
                                  , test.mode = TRUE
                                  )

# if all jobs succeed, save correlation matrix
if (length(findError()) == 0) {
  save(corMatrix, file="corMatrix.RData")
}
