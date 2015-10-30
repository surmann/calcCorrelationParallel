# clear all
rm(list=ls(all=TRUE))

# load packages
library(batchtools)
# Source
source("calcCorrelationParallel.R")

# Settings
# number of columns per chunk
cols.per.part = 40
# name for BatchJobs
reg.name = "calcCorrelationBT"

# sample data
data = matrix(rnorm(1e6), ncol=100)

# calculate correlation matrix
corMatrix <- calcCorrelationParBT(reg.name, data
                                  , cols.per.part
                                  , cor.use = "pairwise.complete.obs"
                                  , test.mode = FALSE
                                  )

# if all jobs succeed, save correlation matrix
if (length(findError()) == 0) {
  save(corMatrix, file="corMatrix.RData")
}
