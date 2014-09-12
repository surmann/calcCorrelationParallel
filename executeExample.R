# clear all
rm(list=ls(all=TRUE))

# load packages
library(BatchJobs)
# Source
source("calcCorrelationParallel.R")

# Settings
# number of columns per chunk
colsPerPart <- 400
# name for BatchJobs
regName <- "calcCorrelation"

# sample data
data <- matrix(rnorm(1e6), ncol=1000)

# delete old BatchJob
if (FALSE) {
  unlink(paste0(regName, "-files"), recursive=TRUE)
}

# create BatchJobs registry
reg <- makeRegistry(regName)
# calculate correlation matrix
corMatrix <- calcCorrelationPar(reg, resources=list(walltime=1800, memory=1024)
                                , data
                                , colsPerPart
                                , corUse="pairwise.complete.obs"
                                , blnTestMode=FALSE
                             )

# if all jobs succeed, save correlation matrix
if (length(findErrors(reg)) == 0) {
  save(corMatrix, file="corMatrix.RData")
}
