# Define Functions
# calculate correlation for one chunk
calcCorrelation <- function(x, y, corUse) {
  dataX <- load2(paste0("chunks/chunk", x, ".RData"))
  if (x != y) {
    dataY <- load2(paste0("chunks/chunk", y, ".RData"))
    cor(dataX, dataY, use=corUse)
  } else {
    cor(dataX, use=corUse)
  }
}

# reduce matrix result to the correct position in the target correlation matrix
redToCorMatrix <- function(aggr, job, res, chunks) {
  aggr[chunks[[job$pars$x]], chunks[[job$pars$y]]] <- res
  if (job$pars$x != job$pars$y)
    aggr[chunks[[job$pars$y]], chunks[[job$pars$x]]] <- t(res)
  # return
  aggr
}

# main function
calcCorrelationPar <- function(reg, resources=list()
                               , data, colsPerPart
                               , corUse="everything"
                               , blnTestMode=FALSE) {
  # load packages
  require(BatchJobs)
  require(Matrix)
  
  # convert data to matrix
  if (!is.matrix(data))
    data <- as.matrix(data)
  
  # split data matrix in parts
  chunks <- chunk(1:ncol(data), colsPerPart)
  nChunks <- length(chunks)
  # save data in chunks
  if (!file.exists("chunks")) {
    dir.create("chunks")
  }
  for (i in 1:nChunks) {
    save2(file=paste0("chunks/chunk", i, ".RData"), dataChunk=data[, chunks[[i]]])
  }
  
  # expand grid, only for one triangular with main diagonal combinations
  dfCombi <- as.data.frame(expand.grid(x=1:nChunks, y=1:nChunks))
  dfCombi <- dfCombi[dfCombi$x >= dfCombi$y, ]
  
  # Registry for BatchJobs
  batchMap(reg, calcCorrelation, x=dfCombi$x, y=dfCombi$y, more.args=list(corUse=corUse))
  submitJobs(reg, resources=resources)
  
  # wait for jobs to be finished
  stopifnot(waitForJobs(reg))
  
  # collect results
  corMatrix <- reduceResults(reg, fun=redToCorMatrix
                             , init=matrix(double(0), nrow=ncol(data), ncol=ncol(data))
                             , chunks=chunks
                             )
  # Test result
  print(paste("symmetric:", isSymmetric(corMatrix)))
  if (blnTestMode)
    print(paste("SSE:", sum((cor(data, use=corUse)-corMatrix)^2)))
  
  # Stats
  showStatus(reg)
  
  # return
  return(corMatrix)
}
