# Define Functions
# calculate correlation for one chunk
calcCorrelation = function(x, y, cor.use) {
  require(BBmisc)
  dataX = load2(paste0("chunks/chunk", x, ".RData"))
  if (x != y) {
    dataY = load2(paste0("chunks/chunk", y, ".RData"))
    cor(dataX, dataY, use=cor.use)
  } else {
    cor(dataX, use=cor.use)
  }
}

# reduce matrix result to the correct position in the target correlation matrix
redToCorMatrix = function(aggr, job, res, chunks) {
  aggr[chunks[[job$pars$x]], chunks[[job$pars$y]]] = res
  if (job$pars$x != job$pars$y)
    aggr[chunks[[job$pars$y]], chunks[[job$pars$x]]] = t(res)
  # return
  aggr
}

# reduce matrix result to the correct position in the target correlation matrix
redToCorMatrixBT = function(init, res, job) {
  browser()
  # noch kein Zugriff auf more.args (chunks)
  params = job$defs$pars[[1]]
  init[chunks[[job$pars$x]], chunks[[job$pars$y]]] = res
  if (job$pars$x != job$pars$y)
    aggr[chunks[[job$pars$y]], chunks[[job$pars$x]]] = t(res)
  # return
  aggr
}

# main function
calcCorrelationPar = function(reg, resources=list()
                               , data, cols.per.part
                               , cor.use="everything"
                               , test.mode=FALSE) {
  # load packages
  require(BatchJobs)
  require(Matrix)
  
  # convert data to matrix
  if (!is.matrix(data))
    data = as.matrix(data)
  
  # split data matrix in parts
  chunks = chunk(1:ncol(data), cols.per.part)
  n.chunks = length(chunks)
  # save data in chunks
  if (!file.exists("chunks")) {
    dir.create("chunks")
  }
  for (i in 1:n.chunks) {
    save2(file=paste0("chunks/chunk", i, ".RData"), dataChunk=data[, chunks[[i]]])
  }
  
  # expand grid, only for one triangular with main diagonal combinations
  dfCombi = as.data.frame(expand.grid(x=1:n.chunks, y=1:n.chunks))
  dfCombi = dfCombi[dfCombi$x >= dfCombi$y, ]
  
  # Registry for BatchJobs
  batchMap(reg, calcCorrelation, x=dfCombi$x, y=dfCombi$y, more.args=list(cor.use=cor.use))
  submitJobs(reg, resources=resources)
  
  # wait for jobs to be finished
  stopifnot(waitForJobs(reg))
  
  # collect results
  corMatrix = reduceResults(reg, fun=redToCorMatrix
                             , init=matrix(double(0), nrow=ncol(data), ncol=ncol(data))
                             , chunks=chunks
                             )
  # Test result
  print(paste("symmetric:", isSymmetric(corMatrix)))
  if (test.mode)
    print(paste("SSE:", sum((cor(data, use=cor.use)-corMatrix)^2)))
  
  # Stats
  showStatus(reg)
  
  # return
  return(corMatrix)
}

# main function
calcCorrelationParBT = function(reg.name, data, cols.per.part
                                 , cor.use = "everything"
                                 , test.mode = FALSE
                                 ) {
  # load packages
  require(batchtools)
  require(BBmisc)
  require(checkmate)
  require(Matrix)
  
  # check arguments
  assertString(reg.name)
  assert(checkMatrix(data), checkDataFrame(data))
  assertNumber(cols.per.part, lower = 2)
  # see parameter 'use' in ?cor
  assertChoice(cor.use, choices = c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs"))
  assertFlag(test.mode)
  
  # delete directory
  stopifnot(unlink(reg.name, recursive = TRUE) == 0)
  # create batchtools registry
  reg = makeRegistry(reg.name)
  reg$default.resources = list(walltime = 1800, memory = 1024)
  saveRegistry(reg)

  # convert data to matrix
  if (!is.matrix(data))
    data = as.matrix(data)
  
  # split data matrix in parts
  chunks = chunk(1:ncol(data), cols.per.part)
  n.chunks = length(chunks)
  # save data in chunks
  if (!file.exists("chunks")) {
    dir.create("chunks")
  }
  for (i in 1:n.chunks) {
    save2(file=paste0("chunks/chunk", i, ".RData"), dataChunk = data[, chunks[[i]]])
  }
  
  # expand grid, only for one triangular with main diagonal combinations
  dfCombi = as.data.frame(expand.grid(x = 1:n.chunks, y = 1:n.chunks))
  dfCombi = dfCombi[dfCombi$x >= dfCombi$y, ]
  
  # Registry for BatchJobs
  batchMap(calcCorrelation, x = dfCombi$x, y = dfCombi$y, more.args = list(cor.use = cor.use))
  submitJobs()
  
  # wait for jobs to be finished
  catf("WaitForJobs: %s", waitForJobs())
  stopifnot(waitForJobs())
  
  # collect results
  corMatrix = reduceResults(fun = redToCorMatrixBT
                             , init = matrix(double(0), nrow = ncol(data), ncol = ncol(data))
                             )
  # Test result
  print(paste("symmetric:", isSymmetric(corMatrix)))
  if (test.mode)
    print(paste("SSE:", sum((cor(data, use=cor.use)-corMatrix)^2)))
  
  # Stats
  getStatus()
  
  # return
  return(corMatrix)
}
