# Define Functions
# calculate correlation for one subset
calcCorrelation = function(x, y, cor.use) {
  require(BBmisc)
  dataX = load2(paste0("dataSubsets/subset", x, ".RData"))
  if (x != y) {
    dataY = load2(paste0("dataSubsets/subset", y, ".RData"))
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
redToCorMatrixBT = function(init, res, job, subsets) {
  # noch kein Zugriff auf more.args (chunks)
  params = job$defs$pars[[1]]
  init[subsets[[params$x]], subsets[[params$y]]] = res
  if (params$x != params$y)
    init[subsets[[params$y]], subsets[[params$x]]] = t(res)
  # return
  init
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
  if (!file.exists("dataSubsets")) {
    dir.create("dataSubsets")
  }
  for (i in 1:n.chunks) {
    save2(file=paste0("dataSubsets/subset", i, ".RData"), dataChunk=data[, chunks[[i]]])
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
calcCorrelationParBT = function(data, cols.per.subset
                                , cor.use = "everything"
                                , test.mode = FALSE
                                , reg.name = "calcCorrelationBT"
                                , resources = list(walltime = 1800, memory = 1024)
                                ) {
  # load packages
  require(batchtools)
  require(BBmisc)     # chunk
  require(checkmate)
  
  # check arguments
  assert(checkMatrix(data), checkDataFrame(data))
  assertNumber(cols.per.subset, lower = 2)
  # see parameter 'use' in ?cor
  assertChoice(cor.use, choices = c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs"))
  assertFlag(test.mode)
  assertString(reg.name)
  assertList(resources, len = 2)
  
  # delete directory
  stopifnot(unlink(reg.name, recursive = TRUE) == 0)
  # create batchtools registry
  reg = makeRegistry(reg.name)
  reg$default.resources = resources
  saveRegistry(reg)

  # convert data to matrix
  if (!is.matrix(data))
    data = as.matrix(data)
  
  # split data matrix in subsets
  subsets = chunk(1:ncol(data), chunk.size = cols.per.subset)
  n.subsets = length(subsets)
  # save data in subsets
  if (!file.exists("dataSubsets")) {
    dir.create("dataSubsets")
  }
  for (i in 1:n.subsets) {
    save2(file=paste0("dataSubsets/subset", i, ".RData"), dataChunk = data[, subsets[[i]]])
  }

  # expand grid, only for one triangular with main diagonal combinations
  dfCombi = as.data.frame(expand.grid(x = 1:n.subsets, y = 1:n.subsets))
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
                            , subsets = subsets
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
