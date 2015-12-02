Calculate Correlation Matrix in Parallel
----------------------------------------------------------------------------------------------

The code provides a function to calculate a correlation matrix out of huge data matrix.
It splits the data matrix into subsets and calculates the correlation between the subsets in parallel using [BatchJobs](https://github.com/tudo-r/BatchJobs) or [batchtools](https://github.com/mllg/batchtools).
The batchtools implementation is marked with 'BT'.
