library(nimble)
source('benchmark_suite.R')

#'  By default, benchmarkSuite() runs three different models:
#'  MCMCtest_1 is the 'seeds' BUGS example model
#'  MCMCtest_2 is the bivariate normal 'birats' BUGS example model
#'  PFtest_1 is a univariate linear gaussian state space model with 100 time points
#'  
#'  To run any subset of these models, the arguments MCMCtests and PFtests.
#'  For example, benchmarkSuite(MCMCtests = "MCMCtest_1", PFtests = NULL) 
#'  will run only the first MCMC test (the 'seeds' model)

## run with basic settings 
tst <- benchmarkSuite(niter = 10000, nparticles = 100000, nreps = 10, summarize = T, debug = F)
  
## only run MCMC benchmarks, using different numbers of iterations for each benchmark
tst <- benchmarkSuite(MCMCtests = NULL, niter = c(100000, 10000), nreps = 10, summarize = T, debug = F)

## specify PF and MCMC options
tst <- benchmarkSuite(MCMCs = 'nimble_RW', PFs = 'auxiliary', niter = 10000, nreps = 10, summarize = T, debug = F)
