#' Executes multiple MCMC algorithms and organizes results.
#'
#' Creates, runs, and organizes output from a suite of MCMC algorithms, all applied to the same model, data, and initial values.
#' This can include WinBUGS, OpenBUGS, JAGS and Stan MCMCs, as well as NIMBLE MCMC algorithms.
#' Trace plots and density plots for the MCMC samples may also be generated and saved.
#'
#' @details
#' Creates and runs an MCMC Suite.
#' By default, this will execute the specified MCMCs, record all samples, generate summary statistics, and create and save trace plots and posterior density plots.
#' This default behavior can ben altered via a variety of arguments.
#' Following execution of the MCMC algorithms, returns a named list containing \code{samples}, \code{summary}, and \code{timing} elements.
#' See the NIMBLE User Manual for more information about the organization of the return object.
#' 
#' @param code The quoted code expression representing the model, such as the return value from a call to \code{nimbleCode}).
#' No default value, this is a required argument.
#'
#' @param constants A named list giving values of constants for the model.
#' This is the same as the \code{constants} argument which would be passed to \code{nimbleModel}.
#' Default value is list().
#'
#' @param data A named list giving the data values for the model.
#' This is the same as the \code{data} argument which would be passed to \code{nimbleModel} or \code{model$setData}.
#' Default value is \code{list()}.
#' 
#' @param inits A named list giving the initial values for the model.
#' This is the same as the \code{inits} argument which would be passed to \code{nimbleModel} or \code{model$setInits}.
#' Default value is \code{list()}.
#'
#' @param monitors A character vector giving the node names or variable names to monitor.
#' The samples corresponding to these nodes will be stored in the output samples, will have summary statistics calculated, and density and trace plots generated.
#' Default value is all top-level stochastic nodes of the model.
#' 
#' @param niter Number of MCMC iterations to run.
#' This applies to all MCMC algorithms in the suite.
#' Default value is 10,000.
#'
#' @param burnin Number of initial, post-thinning, MCMC iterations to discard.
#' Default value is 2,000.
#' 
#' @param thin Thinning interval for the MCMC samples.
#' This applies to all MCMC algorithms in the suite.  The thinning occurs prior to the burnin samples being discarded.
#' Default value is 1.
#' 
#' @param summaryStats A character vector, specifying the summary statistics to calculate on the MCMC samples.
#' Each element may be the character name of an exisiting R function (possibly user-defined) which acts on a numeric vector and returns a scalar (e.g., \code{mean} or \code{sd},
#' or a character string which when parsed and evaluted will define such a function (e.g., \code{function(x) mean(sqrt(x))}).
#' Default value is \code{c('mean', 'median', 'sd', 'CI95_low', 'CI95_upp')}, where the final two elements are functions which calculate the limits of a 95 percent Bayesian credible interval.
#' 
#' @param calculateEfficiency A logical, specifying whether to calculate the efficiency for each MCMC algorithm.  Efficiency is defined as the effective sample size (ESS) of each model parameter divided by the algorithm runtime (in seconds).  Default is FALSE.
#'
#' @param MCMCs A character vector specifying the MCMC algorithms to run.
#' \code{'winbugs'} specifies WinBUGS;
#' \code{'openbugs'} specifies OpenBUGS;
#' \code{'jags'} specifies JAGS;
#' \code{'stan'} specifies Stan; in this case, must also provide the \code{'stan_model'} argument;
#' \code{'nimble'} specifies NIMBLE's default MCMC algorithm;
#' \code{'nimble_noConj'} specifies NIMBLE's default MCMC algorithm without the use of any conjugate Gibbs sampling;
#' \code{'nimble_RW'} specifies NIMBLE MCMC algorithm using only random walk Metropolis-Hastings (\code{'RW'}) samplers;
#' \code{'nimble_slice'} specifies NIMBLE MCMC algorithm using only slice (\code{'slice'}) samplers;
#' \code{'autoBlock'} specifies NIMBLE MCMC algorithm with block sampling of dynamically determined parameter groups attempting to maximize sampling efficiency;
#' Anything else will be interpreted as NIMBLE MCMC algorithms, and must have associated entries in the MCMCdefs argument.
#' Default value is \code{'nimble'}, which specifies NIMBLE's default MCMC algorithm.
#' 
#' @param MCMCdefs A named list of MCMC definitions.  The names of list elements should corespond to any custom MCMC algorithms specified in the \code{MCMCs} argument.
#' The list elements should be quoted expressions, enclosed in {} braces.  When executed, the internal code must return an MCMC configuration object, 
#' specifying the corresponding MCMC algorithm; in particular, setting the appropriate samplers.  The code may assume existance of the R model object \code{Rmodel},
#' and must *return* the MCMC configuration object.  Therefore, the final line of such a code block would frequently be a standalone \code{MCMCconf}, to return this object.
#' 
#' @param winbugs_directory A character string giving the directory of the executable WinBUGS program for the WinBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'C:/WinBUGS14'}.
#' 
#' @param winbugs_program A character string giving the name of the WinBUGS program, for the WinBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'WinBUGS'}.
#'
#' @param openbugs_directory A character string giving the directory of the executable OpenBUGS program for the OpenBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'C:/OpenBUGS323'}.
#' 
#' @param openbugs_program A character string giving the name of the OpenBUGS program, for the OpenBUGS MCMC.
#' This argument will be passed directly to the bugs(...) call, from the R2WinBUGS library.
#' Default value is \code{'OpenBUGS'}.
#' 
#' @param stan_model A character string specifying the location and name of the model file (\code{'modelName.stan'}) for use with the Stan MCMC program.
#' This argument must include the \code{'.stan'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.
#' 
#' @param stan_inits A character string specifying the location and name of the inits file (\code{'modelName.init.R'}) for use with the Stan MCMC program.
#' This argument must include the \code{'.init.R'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.
#' If omitted, it will attempt to locate an inits file in the same directory as the Stan model file.
#'
#' @param stan_data A character string specifying the location and name of the data file (in the form \code{'modelName.data.R'}) for use with the Stan MCMC program.
#' This argument must include the \code{'.data.R'} extension, and must be provided whenever the \code{MCMCs} argument includes \code{'stan'}.
#' If omitted, it will attempt to locate a data file in the same directory as the Stan model file.
#'
#' @param stanNameMaps A list specifying name mappings between Stan and WinBUGS/OpenBUGS.
#' The syntax for list elements is list(BUGS_PARAM_NAME = list(StanSourceName = 'STAN_PARAM_NAME', transform = function(x) TRANSFORMATION_FUNCTION(x))).
#' The transformation is optional.
#' 
#' @param makePlot Logical argument, specifying whether to generate the trace plots and posterior density plots, for each monitored node.
#' Default value is \code{TRUE}.
#' 
#' @param savePlot Logical argument, specifying whether to save the trace plots and density plots.
#' Plots will be saved into the current working directory.
#' Only used when \code{makePlot == TRUE}.
#' Default value is \code{TRUE}.
#' 
#' @param plotName Character string, giving the file name for saving the trace plots and density plots.
#' Only used when \code{makePlot == TRUE} and \code{savePlot == TRUE}.
#' Default value is \code{'MCMCsuite'}.
#'
#' @param setSeed Logical argument, specifying whether to set.seed(0) prior to MCMC sampling.
#' Default value is \code{TRUE}.
#' 
#' @param check Logical argument, specifying whether to check the model object for missing or invalid values.  Default is given by the NIMBLE option 'checkModel', see help on \code{nimbleOptions} for details.
#' 
#' @param debug Logical argument, specifying whether to enter a \code{browser()} at the onset of executing each MCMC algrithm.
#' For use in debugging individual MCMC algorithms, if necessary.
#' Default value is FALSE.
#'
#' @param ... For internal use only
#'
#' @return Returns a named list containing elements:
#' samples: A 3-dimensional array containing samples from each MCMC algorithm.
#' summary: A 3-dimensional array containing summary statistics for each variable and algorithm.
#' timing: A numeric vector containing timing information.
#' efficiency: Minimum and mean sampling efficiencies for each algorithm (only provided if option calculateEfficiency = TRUE).
#' See the NIMBLE User Manual for more information about the organization of the return object.
#'  
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#' })
#' output <- MCMCsuite(code,
#'                     data = list(x=3),
#'                     inits = list(mu=0),
#'                     niter = 10000,
#'                     monitors = 'mu',
#'                     MCMCs = c('nimble', 'nimble_RW'),
#'                     summaryStats = c('mean', 'sd', 'max', 'function(x) max(abs(x))'),
#'                     makePlot = FALSE)
#' }
#' 
#' @author Daniel Turek
#' @export
benchmarkSuite <- function(
            MCMCtests           = "all",
            MCMCs               = c('nimble', 'nimble'),
            niter               = 10000,
            PFtests             = "all",
            PFs                 = 'bootstrap',
            pfControlList       = list(),
            nparticles          = 100000,
            nreps               = 1,
            setSeed             = TRUE,
            check               = getNimbleOption('checkModel'),
            summarize           = TRUE,
            debug               = FALSE
) {
    ## aliased in MCMCsuiteClass
    suite <- nimBenchmarkClass(MCMCtests, MCMCs, niter, PFtests, PFs,
                               pfControlList, nparticles, nreps,
                               setSeed, check, summarize, debug)
    return(suite$output)
}

#' Class \code{MCMCsuiteClass}
#'
#' @aliases MCMCsuiteClass-class
#'
#' @description
#' Objects of this class create, run, and organize output from a suite of MCMC algorithms, all applied to the same model, data, and initial values.
#' This can include WinBUGS, OpenBUGS, JAGS and Stan MCMCs, as well as NIMBLE MCMC algorithms.
#' Trace plots and density plots for the MCMC samples may also be generated and saved.
#'
#' @seealso \link{MCMCsuite}
#' 
#' @author Daniel Turek
#' @export
#' @examples
#' \dontrun{
#' code <- nimbleCode({
#'     mu ~ dnorm(0, 1)
#'     x ~ dnorm(mu, 1)
#' })
#' output <- MCMCsuite(code,
#'                     data = list(x=3),
#'                     inits = list(mu=0),
#'                     niter = 10000,
#'                     monitors = 'mu',
#'                     MCMCs = c('nimble', 'nimble_RW'),
#'                     summaryStats = c('mean', 'sd', 'max', 'function(x) max(abs(x))'),
#'                     makePlot = FALSE)
#' }
#' 
nimBenchmarkClass <- setRefClass(

    Class = 'nimBenchmarkClass',
    
    fields = list(
        ## set in initialize()
        MCMCtests = 'character',   ## parsed expression for the model code; must be contained in { ... }    --- ORIGINAL ARGUMENT
        PFtests = 'character',
        RMCMCmodels = 'list',   ## Rmodel object
        RPFmodels = 'list',
        ## setMonitors()
        monitors = 'list',    ## the original character vector argument to initialize()    --- ORIGINAL ARGUMENT --- SLIGHTLY MODIFIED
        # monitorVars = 'character',    ## character vector of VARIABLE names of parameters to save
        # monitorNodesNIMBLE = 'character',  ## character vector of the monitor node names, with spaces as in nimble: 'y[1, 1]'
        # monitorNodesBUGS = 'character',    ## same as monitorNodes, except for WinBUGS and OpenBUGS: no spaces in node names: 'y[1,1]'
        # nMonitorNodes = 'numeric',   ## number of monitorNodes
        
        ## set in initialize()
        niter = 'numeric',    ## number of MCMC iterations to run    --- ORIGINAL ARGUMENT
        burnin = 'numeric',   ## burn-in period, the number of initial samples to discard, prior to thinning    --- ORIGINAL ARGUMENT
        thin = 'numeric',   ## thinning interval    --- ORIGINAL ARGUMENT
        nkeep = 'numeric',   ## number of samples we'll keep. equal to (niter/thin - burnin)
        burninFraction = 'numeric',  ## fraction of total sampling effort spent on burnin (burnin / (nkeep + burnin))
        
        ##
        nparticles = 'numeric',
        pfControlList = 'list',
        
        ## setMCMCs()
        MCMCs = 'character',   ## character vector of the MCMC analyses.  'winbugs', 'openbugs', 'jags', 'stan', or anything else is nimble    --- ORIGINAL ARGUMENT
        PFs   = 'character',
        latents = 'character',
        ## setMCMCdefs()
        MCMCdefs = 'list',   ## named list of {} expression code blocks, corresponding the setup for nimble MCMCs    --- ORIGINAL ARGUMENT --- SLIGHTLY MODIFIED
        MCMCdefNames = 'list',   ## names of the MCMCdefs list
        PFdefs = 'list',
        PFdefNames = 'list',
        setSeed = 'logical',   ## whether to setSeed(0) prior to running each algorithm    --- ORIGINAL ARGUMENT
        debug = 'logical',   ## whether to enter browser() before running each algorithm    --- ORIGINAL ARGUMENT
        nreps = 'numeric',
        runMCMCs = 'logical',
        runPFs   = 'logical',
        ## set in run()
        output = 'list'   ## list of numeric outputs: samples, summary, timing
    ),
    
    methods = list(
        initialize = function(
          MCMCtestArgs        = 'all',
          MCMCs               = 'nimble',
          niter               = 10000,
          PFtestArgs          = 'all',
          PFs                 = "bootstrap",
          pfControlList       = list(),
          nparticles          = 100000,
          nreps               = 1,
          setSeed             = TRUE,
          check               = getNimbleOption('checkModel'),
          summarize           = TRUE,
          debug               = FALSE) {
          
          runMCMCs <<- if(is.null(MCMCtestArgs)) FALSE else TRUE
          runPFs <<-  if(is.null(PFtestArgs)) FALSE else TRUE
          
          if(debug) browser()
          if(runMCMCs){
            MCMCtests <<- MCMCtestArgs
            if(MCMCtestArgs == 'all')
              MCMCtests <<- c("MCMCtest_1", "MCMCtest_2")  ## need better test names
            if(length(MCMCs == 1)) MCMCs <<- rep(MCMCs, length(MCMCtests))
            if(length(niter) == 1) niter <<- rep(niter, length(MCMCtests))
            
            RMCMCmodels <<- list()
            monitors <<- list() ## in future, make customizable monitors
            for(i in seq_along(MCMCtests)){
              RMCMCmodels[[i]] <<- nimbleBenchmarkModel(MCMCtests[i])
              monitors[[i]] <<- RMCMCmodels[[i]]$getNodeNames(topOnly = TRUE,  ## in future, make monitors only for MCMC benchmarks
                                                              stochOnly = TRUE)
            }
            setMCMCdefs()
          }
            
          if(runPFs){
            PFtests <<- PFtestArgs
            if(PFtestArgs == 'all')
              PFtests <<- c("PFtest_1")  ## need better test names
            
            if(length(PFs == 1)) PFs <<- rep(PFs, length(PFs))
            if(length(pfControlList) <= 1) {
              tmpList <- list()
              for(i in 1:length(PFtests))
                tmpList[[i]] <- pfControlList
              
              pfControlList <<- tmpList
            }
            if(length(nparticles) == 1) nparticles <<- rep(nparticles, length(PFtests))
            
            RPFmodels <<- list()
            for(i in seq_along(PFtests)){
              RPFmodels[[i]] <<- nimbleBenchmarkModel(PFtests[i])
              latents[[i]] <<- nimbleLatentName(PFtests[i])
            }
            
            setPFdefs()
          }
          
            nreps <<- nreps
            
            setSeed <<- setSeed
            debug <<- debug

            ## run
            init_output()
            if(debug)              browser()
            run_nimble()
            if(summarize == T) summarizeTimes()
        },
        
        nimbleBenchmarkModel = function(modelName){
          if(modelName == "MCMCtest_1"){
            dir = nimble:::getBUGSexampleDir('seeds')
            Rmodel <- readBUGSmodel("seeds", dir = dir, data = NULL, inits = NULL, useInits = TRUE,
                                    check = FALSE)
          }
          else if(modelName == "MCMCtest_2"){
            dir = nimble:::getBUGSexampleDir('birats')
            Rmodel <- readBUGSmodel("birats2", dir = dir, data = "birats-data.R", inits = "birats-inits.R", useInits = TRUE,
                                    check = FALSE)
          }
          else if(modelName == "PFtest_1"){
            code <- nimbleCode({
              x[1] ~ dnorm(mean = mu0, sd = sigma_x);
              y[1] ~ dnorm(x[1], sd=sigma_y);
              for(i in 2:N){
                x[i] ~ dnorm(mean = x[i-1], sd = sigma_x);
                y[i] ~ dnorm(mean = x[i], sd=sigma_y);
              }
              sigma_x ~ T(dnorm(1, sd = .5), 0,);
              sigma_y ~ T(dnorm(.1, sd = .5), 0,);
              mu0 <- 0
            })
            
            set.seed(0)
            
            N <- 100
            sigma_x <- 1
            sigma_y <- .1
            x <- rep(NA, N)
            y <- x
            x[1] <- rnorm(1,0,sigma_x)
            y[1] <- rnorm(1,x[1], sigma_y)
            for(i in 2:N){
              x[i] <- rnorm(1,x[i-1], sigma_x)
              y[i] <- rnorm(1,x[i], sigma_y)
            }
            consts <- list(N=N)
            testdata <- list(y=y)
            inits <- list(sigma_x=1, sigma_y=.1, x = x)
            Rmodel <- nimbleModel(code, constants = consts, data = testdata, inits = inits, check = FALSE)
          }
          else stop("not a valid benchmark test name")
          return(Rmodel)
        },
        
        nimbleLatentName = function(modelName){
          if(modelName == 'PFtest_1')
            return('x')
          else
            return(NULL)
        },
        
        setMCMCdefs = function() {
          for(i in seq_along(MCMCtests)){
            MCMCdefs[[i]] <<- list(nimble        = quote(configureMCMC(RMCMCmodels[[MCMCtest]])),
                                   nimble_noConj = quote(configureMCMC(RMCMCmodels[[MCMCtest]], useConjugacy = FALSE)),
                                   nimble_RW     = quote(configureMCMC(RMCMCmodels[[MCMCtest]], onlyRW       = TRUE)),
                                   nimble_slice  = quote(configureMCMC(RMCMCmodels[[MCMCtest]], onlySlice    = TRUE)),
                                   autoBlock     = quote(configureMCMC(RMCMCmodels[[MCMCtest]], autoBlock    = TRUE)))
            MCMCdefNames[[i]] <<- names(MCMCdefs[[i]])
          }
        },
        
        setPFdefs = function(newPFdefs) {
          for(i in seq_along(PFtests)){
            PFdefs[[i]] <<-   list(bootstrap     = quote(buildBootstrapFilter(RPFmodels[[PFtest]], nodes = latents[PFtest], 
                                                                              control = pfControlList[[PFtest]])),
                                   auxiliary     = quote(buildAuxiliaryFilter(RPFmodels[[PFtest]], nodes = latents[PFtest], 
                                                                              control = pfControlList[[PFtest]])),                                   nimble_RW     = quote(configureMCMC(Rmodel, onlyRW       = TRUE)),
                                   enkf          = quote(buildEnsembleKF(RPFmodels[[PFtest]], nodes = latents[PFtest], 
                                                                              control = pfControlList[[PFtest]]))
                                   )
            PFdefNames[[i]] <<- names(PFdefs[[i]])
          }
        },
        
        init_output = function() {
          nMCMCtests <- length(MCMCtests) 
          if(runMCMCs){
            for(i in seq_along(MCMCtests)){
              timing <- rep(NA, nreps+1)
              names(timing) <- c(paste0('MCMC_run_', 1:nreps), 'nimble_compile')
              runParams <- c(niter = niter[i]) 
              initialOutput <- list(timing=timing, runParams = runParams)
              output[[MCMCtests[i]]] <<- initialOutput
            }
          }
          if(runPFs){
            for(i in seq_along(PFtests)){
              timing <- rep(NA, nreps+1)
              names(timing) <- c(paste0('PF_run_', 1:nreps), 'nimble_compile')
              runParams <- c(pfControlList[[i]], pfType = PFs[i])
              initialOutput <- list(timing=timing, runParams = runParams)
              output[[PFtests[i]]] <<- initialOutput
            }
          }
        },
        

        run_nimble = function() {
          if(runMCMCs){
            for(MCMCtest in seq_along(MCMCtests)){
                  mcmcTag <- MCMCs[MCMCtest]
                  mcmcDef <- MCMCdefs[[MCMCtest]][[mcmcTag]]
                  mcmcConf <- eval(mcmcDef)
                  mcmcConf$addMonitors(monitors[[MCMCtest]], print = FALSE)
                  RmcmcFunction <- buildMCMC(mcmcConf)
                  timeResult <- system.time({
                      Cmodel <- compileNimble(RMCMCmodels[[MCMCtest]])
                      CmcmcFunction <- compileNimble(RmcmcFunction, project = RMCMCmodels[[MCMCtest]])
                  })
                  addTimeResult(MCMCtests[MCMCtest], 'nimble_compile', timeResult)
                  if(setSeed) set.seed(0)
                  for(runNum in 1:nreps){
                    timeResult <- system.time({ CmcmcFunction$run(niter[MCMCtest])})
                    addTimeResult(MCMCtests[MCMCtest], paste0('MCMC_run_', runNum), timeResult)
                  }
            }
          }
          if(runPFs){
            for(PFtest in seq_along(PFtests)){
                pfTag <- PFs[PFtest]
                pfDef <- PFdefs[[PFtest]][[pfTag]]
                RpfFunction <- eval(pfDef)
                timeResult <- system.time({
                  Cmodel <- compileNimble(RPFmodels[[PFtest]])
                  CpfFunction <- compileNimble(RpfFunction, project = RPFmodels[[PFtest]])
                })
                addTimeResult(PFtests[PFtest], 'nimble_compile', timeResult)
                if(setSeed) set.seed(0)
                for(runNum in 1:nreps){
                  timeResult <- system.time({ CpfFunction$run(nparticles[PFtest])})
                  addTimeResult(PFtests[PFtest], paste0('PF_run_', runNum), timeResult) 
                }
            }
          }
        },
        
        addTimeResult = function(testNum, timingType, timeResult) {
            output[[testNum]]$timing[timingType] <<- timeResult[[3]]
        },
        
        summarizeTimes = function(){
          for(i in seq_along(output)){
            output[[i]]$timing <<- c(round(summary(output[[i]]$timing[1:nreps]), 1),
                                     St.Dev. = round(sd(output[[i]]$timing[1:nreps]), 2),
                                     round(output[[i]]$timing['nimble_compile'], 1))
          }
        }
    )
)




