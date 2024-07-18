## SAoptX.R: SA algorithm for constructing LHD/OAD using Rcpp
## modified from SAopt.R in library(NMOF) with option thresholds.only 
## Author: Hongquan Xu, 8/18/23

SAoptX <- function(OF, algo = list(), ...) {
    algoD <- list(nD = 2000L, ## random steps for computing acc. prob.
                  nT = 10L,   ## number of temperatures
                  nS = 1000L, ## steps per temperatures
                  nI = NULL,
                  initT = NULL,    ## starting temperature
                  finalT = 0,      ## final temperature
                  initProb = 0.4,  ## initial acceptance probability
                  alpha = 0.9,     ## temperature-reduction rate
                  mStep = 1,       ## step multiplier
                  x0 = NULL,       ## initial solution
                  neighbour = NULL,
                  printDetail = TRUE,
                  printBar = interactive(),
                  storeF = TRUE,
                  storeSolutions = FALSE,
                  classify = FALSE,
                  OF.target = NULL,
                  thresholds.only = FALSE)  # 

#    checkList(algo, algoD)
    algoD[names(algo)] <- algo
 
    if (!is.null(algoD$nI))
        algoD$nS <- ceiling(algoD$nI/algoD$nT)

    ## user *must* specify the following
    if (is.null(algoD$neighbour))
        stop("specify a neighbourhood function ", sQuote("algo$neighbour"))
    if (!is.function(algoD$neighbour))
        stop(sQuote("algo$neighbour"), " must be a function")
    if (is.null(algoD$x0))
        stop("specify start solution ", sQuote("algo$x0"))

	x0 <- eval(algoD$x0)
    OF1 <- function(x)
        OF(x, ...)
    N1 <- function(x)
        algoD$neighbour(x, ...)

    T <- algoD$initT
    alpha <- algoD$alpha

 	makeInteger=function(x, label, min=0L) {x} 
 
    nT <- makeInteger(algoD$nT, "algo$nT")
    nS <- makeInteger(algoD$nS, "algo$nS")
    nD <- makeInteger(algoD$nD, "algo$nD")
    niter <- nS * nT

    if (!is.null(algoD$initT)) {
        T <- algoD$initT
    } else {
  
        xc  <- x0
        xcF <- OF1(xc)
        diffF <- numeric(nD)
        diffF[] <- NA
        for (i in seq_len(nD)){
            xn  <- N1(xc)
            xnF <- OF1(xn)
            diffF[i] <- abs(xcF - xnF)
            xc  <- xn
            xcF <- xnF
        }
        if (any(is.na(diffF)))
            stop("objective function evaluated to NA")

        T <- try(uniroot(function(T)
                         algoD$initProb - sum(exp(-diffF/T))/nD,
                         interval = c(0.00001, 2))$root, silent = TRUE)
        if (inherits(T, "try-error"))
            T <- -mean(diffF)/log(algoD$initProb)
 
     }

    if (algoD$thresholds.only) { # used by 
    	vT = rep(1, nT)
    	for(i in 2:nT)	vT[i] = vT[i-1] * alpha
        ans <- list(xbest = NA,
                    OFvalue = NA,
                    Fmat = NA,
                    xlist = NA,
                    vT = vT,
                    x0 = x0)
        return(ans)
    }

    ## evaluate initial solution
    xc <- x0
    xcF <- OF1(xc)
    xbest <- xc
    xbestF <- xcF

 
    ## counter = total number of iterations
    counter <- 0L
    for (t in seq_len(nT)) {
        for (s in seq_len(nS)) {
            counter <- counter + 1L

            xn <- N1(xc)
            xnF <- OF1(xn)
            if (xnF <= xcF || exp((xcF - xnF)/T) > runif(1L)) {
                xc <- xn
                xcF <- xnF
                if (xnF <= xbestF) {
                    xbest <- xn
                    xbestF <- xnF
                }
            }

         }

        nS <- round(algoD$mStep*nS)
        T <- T*alpha
    } # for t
 
 
    ans <- list(xbest = xbest, OFvalue = xbestF)
    ans
}
