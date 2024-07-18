## TAoptX.R: TA algorithm for constructing LHD/OAD using Rcpp
## simplified from TAopt.R in library(NMOF) 
## Author: Hongquan Xu, 8/18/23

TAoptX <- function(OF, algo = list(), ...) {
    algoD <- list(nD = 2000L, ## random steps for computing thresholds
                  nT = 10L,   ## number of thresholds
                  nS = 1000L, ## steps per threshold
                  nI = NULL,  ## total number of iterations
                  q = 0.5,    ## starting quantile for thresholds
                  x0 = NULL,  ## initial solution
                  vT = NULL,  ## threshold sequence
                  neighbour = NULL,
                  printDetail = TRUE,
                  printBar = interactive(),
                  stepUp = 0L,
                  scale = 1,
                  drop0 = FALSE,
                  storeF = TRUE,
                  storeSolutions = FALSE,
                  classify = FALSE,
                  OF.target = NULL,
                  thresholds.only = FALSE)

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

 	makeInteger=function(x, label, min=0L) {x} 

    nT <- makeInteger(algoD$nT, "algo$nT")
    nS <- makeInteger(algoD$nS, "algo$nS")
    nD <- makeInteger(algoD$nD, "algo$nD")
    stepUp <- makeInteger(algoD$stepUp, "algo$stepUp", 0L)
    niter <- nS * nT * (stepUp+1L)

    ## compute thresholds
    if (is.null(algoD$vT)) {
        if (algoD$q < .Machine$double.eps^0.5) {
            vT <- numeric(nT)
        } else {
#           cat("\n  Computing thresholds ... ")
            xc  <- x0
            xcF <- OF1(xc)
   #         print(xc); print(xcF)
            diffF <- numeric(nD)
            diffF[] <- NA
            for (i in seq_len(nD)){
                xn  <- N1(xc)
                xnF <- OF1(xn)
  #              print(xn); print(xnF)
                diffF[i] <- abs(xcF - xnF)
                xc  <- xn
                xcF <- xnF
            }
            vT <- algoD$q * ((nT - 1L):0)/nT
            if (any(is.na(diffF)))
                stop("objective function evaluated to NA")
            if (algoD$drop0)
                diffF <- diffF[diffF != 0]
            vT <- quantile(diffF, vT, na.rm = FALSE)
            vT[nT] <- 0  ## set last threshold to zero
        }
    } else {
        vT <- algoD$vT
    }
    if (stepUp > 0L)
        vT <- rep.int(vT, stepUp + 1L)
    if (algoD$scale < 0) {
        scale <- 0
        warning(sQuote("scale"), " set to 0 (",
                sQuote("scale"), " must be nonnegative)")
    } else {
        scale <- algoD$scale
    }
    vT <- vT * scale
    nT <- length(vT)
    niter <- nS * nT

    if (algoD$thresholds.only) {
        ans <- list(xbest = NA,
                    OFvalue = NA,
                    Fmat = NA,
                    xlist = NA,
                    vT = vT,
                    x0 = x0)
        if (algoD$classify)
            class(ans) <- "TAopt"
        return(ans)
    }

    ## evaluate initial solution
    xc <- x0
    xcF <- OF1(xc)
    xbest <- xc
    xbestF <- xcF

    counter <- 0L
    for (t in seq_len(nT)) {
        for (s in seq_len(nS)) {
            ## counter = total number of iterations
            counter <- counter + 1L

            xn <- N1(xc)
            xnF <- OF1(xn)
            if (xnF <= xcF + vT[t]) {
                xc <- xn
                xcF <- xnF
                if (xnF <= xbestF) {
                    xbest <- xn
                    xbestF <- xnF
                }
            }
            
		} # for s
   }		# for t


    ans <- list(xbest = xbest, OFvalue = xbestF,
                 vT = vT)
    ans
}

