## Meta-SWX.R: meta-heuristic algorithms (DE, PSO, GA, TA and SA) 
## for constructing Order-of-Addition Designs, Latin hypercube designs and multi-level space-filling designs.
##
## References: 
#   Stokes, Z., Wong, W.-K. and Xu, H. (2023). Metaheuristic Solutions for Order-of-Addition Design Problems. Journal of Computational and Graphical Statistics.
#   Yin, Y. and Xu, H. (2023). Bayesian-Inspired Space-Filling Designs. \emph{Manuscript}.
##
## main functions: MetaX(), cmpMetaX()
## see examples.R: test(), run.all()
##
## modified for R package Meta4Design
## version: 0.1.0
## Author: Hongquan Xu (2023/9/15)
## version: 0.1.2 (2023/9/26, add type=-1)
## version: 0.1.3 (2023/10/18, add type=-5)
## version: 0.1.4 (2024/4/13, add UniPro for Stats 143)
## version: 0.1.5 (2024/6/11, add UniProX and early stop)

# History
# 9/15/23: version 0.1.0 finalized for JCGS.
# 9/16/23: version 0.1.1 fix ploting for type=0
# 9/18/23: version 0.1.2 add mySAxy() and type=-1 for bid criterion, with beta=1.0
# add SAxy.cpp and mySAxy() with method="SAx", "SAx1", etc
# 9/18/23, modify TA1 ans SA1 in myTASAx()
# it is better to use a larger nT than default 10
# SA1 is significantly better than SA, TA1 and TA are similar
# 9/20/23: add GAxy.cpp and myGAxy() with method="GAx", "GAx0", "GAx1"
# 9/25/23: use q=NULL in myMetaOpt, myTASAx, etc so that the default q=n for type<=0 and q=m for type >0
# 9/28/23: fix objFn.Rcpp and objFn.out for BID criterion
# 10/17/23: add pMut=0.1 and pCR=0.5 in myMetaX and MetaX
# 10/18/23: add uniform projection criterion/design (type=-10)
# 10/24/23: add bid.lb() and bid(x,beta=0)=dmaxpro(x), 
#	add neighbor.lhd for TASAx when type <= 0 
# 12/19/23: add type=-5, maximin L1-distance 
# 4/13/24: add UniPro for Stats 143 project (with method="DE" and type=-10)
# install.packages("http://www.stat.ucla.edu/~hqxu/pub/Meta4Design_0.1.4.tar.gz", repos = NULL, type="source")
# 7/8/24: change UniPro(pSelf=0) to pSelf=(1-pGBest)/2 to be consistent with HPO-DE.R add_yy()
# 6/10/24: add maxTime=1800, maxNoImprove=1000 in myDE, myPSO, myMetaX() and UniPro(), add UniProX()
# 6/10/24: modify objFn.out() for type==0 and -5 using myphip() instead of min(dist)
#		set UniProX(maxTime=60)

#library(gtools) # permutations
#library(laGP) # distnace
#library(NMOF) # TA and SA algorithms
#library(parallel) # mclapply, works for MacBook Pro; to turn off, use cmp.Kn(mcl=F)

sr=function(dir=".")
{
	library(Meta4Design)
	setwd(dir) # set working directory
	source("Meta-SWX.R") # this file
}

# library(Rcpp)
# sourceCpp("DEx.cpp") # function DEx_Rcpp, TAx_Rcpp
# sourceCpp("PSOx.cpp") # function PSOx_Rcpp
# sourceCpp("GAx.cpp") # function GAx_Rcpp, SRSx_Rcpp

# source("ADopt.R")	# GeomAeff(), GeomAeff_SO(), Deff.PWO()
# source("SAoptX.R")	# SAoptX(), used in myTASAx
# source("TAoptX.R")	# TAoptX(), used in myTASAx

## design criterion and objective function based on type
type2criterion=function(type)
{
	# type=0 (maximin LHD), 1 (maximin OAD), 11 (minimax OAD) are implemented in cpp

	if(type == -10)	criterion = "UPD"
	else if(type == -1)	criterion = "BID"
	else if(type == -5)	criterion = "maximinL1"		# multi-level
	else if(type == 0)	criterion = "maximinL2"		# multi-level
	else if(type < 10)	criterion = "Maximin"
	else if(type >=10 && type < 20)	criterion = "Minimax"
	else if(type == 21)	criterion = "A1opt"	# position models
	else if(type == 22)	criterion = "A2opt"
	else if(type == 23)	criterion = "A3opt"
	else if(type >= 25 && type <= 30)	criterion = "D-optimal"	# PWO model
#	else if(type >= 25 && type <= 30)	criterion = "Dpwo"	# PWO model
	else criterion = "maximin"	# default

	criterion
}

type2ylab=function(type)
{	# called by plot.Kn
	if(type == -10)	ylab = "1000*PCD value" # "Negation of PCD value"
	else if(type == -1)	ylab ="BID Efficiency" # "Negation of BID value"
	else if(type == -5)	ylab = "Minimum L1-Distance"
	else if(type == 0)	ylab = "Minimum L2-Distance"
	else if(type < 10)	ylab = "Minimum Distance"
	else if(type >=10 && type < 20)	ylab = "Negation of Maximum Distance" # 8/29
	else if(type >= 21 && type <= 23)	ylab = "A-Efficiency"
	else if(type >= 25 && type <= 30)	ylab = "D-Efficiency"
	else ylab = "Minimum Distance"	# default

	ylab
}

objFn.Rcpp =function(x, type, full, ...)
{	## objective function called by Rcpp, to be mininized
	## which is useful for self-defined objective function
	# x is an OAD if type >0; or an LHD if type <=0
	# type=0 (maximin LHD), 1 (maximin OAD), 11 (minimax OAD) are implemented in cpp

	if(type == -10) 	 return (pcd(x)) # UPD
	else if(type == -1) 	 return (bid_log(x, ...)) # BID
	else if(type == -5) 	 return (myphip(x, p=1, ...)) # L1-maximin
	else if(type < 10) 	 return (myphip(x, ...)) # L2-maximin
	else if(type >=10 && type < 20) 	return( maxDist(x, full)) 	#
	else if(type == 21) return( -GeomAeff(x, deg=1) )	# A1opt
	else if(type == 22) return( -GeomAeff(x, deg=2) )	# A2opt
	else if(type == 23) return( -GeomAeff(x, deg=3) )	# A3opt: 2nd-order model
	else if(type >= 25 && type <= 30) return( -Deff.PWO(x) )	# PWO model
	else return(myphip(x, ...)) # maximin
}

objFn.out =function(x, type, full, ...)
{	## objective function used for output and display
	# type=0 (maximin LHD), 1 (maximin OAD), 11 (minimax OAD) are implemented in cpp

	if(type == -10) 	 return ( 1000*pcd(x) ) # UPD
	else if(type == -1) 	 return ( bid.eff(x, ...) ) # BID efficiency
#	else if(type == -5)	 return (min(dist(x, method="manhattan"))) # L1maximin
	else if(type == -5) 	 return (myphip(x, p=1, ...)) # L1-maximin, 6/10/24
	else if(type == 0) 	 return (myphip(x, ...)) # L2-maximin, 6/10/24
	else if(type < 10)	 return (min(dist(x))) # L2maximin
	else if(type >=10 && type < 20)	return( maxDist(x, full))
	else if(type >= 21 && type<=30) 	return( -objFn.Rcpp(x, type, full, ...) )	# A- or D-optimality
	else 	return( objFn.Rcpp(x, type, full, ...) )
}

pcd=function(x)
{ # from upd.R 7/11/17
	n=nrow(x); m=ncol(x)
	s=max(x)-min(x)+1	# levels
	if(s%%2) 	Cnms=(4*(5*m-2)*s^4+30*(3*m-5)*s^2+15*m+33)/(720*(m-1)*s^4)
	else Cnms=(8*(5*m-2)*s^4+60*(3*m-5)*s^2+75*m+21)/(1440*(m-1)*s^4)
	x = x - min(x) + 1 # coded as 1,..., s
	x = (x-(s+1)/2)/s  # code for PCD
	d1ij <- dist(x, method = "manhattan", diag=TRUE, upper=TRUE)
	d1ij=as.matrix(d1ij)
	d1i=apply(d1ij,1,sum)
	(sum(d1ij^2) - (2/n)*sum(d1i^2))/(4*m*(m-1)*n^2)+Cnms
}


dmaxpro2=function(x, scale=T)
{ # criterion defined by Joseph et al. (2015), eq. (5)
	n=nrow(x); m=ncol(x); 
	if(scale)	x = x/n		# scale to [0,1] for LHD
	dmp = 0
  for(i in 1:(n-1)) for(j in (i+1):n) {
		dmp = dmp + 1/prod((x[i,]-x[j,])^2)
		}
	return( (dmp/choose(n,2))^(1/m) )
}

bid.lb=function(n,m,s, beta)
{ # lower bound of BID similar to Theorem 5 of Tang and Xu (2021)
	n.k=((s-1):1)*2*n/(n-1)/s^2
	dist.k = (1:(s-1))/s
	if(beta == 0)	beta.k = 1/( 0+dist.k^2 )	# maxpro 
	else beta.k = beta/( beta+dist.k^2 )	# normalized L2-distance
	prod(beta.k^n.k)
}
bid.eff=function(x, beta)
{ # need to update if s<n
	s= max(x)-min(x)+1
	bid.lb(nrow(x), ncol(x), s, beta)/bid(x, beta)
}


randomLHD = function(n, k)
{ # return a random LHD(n, k)
	x = matrix(0, n, k)
	for(i in 1:k) x[,i] = sample(n)
	x
}

getFull = function(m, q=m)
{
#m is the number of components
#q is the number of positions
	gtools::permutations(m,q)
}

fullsize=function(m,q) choose(m,q)*prod(1:q)

# myphip: surrogate of Maximin objective function
myphip = function (x, pow=15, p=2)
{  # Lp-norm
  design <- as.matrix(x)
  D = dist(design, method = "minkowski", diag = FALSE, upper = FALSE, p = p)
  D <- D^(-pow)
  fi_p <- (sum(D))^(1/pow)
  if(fi_p==Inf)   fi_p=999999
  return(fi_p)
}

# Minimax Criteria using full design
maxDist=function(X, full)
{ # return the max min.dist
	disXF = laGP::distance(X, full) # Calculate the squared Euclidean distance between pairs of points and return a distance matrix

	disF = apply(disXF, 2, min) # min dist for each point in Full
#	sqrt(max(disF)) # max min.dist
	max_dist = max(disF)
	idx = sum(disF == max_dist)
	if(idx >= 1000)  cat("max_dist", max_dist, "idx=", idx, "\n")
	sqrt( max_dist + (idx/1000) )
}

## DE, PSO and GA functions using cpp codes
myDEx = function (n, k, q=n, NP=50, itermax=1000, pMut=0.1, pCR=0.5, pGBest=1, pSelf=0, beta=1.0,  type=1, seed=NULL, maxTime=1800, maxNoImprove=1000)
{
#  return an n*k LHD (type=0) or k*n OAD (type>0) by calling DEx_Rcpp() in DEx.cpp
#  Note: n should be total # of components, k should be # of runs, q is # of positions
# pGBest=0.5 or 1.0: probability to use GBest as a donor
# pSelf=0 or 1: probability to use self as a donor
# (pRandom=1-pGBest-pSelf: donor selection probability of using a random donor)

	if(is.null(seed)) seed = sample(100000, 1) # set a seed
	PValue=15;
	if(type == -1)	PValue = beta		# for BID
	gbest=DEx_Rcpp(n, k, q, NP, itermax, PValue, pMut, pCR, pGBest, pSelf, type, seed, maxTime, maxNoImprove);
     gbest
}

myPSOx = function (n, k, q=n, NP=50, itermax=1000, SameNumP = 0, SameNumG = n/4, pMut = 1/(k - 1), pGBest=0.75, pPBest=0.25, beta=1.0,  type=1, seed=NULL, maxTime=1800, maxNoImprove=1000)
{ #  return an n*k LHD (type=0) or k*n OAD (type>0) by calling PSOx_Rcpp() in PSOx.cpp
#	Note: n should be total # of components, k should be # of runs, q is # of positions
#  pMut: mutation probability for each column
#  pGBest and pPBest: probability of taking global and personal best position for each factor.

	if(is.null(seed)) seed = sample(100000, 1) # set a seed
	PValue=15; SwapNumR=1;
	if(type == -1)	PValue = beta		# for BID
	gbest=PSOx_Rcpp(n, k, q, NP, itermax, PValue, SameNumP, SameNumG, pMut, SwapNumR, pGBest, pPBest, type, seed, maxTime, maxNoImprove);
     gbest
}


myGAx = function (n, k, q=n, NP=50, itermax=1000, pMut=0.1, beta=1.0, type=1, seed=NULL, maxTime=1800, maxNoImprove=1000)
{
#  return an n*k LHD (type=0) or k*n OAD (type>0) by calling GAx_Rcpp() in GAx.cpp
# NP=50, pMut=0.1 perturbation probability
# Reference: Liefvendahl, M., and Stocki, R. (2006) A study on algorithms for optimization of Latin hypercubes. Journal of Statistical Planning and Inference, 136, 3231-3247.

	if(is.null(seed)) seed = sample(100000, 1) # set a seed
	PValue=15;
	if(type == -1)	PValue = beta		# for BID
	gbest=GAx_Rcpp(n, k, q, NP, itermax, PValue, pMut, type, seed, maxTime, maxNoImprove);
    gbest
}


mySRSx = function (n, k, q=n, NP=50, itermax=1000, beta=1.0, type=1, seed=NULL, maxTime=1800, maxNoImprove=1000)
{
#  return an n*k LHD (type=0) or k*n OAD (type>0) by calling SRSx_Rcpp() in GAx.cpp

	if(is.null(seed)) seed = sample(100000, 1) # set a seed
	PValue=15;
	if(type == -1)	PValue = beta		# for BID
	gbest=SRSx_Rcpp(n, k, q, NP, itermax, PValue, type, seed, maxTime, maxNoImprove);
    gbest
}

myDE = function(m, n, q=m, NP=50, itermax=1000, method=c("DE1","DE2", "DE3", "DE4")[1],  pMut = 0.1, pCR=0.5, pGBest=0.75, pSelf=0, beta=1.0, type=1, seed=NULL, maxTime=1800, maxNoImprove=1000)
{
#  	DE strategies based on pGBest, pSelf (and pCR, pMut)
	if(method == "DE1") { pGBest=1; pSelf=0 } # use GBest as a donor
	else if(method == "DE2") { pGBest=0; pSelf=1 } # use current/self as a donor
	else if(method == "DE3") { pGBest=0; pSelf=0 } # use a random donor,
	else if(method == "DE4") { pGBest=0.5; pSelf=0.25 } # use Gbest/self/random as a donor, 50/25/25
	else if(method == "DE5") { pGBest=.75; pSelf=0.25 } # use 75/25/0
	else if(method == "DE6") { pGBest=1; pSelf=0; pCR=pCR/4 } #  simialr to DE1 with pCR/4
	else if(method == "DE7") { pGBest=1; pSelf=0; pCR=pCR/2 } # simialr to DE1 with pCR/2
	else if(method == "DE8") { pGBest=1; pSelf=0; pCR=1 }  # simialr to DE1 to minic PSO3 with pCR=1
	else if(method == "DE9") { pGBest=0; pSelf=1; pMut=1; pCR=0 } # to minic TA (with pMut=1 and pCR=0)
	else { pGBest=pGBest; pSelf=pSelf } # default use specific pGBest and pSelf

	xbest <- myDEx(m, n, q=q, NP=NP, itermax =itermax, pMut=pMut, pCR=pCR, pGBest=pGBest, pSelf=pSelf, beta=beta, type=type, seed=seed, maxTime=maxTime, maxNoImprove=maxNoImprove)
	return (xbest)
}

myPSO = function(m, n, q=m, NP=50, itermax=1000, method=c("PSO1", "PSO2","PSO3", "PSO4")[1],  pMut = 0.1, pGBest=0.5, pPBest=0.5, beta=1.0, type=1, seed=NULL, maxTime=1800, maxNoImprove=1000)
{
#  PSO strategies: based on pGbest/pPbest/SameNumG/SameNumP
#		1=best/1 	pGBest=1, pPBest=0,  SameNumG=-1, SameNumP=-1,
# 		2=LaPSO-G   pGBest=1, pPBest=0, SameNumG=m/4
#		3=LaPSO-P   pGBest=0, pPBest=1, SameNumP=m/2
# 		4=hybrid    pGBest=0.5, pPBest=0.5, SameNumG=-1, SameNumP=-1,
#		5-6=hybrid  other combinations 0< pGBest+pPBest <=1
#  SameNumG=-1 has the same effect as SameNumG == m.

	SameNumG=SameNumP=-1 	# default: use the whole variable
	nG=round(m/4); nP=round(m/2); # recommended by Chen et al. (2013)

	if(method == "PSO1"){ pGBest=1; pPBest=0;  }	 # best/1, 100/0/0; cf DE8
	else if(method == "PSO2") { pGBest=1; pPBest=0; SameNumG=nG; SameNumP=0 } # LaPSO-G, 100/0/0
	else if(method == "PSO3") { pGBest=0; pPBest=1; SameNumG=0; SameNumP=nP }	# laPSO-P
	else if(method == "PSO4") { pGBest=.5; pPBest=.5;  } # hybrid, 50/50/0
	else if(method == "PSO5") { pGBest=.75; pPBest=.25;  } # similar to PSO1
	else if(method == "PSO6") { pGBest=.5; pPBest=.5; SameNumG=nG; SameNumP=nP } # similar to PSO3
	else { pGBest=pGBest; pPBest= pPBest;  }  # default use specific pGBest and pPBest

	xbest <- myPSOx(m, n, q=q, NP=NP, itermax=itermax, pMut=pMut,  SameNumG=SameNumG, SameNumP=SameNumP, pGBest=pGBest, pPBest=pPBest, beta=beta, type=type, seed=seed, maxTime=maxTime, maxNoImprove=maxNoImprove)
	return (xbest)
}

myMetaX = function(n, m, q=NULL, NP=50, itermax=1000, method=c("GA", "DE1","DE2", "DE4", "PSO1", "PSO2","PSO4", "SRS")[1],  pMut = 0.1, pCR=0.5, pGBest=0.75, pSelf=0, beta=1.0, type=1, seed=NULL, maxTime=1800, maxNoImprove=1000)
{ #	return an n*q OAD (if type>0) or n*m LHD (if type=0).
	if(is.null(q)){
		if(type > 0)	q = m	# number of positions
		else		q = n		# number of levels 
	}
	if(type > 20)	q=m		# q must be m for A-optimality
	if(type <= 0){	# maximin or bid n*m LHD,  q is the number of levels
		tmp = n; n = m; m = tmp	# switch n and m
	}
		
	DEmethods = c(paste0("DE", 0:9), "DE")
	PSOmethods = c(paste0("PSO", 0:9), "PSO")
	if(is.null(seed)) seed = sample(100000, 1) # set a seed
	ptm <- proc.time()
	if(!is.na(match(method, DEmethods)))
		xbest <- myDE(m, n, q=q, NP=NP, itermax=itermax, method=method, pMut=pMut, pCR=pCR, pGBest=pGBest, pSelf=pSelf, beta=beta,  type=type, seed=seed, maxTime=maxTime, maxNoImprove=maxNoImprove)
	else if(!is.na(match(method, PSOmethods)))
		xbest <- myPSO(m, n, q=q, NP=NP, itermax=itermax, method=method, pMut=pMut, pGBest=pGBest, pPBest=pSelf, beta=beta,  type=type, seed=seed, maxTime=maxTime, maxNoImprove=maxNoImprove)
	else if(method == "GA")
		xbest <- myGAx(m, n, q=q, NP=NP, itermax=itermax, pMut=pMut,  beta=beta, type=type, seed=seed, maxTime=maxTime, maxNoImprove=maxNoImprove)
	else if(method == "SRS")
		xbest <- mySRSx(m, n, q=q, NP=NP, itermax=itermax,  beta=beta, type=type, seed=seed, maxTime=maxTime, maxNoImprove=maxNoImprove)
	else stop(paste(method, "is not implemented in myMetaX."))
	time = proc.time()	- ptm

	# compute bestvalue based on criterion
	full=NULL
	if(type >= 10 && type < 20) full = getFull(m, q) # used for minimax
	if(type <= -1)	bestvalue = objFn.out(xbest, type, full, beta)		# BID, 9/28/23
	else 	bestvalue = objFn.out(xbest, type, full)

	list(opt.value=bestvalue, xbest=xbest, time=time[3], method=method, type=type, seed=seed)
}

UniPro = function(n, m, q=NULL, NP=50, itermax=1000,  pMut = 0.1, pCR=0.5, pGBest=0.75, pSelf=(1-pGBest)/2, beta=1.0, type=-10, seed=NULL, maxTime=1800)
{ # DE algorithm for constructing n*m uniform projection designs with q levels
 # n: run size; m: # of factors; q: # of levels 
 # use q=n or q=NULL to construct Latin hypercube designs
 # The quality of the designs constructed depends on the setting of the hyperparameters.
 # Goal: to find an optimal/robust setting of the hyperparameters
 # Used for Stats 143 project (with method="DE" and type=-10), 4/13/24
 # change pSelf=0 to pSelf=(1-pGBest)/2 to be consistent with HPO-DE.R add_yy(), 6/8/24
 # stop at maxTime=1800 seconds or after itermax iterations, 6/11/24
  	a= myMetaX(n=n, m=m, q=q, NP=NP, itermax= itermax, method="DE", 
  	           pMut=pMut, pCR=pCR, pGBest=pGBest, pSelf=pSelf, 
  	           beta=beta, type=type, seed=seed, maxTime=maxTime, maxNoImprove=itermax)		
	a
}

UniProX = function(n, m, q=NULL, NP=50, itermax=1000,  pMut = 0.1, pCR=0.5, pGBest=0.75, pSelf=(1-pGBest)/2, beta=1.0, type=-10, seed=NULL, maxTime=60)
{ # same as UniPro, except set itermax= 1e6, maxTime=60, maxNoImprove=itermax, 6/11/24
 	# stop at maxTime=60 seconds or no improvement for itermax iterations
  	a= myMetaX(n=n, m=m, q=q, NP=NP, itermax= 1e6, method="DE", 
  	           pMut=pMut, pCR=pCR, pGBest=pGBest, pSelf=pSelf, 
  	           beta=beta, type=type, seed=seed, maxTime=maxTime, maxNoImprove=itermax) 
	a
}


### implement TA and SA with library(NMOF)
## a design is an n*m matrix where each row is a permutation of 1,..., m
##
SwapQ=function(x, q=n)
{ # swap two elements in x, at least one in the first q positions
	n=length(x)
	p = sample(n, 2)
	while(min(p)>q) p = sample(n, 2) # make sure at least one number <=q
	# exchange x[p[1]] and x[p[2]]
	x[p[2:1]] = x[p]
	return(x)
}
## a neighbor() of a design differs in a random row and two random positions
neighbor = function(x, q=m, ...)
{ # assume each row is a permutation of 1-m.
# if q<m, only first q columns will be used
	n=nrow(x); m=ncol(x)
	r = sample(n, 1)  # pick a row
	# exchange two elements in x[r,]
	x[r, ] =SwapQ(x[r,], q)
	return(x)
}

neighbor.lhd = function(x, q=n, ...)
{ # assume each column is a permutation of 1-n
	n=nrow(x); m=ncol(x)
	r = sample(m, 1)  # pick a col
	# exchange two elements in x[,r]
	x[, r] =SwapQ(x[,r], n)
	return(x)
}

## myTASAx using TAx_Rcpp in DEx.cpp
myTASAx=function(n, m, q=NULL, NP=50, itermax=1000, method=c("TA", "SA")[1], beta=1.0, type=1, pMut=NULL, seed=NULL, maxTime=1800, maxNoImprove=1000)
{ #	return an n*q OAD (if type>0) or n*m LHD (if type=0) using TAx_Rcpp in DEx.cpp.
# pMut is not used
	if(is.null(q)){
		if(type > 0)	 q = m	# number of positions
		else		q = n		# number of levels 
	}
	if(type > 20)	q=m		# q must be m for A-optimality
		
	if(is.null(seed)) seed = sample(100000, 1) # set a seed
	set.seed(seed)

	algo <- list() # total #iterations = nT*nS
#	algo$nT <- NP # thresholds or temperatures # use NP instead of nT
#	algo$nS <- itermax # steps
	algo$nI <- itermax # total # of iterations 8/9/23
	algo$neighbour <- neighbor # neighbor function
	algo$printDetail <- FALSE #
	algo$printBar <- FALSE #

	full=NULL
	if(type >= 10 && type < 20) full = getFull(m, q) # used for minimax

	# an initial random n*m design is used to generate thresholds or temperatures only
	algo$thresholds.only = TRUE
	if(type <= 0 ){
		algo$x0 <- generate_rbdesign(n, m, q) # balanced design in SAxy.cpp
		algo$neighbour <- neighbor.lhd # neighbor function for lhd
		OF = function(X, q, full) objFn.Rcpp(X, type, full, beta)  # LHD
		tmp = n; n = m; m = tmp	# switch n and m for TAx_Rcpp
	}
	else{
		algo$x0 <- t(randomLHD(m, n)) # OAD
		# define object_function using the first q columns only
		OF = function(X, q, full) objFn.Rcpp(X[,1:q], type, full)
	}

	# 9/18/23, it is better to use a larger nT than default 10
	# SA1 is significantly better than SA, TA1 and TA are similar
	if(!is.na(match(method, c("TA1",  "SA1")))) { 
		algo$nT <- min(1000, itermax)		
		if(method == "SA1")		algo$alpha = 0.95 	# default 0.9
		}

	ptm = proc.time()
	SA = 0;
	if(!is.na(match(method, c("SA",  "SA1")))){
		SA = 1
#		algo$initProb <- 0.4 ## initial acceptance probability, default 0.4
		a = SAoptX(OF=OF, algo = algo, q=q, full=full)
		}
	else	 		a = TAoptX(OF=OF, algo = algo, q=q, full=full);

#	print(round(a$vT,6))		# examine thresholds or temperatures

	PValue = 15  #  it is very bad to set NP=1 as they convergence to local minima 8/20/23
	if(type == -1)	PValue = beta	# for BID
	xbest = TAx_Rcpp(m, n, q, NP, itermax, PValue, a$vT, SA, type, seed, maxTime, maxNoImprove)
	time = proc.time()	- ptm

	# compute bestvalue based on criterion
	if(type == -1)	bestvalue = objFn.out(xbest, type, full, beta)		# BID, 9/28/23
	else 	bestvalue = objFn.out(xbest, type, full)

	list(opt.value=bestvalue, xbest=xbest, time=time[3],  method=method, type=type, seed=seed)
}

##
## mySAx using MultiSimulatedAnnealing() in SAyy.cpp, 9/18/23
mySAxy=function(n, m, q=n, NP=50, itermax=1000, method=c("SAx", "SAx0", "SAx100")[1], type=-1, beta=1.0,  maxTime=1800, T0=100, Npairs=100,  Imax=100, nFail=3, seed=NULL, pMut=NULL, verbose=0)
{ # mySAxy using SA_MM() and SA_XY() in SAxy.cpp, 9/18/23
 # q is the number of levels
 # Npairs and Imax are tuning parameters to end the search earlier than itermax.
 # pMut is not used
 	if(q < 2 || q > n) stop("q is the number of levels and must be between 2 and n.\n")
	if(type > 0)	stop("type > 0 is not implemented in mySAxy.\n")
	if(type == 0)	beta = 15;		# use p=15 for phip 
	
	if(is.null(seed)) seed = sample(100000, 1) # set a seed
	set.seed(seed)

# all SAx uses up to NP*itermax function evaluations
# SAx0=SA_MM is only good for small n and m.  SAx is better.
# SAx=SAx100=SA_XY is good for small m<=6, but PSO1 and DE1 is better for large m>=8 for BID
# SAx0 is better than SAx given a large budget to allow the algorithm converge. 9/22/23 
# The performances of SAx0 and SAx depend on the budget and how long it takes them to converge.
# Conclusion: DE/PSO is better than SAx0, SAx, and SLHD. 
	ptm = proc.time()
	if(method == "SLHD")	{	# SLHD is similar to SAx0. They differ in T0 and measure/phip
		if(type != 0 || q != n)	stop("SLHD only works for type=0 and q=n")
		a = SLHD::maximinSLHD(t=1, n, m, nstarts=NP, itermax=Imax, total_iter = NP*itermax)
		xbest = as.matrix(a$Design-1); # print(a$measure)
		if(verbose)	 cat("SLHD: T0=", a$temp0, "\n")
	}	
	else if(method == "SAx0"){		# comparable with SAx when Imax=Npairs
		xbest = SA_MM(n, m, q, beta=beta, Ntrials=NP, itermax=itermax, maxTime=maxTime, T0=T0, Imax=Imax, nFail=nFail, type=type, verbose= verbose)
		}
	else{ 	
		no = as.numeric(substring(method,4))
		if(!is.na(no))	Npairs = max(no, 1)		# reset Npairs 
		xbest = SA_XY(n, m, q, beta=beta, Ntrials=NP, itermax=itermax, maxTime=maxTime, T0=T0,  Npairs= Npairs, nFail=nFail, type=type, verbose= verbose)
		}
	time = proc.time()	- ptm

	# compute bestvalue based on criterion
	if(type == -1)	bestvalue = objFn.out(xbest, type, full, beta)		# BID, 9/28/23
	else 	bestvalue = objFn.out(xbest, type, full)

	list(opt.value=bestvalue, xbest=xbest, time=time[3],  method=method, type=type, seed=seed)
}


myGAxy=function(n, m, q=n, NP=50, itermax=1000, method=c("GAx", "GAx0")[1], type=-1, beta=1.0, seed=NULL, maxTime=1800, pElite=0.1, pCR=0.2, pMut=0.1, verbose=0)
{ ## myGAxy using GeneticAlg() in GAxy.cpp, 9/20/23
 # q is the number of levels
 	if(q < 2 || q > n) stop("q is the number of levels and must be between 2 and n.\n")
	if(type > 0)	stop("type > 0 is not implemented in myGAxy.\n")
	if(type == 0)	beta = 15;		# use p=15 for phip 
	
	if(is.null(seed)) seed = sample(100000, 1) # set a seed
	set.seed(seed)

	# new GAx by Yin is not as good as GAx0, which is the same as GA in myGAx
	ptm = proc.time()
	if(method == "GAx0")		## base GA, as Liefvendahl, M. and Stocki, R. (2006).
		xbest = GeneticAlg_base(n, m, q, beta, NP, pMut, itermax, type, verbose)
	else  ## default: "GAx", new GA by Yin using itermax
		xbest = GeneticAlg_YY(n, m, q, beta,  NP, pElite, pCR, pMut, itermax, type, verbose)
	time = proc.time()	- ptm

	# compute bestvalue based on criterion
	if(type == -1)	bestvalue = objFn.out(xbest, type, full, beta)		# BID, 9/28/23
	else 	bestvalue = objFn.out(xbest, type, full)

	list(opt.value=bestvalue, xbest=xbest, time=time[3],  method=method, type=type, seed=seed)
}

test.bid=function()
{ # effects of beta, Fig 1-2

	# Fig 1
	n=80; m=2; q=n;
	op = par(mfrow=c(2,4))
	## MaxPro design
	a=MaxPro::MaxProLHD(n, m)
	plot(a$Design[,1:2], main="MaxProLHD", sub=bid.eff(a$Design*q+0.5, 0))

	## BID
	betas=c(0, 1e-5, 1e-2, 0.05, 0.1, 0.2, 1)
	for(beta in betas){ # DE is better than MaxProLHD when beta=0
		a=MetaX(n, m, q, NP=50, itermax = 5000, beta = beta, method="DE4", type=-1); 
		plot(a$xbest, main=beta, sub=a$opt.value); print(c(beta, a$opt.value))
	}
	par(op)
	
	## Fig 2
	n=64; m=2; q=8;
	op = par(mfrow=c(2,4))
	## MaxProQQ design
	x0=generate_rbdesign(n,m,q)/(q-1)	# needs to be in [0,1]
	a=MaxPro::MaxProQQ(x0)
	plot(a$Design[,1:2], main="MaxProQQ", sub=bid.eff(a$Design*(q-1), 1e-5))
	## BID
	betas=c(1e-5, 1e-3, 0.003, 0.005, 0.007, 0.01, 0.05)
	for(beta in betas){	# use DE4, DE2, TA1, SA1 or SAx0=SA_MM, instead of DE1
		a=MetaX(n, m, q, NP=50, itermax = 5000, beta = beta, method="DE4", type=-1); 
		plot(a$xbest[,1:2], main=beta, sub=a$opt.value); print(c(beta, a$opt.value))
	}
	par(op)
	
	# SAx0 (MM) < SAx100=SAx < SAx10 < PSO1, 9/21/23; SA1 is comparable with SAx10
	methods=c("PSO1", "SAx100", "SAx0", "SAx10", "SA1")#[1:3]
	m=8; nvec=seq(m, len=5, by=2*m)
	df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods, itermax= 1000, type=-1)
	
	## 9/24/23
	# The performances of SAx0 and SAx depend on the budget and how long it takes them to converge.

	a=mySAxy(40,10,method="SAx", verbose=1, NP=40, itermax = 1000, T0=1e2, type=0); a$opt.value
	a=mySAxy(40,10,method="SAx0", verbose=1, NP=40, itermax = 1000, T0=1e2, type=0); a$opt.value
	
	df=cmpMetaX(K=10, n=40, m=10, methods = c( "DE1", "SAx", "SAx0", "PSO1"), NP=10, itermax=1000, type=0)
	df=cmpMetaX(K=10, n=40, m=10, methods = c( "DE1", "SAx", "SAx0", "PSO1"), NP=10, itermax=5000, type=0)
	df=cmpMetaX(K=10, n=40, m=10, methods = c( "DE1", "SAx", "SAx0", "PSO1"), NP=10, itermax=10000, type=0)
	
	# SAx0 is better than SLHD; DE1/PSO1/GA are better than SA and TA.
	methods=c("DE1", "PSO1",  "GA", "SAx0", "SLHD", "TA")
	m=10; nvec=seq(2*m, len=5, by=2*m)
	df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods, NP=20, itermax= 200, type=0)
}

myMetaOpt=MetaX=function(n, m, q=NULL, NP=2*n, itermax=1000, method="DE1", pMut=0.1, pCR=0.5, beta=1.0, type=1, seed=NULL)
{ # method=c("DE1","DE2","DE3", "DE4","PSO1", "PSO2", "PSO3","PSO4","GA", "SA", "TA", "SRS")[1]

	if(is.null(q)){
		if(type > 0)	 q = m	# number of positions
		else		q = n		# number of levels 
	}
	if(is.null(seed)) seed=sample(100000,1)

	if(!is.na(pmatch( c("SAx"), method)) ) myFun=mySAxy # call SA by Yuhao Yin
	else if(!is.na(pmatch( c("SLHD"), method)) ) myFun=mySAxy # call SLHD 
	else if(!is.na(pmatch( c("GAx"), method)) ) myFun=myGAxy # call GA by Yuhao Yin
	else if(!is.na(match(method, c("TA", "SA", "TA1", "SA1")))) myFun=myTASAx # call TAx_Rcpp
	else{ # myFun= myMetaX  # discrete DE, PSO, GA, SRS
		return (myMetaX(n, m, q, NP=NP, method=method, itermax=itermax, pMut=pMut, pCR=pCR, beta=beta, type=type, seed=seed) )
	}

	a=myFun(n, m, q, NP=NP, method=method, itermax=itermax, pMut=pMut, beta=beta, type=type, seed=seed)
	a
}

cmpMeta1=function(n, m, q=NULL, NP=2*n, methods=c("DE1","DE2", "PSO1", "PSO2", "GA", "SA", "TA", "SRS", "SAx"), itermax =1000, pMut=0.1, pCR=0.5, beta=1.0, type=1, seed=NULL)
{ #
	if(is.null(q)){
		if(type > 0)	q = m	# number of positions
		else		q = n		# number of levels 
	}
	if(type > 20)	q=m		# q must be m for A-optimality
	
	res = NULL
	for(i in 1:length(methods) ){
		a=myMetaOpt(n, m, q, NP=NP, method=methods[i], itermax= itermax, pMut=pMut, pCR=pCR, beta=beta, type=type, seed=seed)
		res=rbind(res, c(n, m, q, NP, itermax, a$opt, a$time, a$type, a$seed))
	}
	df = as.data.frame(round(res,5))
	dimnames(df)[[2]]=c("n", "m", "q", "NP", "itermax", "opt.value", "time", "type", "seed")
	df=cbind(method=methods, df)
	df
}


cmpMeta.K1=function(K=2, n, m, q=NULL, NP=2*n, methods=c("DE1", "PSO1", "GA", "SA", "TA", "SRS"), itermax =2000, pMut=0.1, pCR=0.5, beta=1.0, type=1)
{ # repeat cmpMeta1 K times
#	a=replicate(K, cmpMeta1(n, m, q, NP=NP, methods=methods, itermax= itermax, type=type ), simplify=FALSE)
	df=NULL
	for(i in 1:K){
		cat("replicate=", i, date(), "\n")
	#	df=rbind(df, a[[i]])  # replicate
		b = cmpMeta1(n, m, q, NP=NP, methods=methods, itermax= itermax, pMut=pMut, pCR=pCR, beta=beta, type=type)
		df=rbind(df, b)
	}
	df
}

cmpMeta.K=function(K=2, n, m, q=NULL, NP=2*n, methods=c("DE1", "PSO3", "GA", "SA", "TA", "SRS"), itermax =2000, pMut=0.1, pCR=0.5, beta=1.0, type=1, verbose=F)
{ # repeat cmpMeta1 K times with mclapply
 # one common seed is used for all methods for each replicate

	# define a function to be called by mclapply
	fun4mcl=function(seed) cmpMeta1(n, m, q, NP=NP, methods=methods, itermax= itermax, pMut=pMut, pCR=pCR, beta=beta, type=type, seed=seed)

#	require(parallel) # mclapply
	mc.cores = min(parallel::detectCores(), K) # cores. total 12 for 15-inch and 8+2 for 16-inch macbook pro
	if(verbose) cat("cmpMeta.K using mclapply with mc.cores=", mc.cores,  date(), "\n")

	ptm = proc.time()
	seed.K=sample(100000, K)	# use K different seeds for K replicates
	res = parallel::mclapply(seed.K, fun4mcl, mc.cores = mc.cores )
	time = proc.time()	- ptm
	if(verbose)	print(time[1:5])
	df=NULL
	for(i in 1:K){
		df=rbind(df, res[[i]])  # replicate
	}
	df
}

label2=function(nvec)
{
	nn2 = nvec[1]
	if(max(nvec) > min(nvec)) nn2= paste(range(nvec), collapse="-")
	nn2
}

make.filename=function(K, nvec, m, q, NP, itermax, type)
{
	if(is.null(q)){
		if(type > 0)	q = m
		else 	q=nvec
	}
	

	nn2 = label2(nvec)
	m2 = label2(m); q2=label2(q)
	NP2 = label2(NP)
	itermax2 = label2(itermax)
	K2 = label2(K)

	criterion=type2criterion(type)
	if(type <=0)	# LHD, 
		title=paste0(criterion, "_n", nn2, "_m", m2, "_q", q2, "_NP", NP2, "_it", itermax2, "_K", K2)
	else if(type >=21)	#  A- or D-optimality, ignore q
		title=paste0(criterion,  "_m", m2, "_n",nn2,"_NP", NP2, "_it", itermax2, "_K", K2)
	else		# maximin or minimax, include q
		title=paste0(criterion, "_m", m2,"_q",q2, "_n", nn2, "_NP", NP2, "_it", itermax2, "_K", K2)
	title
}

make.title=function(m, q, type)
{ # 9/25/23
	criterion=type2criterion(type)
	if(is.null(q)) paste0(criterion, " (m=",m,")")	
	else if(length(q) == 1 && m == q[1] && type >0) paste0(criterion, " (m=",m,")")	
	else 	 paste0(criterion, " (m=",m,", q=",label2(q), ")")  
}


cmp.Kn = cmpMetaX =function(K=5, nvec, m, q=NULL, NP.factor=2, itermax=1000, methods=c("DE1","PSO1",  "GA", "SA", "TA"), type=1, pMut=0.1, pCR=0.5, beta=1.0,file=FALSE, plot=TRUE, mcl=TRUE)
{ # NP.factor >= 10 as treated as NP
	# create output dir
	criterion = type2criterion(type)
	if(file){
		output.dir = "output"
		if(!dir.exists(output.dir)) dir.create(output.dir)
		output.dir = paste0(output.dir, "/type", type, methods[1])	# 8/21/23
		if(!dir.exists(output.dir)) dir.create(output.dir)
	}

	# set up NP
	if(NP.factor >= 1 && NP.factor < 10){	# changed on 9/22/23
		if(type > 0)		NP = NP.factor * nvec  # OAD
		else 	NP = rep(m*(NP.factor), length(nvec))	# maximin LHD
	}
	else{
		 NP = rep(abs(NP.factor), length(nvec))  # fixed NP
#		 for(i in 1:length(nvec)) NP[i] = max(NP[i], 2*nvec[i]) # NP >= 2n
#		 for(i in 1:length(nvec)) NP[i] = min(NP[i], 10*nvec[i]) # NP <= 10n
		}


	ptm = proc.time()
	df=NULL
	if(mcl){	# use mclapply for each n, each replicate uses the same seed for all methods
		for(i in 1:length(nvec)){
			cat("n=", nvec[i], date(), "\n")
			a =cmpMeta.K(K, nvec[i], m, q, NP=NP[i], itermax=itermax, methods=methods, pMut=pMut, pCR=pCR, beta=beta, type=type);
			df = rbind(df, a)
			}
	}
	else{	# use single cpu, do one replicate at a time
		for(k in 1:K){
			for(i in 1:length(nvec)){
			cat("replicate=", k, "n=", nvec[i], date(), "\n")
			a =cmpMeta1(nvec[i], m, q, NP=NP[i], itermax=itermax, methods=methods, pMut=pMut, pCR=pCR, beta=beta, type=type);
			df = rbind(df, a)
			}
		}
	}
	time = proc.time()	- ptm
	print(time[1:5])

	if(file){
		title = make.filename(K, nvec, m, q, NP, itermax, type)
		filename=paste0(output.dir,"/", title, ".csv")
		write.csv(df,filename)
		cat("===> Results are saved in", filename, "\n")
	}

	if(plot){
		## creat a 2-way table and make a plot
		cat("Average time (in seconds) used:\n")
		a=xtabs(time~n+method, data=df)/K; print(t(a)) ## time used
		cat("Average optimum ", criterion, " value:\n")
		a=xtabs(opt.value ~n+method, data=df)/K; print(t(a))
		plot.Kn(df, fun=mean, type=type)
	}
	df
}


plot.Kn=function(df, fun=mean,  type=1, title=NULL)
{ # df: output from cmp.Kn
	criterion = type2criterion(type)
	if(type >= 10 && type < 20)		df$opt.value = -df$opt.value  # minimax, 8/29
	a=by(df, ~n+method, function(x) fun(x$opt.value))
	df1=array2DF(a, response="opt.value")
	df1$n=as.numeric(df1$n)
	a=xtabs(opt.value ~n+method, data=df1); # print(a)

	# ploting
	nvec = unique(df$n)
	methods= unique(df$method)
	K = nrow(df)/length(nvec)/length(methods)  # replicates
	if(is.null(title))
		title = make.title(unique(df$m), unique(df$q), type)
	ylab = type2ylab(type)

	if(length(nvec)>1){
		lwd = 1.4
		pch=lty=col=1:ncol(a)
	#	ylab=paste("Average", ylab)
		matplot(nvec, a, type="b", xlab="Run Size (n)", ylab=ylab, axes=F, lty=lty, col=col, pch=pch, lwd=lwd, main=title)
		axis(side=1,at=nvec); axis(side=2); box()

		leg.pos = "topright"
		if(type >= 10) leg.pos ="bottomright"  # minimax, A- or D-optimal
#		else if(type ==0) leg.pos ="topleft" # maximinLHD, 9/16/23
		legend(leg.pos, legend=dimnames(a)[[2]], lwd=lwd, lty=lty, col=col, pch=pch)
	}
	else{ # one n
		plot(opt.value~as.factor(method), dat=df, ylab=ylab, xlab="Method", main=title)
		}
	t(a)
}


plot.file1=function(subdir, filename, fun=mean, type=NULL, title=NULL)
{
	df=df0=read.csv(paste0(subdir, "/", filename))

	if(is.null(type)){	# use filename to determine type
		if(!is.na(pmatch("minimax", filename))) type=11
		else if(!is.na(pmatch("Minimax", filename))) type=11
		else if(!is.na(pmatch("maximin", filename))) type=1
		else if(!is.na(pmatch("A1opt", filename))) type=21
		else if(!is.na(pmatch("A2opt", filename))) type=22
		else if(!is.na(pmatch("A3opt", filename))) type=23
		else if(!is.na(pmatch("Dopt", filename))) type=25
		else if(!is.na(pmatch("Dpwo", filename))) type=25
		else type=1 	#  maximin
	}
	a=plot.Kn(df, fun=fun, type=type, title=title)
	a
}

plot.files=function(subdir, fun=mean, type=NULL, title=NULL)
{
	## all files
	files = dir(subdir) #
	for(filename in files ){
		a=plot.file1(subdir, filename, fun=fun, type=type,  title=title)
		cat(paste0("--> file: ", subdir,"/", filename, "\n"))
		print(a)
		}

}
