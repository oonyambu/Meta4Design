\name{MetaX, cmpMetaX}
\alias{MetaX}
\alias{cmpMetaX}
\title{Metaheuristic Solutions for Experimental Design Problems}
\concept{minimization, experimental design, metaheuristic algorithms}
\description{
Uses metaheuristic algorithms for constructing order-of-addition designs, Latin hypercube designs and multi-level space-filling designs.

Algorithms include Differential Evolution (DE), Particle Swarm Optimization (PSO), Genetic Algorithm(GA), Thresholding Accepting (TA), and Simulated Annealing (SA).
}
\usage{
MetaX(n, m, q=m, NP=2*n, itermax=1000, method="DE1",  type=1, seed=NULL)

cmpMetaX(K=5, nvec, m, q=m, NP.factor=2, itermax=1000, methods=c("DE1", "PSO1", "GA", "SA", "TA"), type=1, file=TRUE, plot=TRUE, mcl=TRUE)

}
\arguments{
  \item{n}{number of runs.}
  \item{m}{number of components if type>0 or number of factors if type<=0.}
  \item{q}{number of positions if type>0 or number of levels if type <=0. q should be less than or equal to m or n.}
  \item{NP}{population size.}
  \item{itermax}{number of iterations.}
  \item{method, methods}{metaheurisitc algorithm. method should be one of the methods= c("DE1","DE2","DE3", "DE4","PSO1", "PSO2", "PSO3","PSO4","GA", "SA", "SA1", "TA", "SRS"). If type<=0, additional methods include c("SAx", "SAx0", "GAx", "GAx0"). GAx0 is equivalent to GA. }
  \item{type}{type of target designs: order of addition designs if type > 0 or Latin hypercube or multi-level designs if type<=0.
  
  Order of addition designs:  type=1 (maximin), 11 (minimax), 21 (A-optimal, first-order model), 22 (A-optimal, quadratic model), 23 (A-optimal, second-order model), 25 (D-optimal, PWO model)
  
  Latin hypercube or multi-level designs: type=0 (L2-maximin),  type=-1 (BID), type=-5 (L1-maximin), type=-10 (UPD)
  
  }
  
  \item{seed}{seed for random number generators.}
 \item{K}{number of replicates.}
  \item{nvec}{a vector of run sizes.}
  \item{NP.factor}{multiplication factor of the population size. If NP.factor < 10, NP=NP.factor * n; otherwise, NP=NP.factor for all n.}
  \item{file}{save the result in a csv file if it is true.}
  \item{plot}{plot the result if it is true.}
  \item{mcl}{use mclapply if it is true.}
}

\value{
  A list or data.frame with the following components:
 \item{opt.value}{The best value found.}
  \item{xbest}{The best design found (for MataX).}
  \item{time}{Time in seconds used by the algorithm.}
  \item{method}{The algorithm used.}
  \item{type}{The type of design constructed.}
 \item{seed}{The seed used (to reproduce the design).}
}

\note{
 version 0.1.2 (9/26/23): include type=-1 and type=0.
 Add new methods:
 SA1 and TA1 (in MetaX.R),
 GAx=GA_YY, GAx0=GA_base (in GAxy.cpp),   GAx0 = GA (in GAx.cpp).
 and SAx0=SA_MM, SAx=SA_XY (in SAxy.cpp).

Note: GAx is not as good as GAx0=GA.
SAx0 is similar to SLHD (type=0), with different T0 and phip.
SA1 is better than SA; while TA1 is similar to TA.
The performances of SAx and SAx0 highly depend on the buget (=NP*itermax).
In general, DE/PSO/GA are better than SA and TA.

9/22/23: If NP.factor >= 10, NP=NP.factor; otherwise, NP=Np.factor * n. 
9/25/23: q is the number of levels if type <=0.
11/4/23: end search if there is no improvement for maxNoImprove=100 consecutive iterations for DE/PSO/GA/TA/SA/SRS.

version 0.1.3 (10/18/23): add uniform projection criterion/design (type =-10).

version 0.1.4 (12/19/23): add maximin L1-distance criterion/design (type =-5).
}


\references{
  Stokes, Z., Wong, W.-K. and Xu, H. (2023). Metaheuristic Solutions for Order-of-Addition Design Problems. \emph{Journal of Computational and Graphical Statistics}.
  
  Yin, Y. and Xu, H. (2023). Bayesian-Inspired Space-Filling Designs. \emph{Manuscript}.
}

\examples{
# search for maximin design (type=1)
	a = MetaX(n=30, m=6, itermax= 100, method="DE1", type=1);
	a$opt # min pairwise distance
	min(dist(a$xbest)) # = a$opt

# search for minimax design (type=11)
	a = MetaX(n=30, m=6, itermax= 100, method="PSO1", type=11);
	a$opt # max distance

# search for D-optimal design (type=25) under the PWO model
	a = MetaX(n=19, m=4, itermax= 100, method="DE4", type=25);
	a$opt # D-efficiency

# Compare DE1, PSO1 and SRS for constructing maximin designs (type=1)
  methods=c("DE1","PSO1", "SRS")
  m=6; nvec=seq(m, len=5, by=m); K=3
  df = cmpMetaX(K=K, nvec=nvec, m=m, itermax= 100, methods=methods, type=1, mcl=FALSE, plot=FALSE); # without using mclapply
  xtabs(time ~ method + n, data=df)/K  # average time
  xtabs(opt.value ~ method + n, data=df)/K # average opt.value

\dontrun{
# Compare 5 DE variants for constructing maximin designs (type=1)
  methods=c("DE1","DE2", "DE3", "DE4", "SRS")
  m=6; nvec=seq(m, len=5, by=m)
  df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods, itermax= 200, type=1)

# Compare 5 PSO variants for constructing maximin designs (type=1)
  methods=c("PSO1", "PSO2", "PSO3", "PSO4", "SRS")
  m=6; nvec=seq(m, len=5, by=m)
  df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods, itermax= 200, type=1)

# compare 5 algorithms for constructing maximin designs (type=1)
  methods=c("DE1", "PSO1", "GA", "TA", "SA")
  m=6; nvec=seq(m, len=5, by=m)
  df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods, itermax= 200, type=1)

## evaluate the effect of beta
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

## multi-level designs
	n=64; m=2; q=8;
	op = par(mfrow=c(2,4))
	## MaxProQQ design
	x0=generate_rbdesign(n,m,q)/(q-1)	# needs to be in [0,1]
	a=MaxPro::MaxProQQ(x0)
	plot(a$Design[,1:2], main="MaxProQQ", sub=bid.eff(a$Design*(q-1), 1e-5))
	## BID
	betas=c(1e-5, 1e-3, 0.003, 0.005, 0.007, 0.01, 0.05)
	for(beta in betas){	# use DE4, DE2, TA1, SA1 or SAx0=SA_MM, instead of DE1
		a=MetaX(n, m, q, NP=50, itermax = 5000, beta = beta, method="DE2", type=-1); 
		plot(a$xbest[,1:2], main=beta, sub=a$opt.value); print(c(beta, a$opt.value))
	}
	par(op)


# compare algorithms for constructing bidLHDs (type=-1)
  methods=c("DE1",  "DE4", "GAx", "GA")   # GA = GAx0
  m=8; nvec=seq(2*m, len=5, by=2*m)
  df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods, itermax= 200, type=-1)

# compare algorithms for constructing bidLHDs (type=-1)
  methods=c("DE4", "SAx", "SAx0",  "SA1") 
  m=7; nvec=seq(2*m, len=5, by=2*m)
  df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods, itermax= 1000, type=-1)

# compare algorithms for constructing bidLHDs (type=-1)
  methods=c("DE1", "DE4", "GA", "TA", "SAx") 
  m=7; nvec=c(2,5,10)*m
  df = cmpMetaX(K=10, nvec=nvec,  m=m, beta=0.1, methods=methods, itermax= 1000, type=-1)
 
# compare algorithms for constructing bidLHDs (type=-1)
 methods=c("GA", "DE4", "SAx", "SAx0", "SA1")[3:5]
m=10; nvec=seq(2*m, len=5, by=2*m)
df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods,  itermax= 1000, type=-1) 

methods=c("DE1", "DE4", "GAx0", "GAx", "SAx", "SAx0")
m=8; nvec=seq(m, len=5, by=m)[2]
df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods, itermax= 1000, type=-1)

# compare algorithms for constructing maximin LHDs (type=0)
  methods=c("DE1", "PSO1", "GA", "SAx", "SAx0")
  m=6; nvec=seq(m, len=5, by=m)
  df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods, itermax= 200, type=0)

methods=c("DE1", "PSO1",  "GA", "SAx0", "SLHD", "TA")
m=10; nvec=seq(2*m, len=5, by=2*m)
df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods, NP=20, itermax= 200, type=0)

methods=c("DE1", "PSO1",  "GA", "SAx0", "SLHD")
m=8; nvec=seq(2*m, len=5, by=2*m)[5]
df = cmpMetaX(K=10, nvec=nvec,  m=m, methods=methods, itermax= 1000, type=0)

## multi-level designs for type <= 0, q=number of levels 
methods=c("DE1", "PSO1",  "GA", "SA1", "TA")
m=8; nvec=seq(2*m, len=5, by=2*m)
df = cmpMetaX(K=10, nvec=nvec,  m=m, q=m, methods=methods, NP=50, itermax= 200, type=0)

}
