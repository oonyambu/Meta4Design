\name{Meta4Design-package}
\alias{Meta4Design-package}
\alias{Meta4Design}
\docType{package}
\title{Metaheuristic Solutions for Experimental Design Problems}
\description{
Uses metaheuristic algorithms for constructing order-of-addition designs, Latin hypercube designs and multi-level space-filling designs.

Algorithms include Differential Evolution (DE), Particle Swarm Optimization (PSO), Genetic Algorithm(GA), Thresholding Accepting (TA), and Simulated Annealing (SA).
}
\details{See references and examples.
 }
 
\author{
Hongquan Xu

Maintainer: Hongquan Xu
}
\references{
  Stokes, Z., Wong, W.-K. and Xu, H. (2023). Metaheuristic Solutions for Order-of-Addition Design Problems. \emph{Journal of Computational and Graphical Statistics}.

  Yin, Y. and Xu, H. (2023). Bayesian-Inspired Space-Filling Designs. \emph{Manuscript}.

}

\note{
  Please cite the package and the reference in publications.   Use \code{citation("Meta4Design")}.
}

\keyword{ package }
\seealso{
  \code{\link{MetaX}}, \code{\link{cmpMetaX}}.
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
  
}
