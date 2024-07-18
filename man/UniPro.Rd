\name{UniPro, UniProX}
\alias{UniPro}
\alias{UniProX}
\title{Differential Evolution (DE) Algorithm for Experimental Design Problems}
\concept{minimization, experimental design, metaheuristic algorithms}
\description{
Differential Evolution (DE) algorithm for constructing uniform projection designs.
}
\usage{
UniPro(n, m, q=NULL, NP=50, itermax=1000,  pMut = 0.1, pCR=0.5, pGBest=0.75, pSelf=0, type=-10, seed=NULL, maxTime=1800)

UniProX(n, m, q=NULL, NP=50, itermax=1000,  pMut = 0.1, pCR=0.5, pGBest=0.75, pSelf=0, type=-10, seed=NULL, maxTime=60)

}
\arguments{
  \item{n}{number of runs.}
  \item{m}{number of factors.}
  \item{q}{number of levels. q should be less than or equal to n. If q=NULL, Latin hypercube designs will be constructed.}
  \item{NP}{population size. Suggested range: [10, 100].}
  \item{itermax}{number of iterations. Suggested range: [500, 1500].}
  \item{pMut}{mutation probability. Range: [0, 1]. }
  \item{pCR}{crossover probability. Range: [0, 1]. }
  \item{pGBest}{probability of selecting the global best donor. Range: [0, 1]. }
  \item{pSelf}{probability of selecting self. Range: [0, 1-pGBest].  }
  \item{type}{type of target designs. type=-10 for uniform projection design.}
  \item{seed}{seed for random number generators.}
  \item{maxTime}{max time (in seconds) allowed.}
 }

\value{
  A list with the following components:
  \item{opt.value}{The best value found for the criterion. The smaller the better.}
  \item{xbest}{An n*m matrix containing the best design found.}
  \item{time}{Time in seconds used by the algorithm.}
  \item{method}{The algorithm used. method="DE".}
  \item{type}{type=-10 for uniform projection design.}
 \item{seed}{The seed used (to reproduce the design).}
}

\note{
 This is a wrapper of myMetaX with method="DE" and type=-10 (version 0.1.4; 4/13/24).
 
 For Stats 143 projects: 
 
 1. How is the performance of the algorithm related to the hyperparameters (NP, itermax, pMut, pCR, pGBest, pSelf), as well as the design size (n, m, q)?
 
 2. How to tune or set the hyperparameters (NP, itermax, pMut, pCR, pGBest, pSelf) to optimize the performance of the algorithm for any design size (n, m, q)?
 
New in Version 0.1.5; 6/11/24
  
 UniPro stops at maxTime=1800 seconds or after itermax iterations.

 UniProX stops at maxTime=60 seconds or no improvement for itermax iterations.
 
}


\references{
  Stokes, Z., Wong, W.-K. and Xu, H. (2023). Metaheuristic Solutions for Order-of-Addition Design Problems. \emph{Journal of Computational and Graphical Statistics}, \url{https://doi.org/10.1080/10618600.2023.2277878}.
  
  Sun, F., Wang, Y. and Xu, H. (2019). Uniform Projection Designs. \emph{Annals of Statistics}, 47, 641-661. \url{https://doi.org/10.1214/18-AOS1705}.

}

\examples{
# search for uniform projection design with default settings
a = UniPro(n=25, m=3); a
pairs(a$xbest)

# try different hyperparameters
a=UniPro(n=25, m=3, NP=100, itermax=500, pMut = 0.1, pCR=0.9, pGBest=0.05, pSelf=0); 
c(opt.value=a$opt, time=a$time)

a=UniPro(n=25, m=3, NP=10, itermax=1500, pMut = 0.9, pCR=0.1, pGBest=0.95, pSelf=0); 
c(opt.value=a$opt, time=a$time)

# repeat multiple times to see the variation of the algorithm
b=replicate(10, UniPro(n=25, m=3, NP=10, itermax=1500, pMut = 0.9, pCR=0.1, pGBest=0.95, pSelf=0)$opt)
summary(b)
boxplot(b)

}
