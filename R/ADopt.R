## ADopt.R : A-optimality for position-based models and D-optimality for PWO models
## Functions:  GeomAeff(x, deg), and Deff.PWO(x) 
## sourced by MetaSWX.R 
## Reference: Stokes, Z., Wong, W.-K. and Xu, H. (2023). Metaheuristic Solutions for Order-of-Addition Design Problems. Journal of Computational and Graphical Statistics.

## Authors: Stokes, Wong and Xu (2023/9/13) 

#library(gtools) # permutations, combinations

## A-optimality of order-of addition designs
Aopt=function(X, Eps=1e-10)
{ # X is the model matrix
# A-optimality, the smaller the better: return tr(solve(t(X)%*%X)), 
	p = ncol(X)
	xx = t(X)%*%X
	a=eigen(xx, sym=T, only=T)
	if(a$value[p] < Eps) return(0)  # singular
	return(nrow(X) * sum(1/a$value[1:p]))	# normalize with run size
}

  
get.interaction=function(x)
{	# interaction of x
	sub = gtools::combinations(ncol(x), 2)	
	xx = apply(sub, 1, function(s) x[, s[1]]*x[, s[2]] ) 
	if(nrow(sub) == 1)	return(matrix(xx, ncol=1))	# as a matrix
	(xx)
}

pos.model=function(x0, deg=1)
{ # x0 is a component matrix where each row is a permutation of 1:m
	m=ncol(x0)
	x = t(apply(x0, 1, order)) 	# change to position matrix
	x = x - (m+1)/2	 # center
	c1 = sqrt(12/(m^2-1)); 	
	x1 = c1 * x
	if(deg == 1) return(x1) # FO model
	
	c2 = sqrt(15/(m^2-4)) 
	x2 = (c2/c1) * (x1^2 -1) # see Mee (2024)
	if(deg == 2) 	return ( cbind(x1, x2) ) # PQ model

	if(deg > 3) return (NULL)
	
	# 2nd-order model
	xx = get.interaction ( x1 )
	return ( cbind(x1, x2, xx) )	 # SO model with all terms
}


## global Aopt for full design, used in computing Aeff
getAopt.full=function(m)
{ # compute just once
	full = gtools::permutations(m, m)
	x = pos.model(full, deg=2)
	AL = Aopt(cbind(1, x[,2:m]))	 # linear
	AQ = Aopt(cbind(1, x[,-c(1,m+1)]))	# remove 2 positions
#	x = pos.model(full, deg=3)
	xx = get.interaction ( x[, 2:m] )
	xx1 =cbind(1, x[,-c(1, m+1, m+2)], xx) # remove an extra quadratic term
	A2 = Aopt(xx1)	
	c(AL, AQ, A2)
}

## pre-calculated Aopt.full for m=3:8 once
Aopt.full  = matrix(NA, 10, 3)
for(m in 3:8) Aopt.full[m, ] = getAopt.full(m)

GeomAeff_SO=function(x0)
{ # x0 is a component matrx, where each row is a permutation of 1:m
 # require Aopt.full[,]
 	#  return geometrical mean of relative A efficiency 
	#  deg=3: SO model excludes one position and one extra quadratuc term 
 	m = ncol(x0)
 	x = pos.model(x0, deg=2)
 	Aval = matrix(1, m, m)		# set Aval[i,i] = 1 so it does not affect prod(Aval)
 	for(i in 1:m){
 		idx = (1:m)[-i] # exclude ith position
 		xx = get.interaction ( x[, idx] )
 		for(j in idx){
 			xx2 = cbind(1, x[, -c(i, m+i, m+j)], xx)		# exclude ith position and jth quadratic effect
 			Aval[i, j] = 	Aopt(xx2)
 		}
 	}		
 	GeomA = prod(Aval)^(1/m/(m-1))  # geometrical mean of A-optimality
 	if(is.na(GeomA))	return(0)	# non-estimable
 	if(GeomA > 0){
 		if(is.na(Aopt.full[m, 3]))	return(1/GeomA)
 		return( Aopt.full[m, 3]/GeomA )	# relative efficiency to full design
 	}	
 	else return (0)		# non-estimable
}

GeomAeff=function(x0, deg=1)
{ # x0 is a component matrx, where each row is a permutation of 1:m
 # require Aopt.full[,]
 	#  return geometrical mean of relative A efficiency 
 	#  deg=1: FO model excludes one position
	#  deg=2: PQ model excludes one position
	#  deg=3: SO model excludes one position and one extra quadratuc term 

	if(deg==3) return ( GeomAeff_SO(x0) ) # 2nd-order model
	
 	m = ncol(x0)
 	x = pos.model(x0, deg=deg)
 	Aval = rep(NA, m)
 	if(deg == 1) for(i in 1:m)		Aval[i] = 	Aopt(cbind(1, x[,-i]))
 	else if(deg == 2) for(i in 1:m)		Aval[i] = 	Aopt(cbind(1, x[,-c(i,m+i)]))
	GeomA = prod(Aval)^(1/m)  # geometrical mean of A-optimality
  	if(is.na(GeomA))	 return(0)	# non-estimable
	if(GeomA > 0){
		if(is.na(Aopt.full[m, deg]))	return(1/GeomA)
		return( Aopt.full[m, deg]/GeomA )	# relative efficiency to full design
	}
	else return (0)		# non-estimable
}

## D-optimality for PWO model
perm.pair2=function(perm, m=q, debug=F, ...)
{ # construct a pair contrast for a permutation, perm
 # include distance	
	q =length(perm) # number of components
	x =matrix(0,m,m)
	for(i in 1:(q-1)) for(j in (i+1):q) # all pairs, extending Lin's pairwise model
	{ # extending Lin's model, using a function of (j-i)
		if(perm[i]<perm[j])  x[perm[i], perm[j]] = (j-i)
		else x[perm[j], perm[i]] = -(j-i)
		}	
	# convert x to a contrast matrix
	if(debug)print(x)
	z=x[row(x)<col(x)]  # upper diagonal
		# add label
	lab=outer(1:m, 1:m, paste, sep=".")
#	lab=t(lab)
	lab=c(lab[row(x)<col(x)])
	lab=paste("I",lab, sep="")
	names(z)=lab
	z
}

perm.pwo1=function(perm, taper=0, ...)
{ # PWO
	d = perm.pair2(perm, ...)	# signed distance between components
	if(taper==1)		1/d		# tapered PWO with c_h = 1/h
	else if(taper==2)	sign(d) * (1/2)^(abs(d)-1)	# tapered PWO with c_h = (1/2)^(h-1)
	else		sign(d)		# original PWO model
}


PWO.model = function(perms,  taper=0, ...)
{ # PWO.model(c(1,2,5), m=5)
	if(is.vector(perms)) perm.pwo1(perms, taper=taper, ...)
	else t(apply(perms, 1, perm.pwo1, taper=taper,...))
}

Deff.PWO=function(x, taper=0)
{ # x is the component matrix
	m = ncol(x)	# number of components
	X = cbind(1, PWO.model(x, taper=taper) )
	p = ncol(X)
	XX = t(X)%*%X
	De = det(XX)^(1/p) / nrow(X)
	if(is.na(De))	return(0)	# singular
	if(taper == 0) 	De / ( (m+1)^(m-1)/3^(p-1) )^(1/p)		# D-efficiency, 
	else De  # need to be adjusted for tapered models
}

