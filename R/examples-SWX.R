## examples-SWX.R : codes to generate Fig 1-5 in Stokes, Wong and Xu (2023)
## Reference: Stokes, Z., Wong, W.-K. and Xu, H. (2023). Metaheuristic Solutions for Order-of-Addition Design Problems. Journal of Computational and Graphical Statistics.
## 
## examples: see test(), 
## to generate Fig 1-5, use run.all(), which takes about 3h30m
## default using mclapply, to turn it off, set option mcl=F in cmpMetaX()
## 
## Author: Hongquan Xu (2023/9/15) 

sr=function(dir=".")
{
	setwd(dir) # set working directory
	source("Meta-SWX.R") # main codes for MetaX(), cmpMetaX()
	source("examples-SWX.R") # this file
}

test=function()
{ 
## available types and methods
## type=1 (maximin), 11 (minimax), 21 (Aopt, 1st-oder), 22 (Aopt, quadratic), 23 (Aopt, 2nd-order), 25 (Dopt, PWO model)
## methods=c("DE1","DE2","DE3", "DE4","PSO1", "PSO2", "PSO3","PSO4","GA", "SA", "TA", "SRS")
	
	library(Meta4Design)		# 9/15/23
	# search for maximin design (type=1)
	a = MetaX(n=30, m=6, itermax= 1000, method="DE1", type=1);
	a$opt # min pairwise distance	
	# search for minimax design (type=11)
	a = MetaX(n=30, m=6, itermax= 1000, method="PSO1", type=11);
	a$opt # max distance	
	# search for D-optimal design (type=25) under the PWO model
	a = MetaX(n=19, m=4, itermax= 1000, method="DE4", type=25);
	a$opt # D-efficiency	
	
	# compare 5 algorithms for constructing maximin designs (type=1) 
	methods=c("DE1", "PSO1", "GA", "TA", "SA")
	m=6; df = cmpMetaX(K=10, nvec=seq(m, len=5, by=m), m=m, itermax= 200, methods=methods, type=1, file=F, mcl=T); # using mclapply

	m=6; df = cmpMetaX(K=2, nvec=seq(m, len=5, by=m), m=m, itermax= 200, methods=methods, type=1,  file=F, mcl=F); # without using mclapply

}

## codes to general outputs for making plots in the paper
run.fig23=function(K=10, type=1, it=200)
{ 
 	# type = 1 or 11 for maximin or minimax designs
	DE.methods=c("DE1","DE2", "DE3", "DE4", "SRS")
    PSO.methods=c("PSO1", "PSO2", "PSO3", "PSO4", "SRS")

#	type = 1; K=5	 # maximin
#	type = 11; K=5	 # minimax
	m=6; q=m; (nvec=seq(m, len=5, by=m)) # 6-30
	m=7; q=m; (nvec=seq(m, len=5, by=m)) # 7-35, 13min each for minimax
	m=8; q=m; (nvec=seq(m, len=5, by=m)) # 8-40, 1h47m each for minimax
	
	#	par(mfcol=c(2,3)); type = 1; K=5; it=300	
	for(m in 6:8){
	  q=m; (nvec=seq(m, len=5, by=m)) # 
	  df = cmpMetaX(K=K, nvec=nvec,  m=m, q=q, methods=DE.methods, it=it, type=type)
	  df = cmpMetaX(K=K, nvec=nvec,  m=m, q=q, methods=PSO.methods, it=it, type=type)
	}
}

run.fig2=function(K=10,  it=200)
{ # type 1: maximin
	par(mfcol=c(2,3)); 
	run.fig23(K=K, type=1, it=it)
}
run.fig3=function(K=10,  it=200)
{ # type 11: minimax
	par(mfcol=c(2,3)); 
	run.fig23(K=K, type=11, it=it)
}

run.fig4=function(K=10, type=1, it=1000)
{ # compare methods for large run sizes with nvec=seq(2*m, m*10, 2*m)
# maximin: type 1, small difference among TA/GA/DE1/PSO3 for m=4-7
#	type = 1; K=10
#	par(mfrow=c(1,3))
  	methods=c("GA", "SA", "TA", "DE1", "PSO1")
	m=8; q=m; df = cmpMetaX(K=K, nvec=seq(2*m, m*10, 2*m),  m=m, q=q, methods=methods, it=it, type=type) # 
	m=9; q=m; df = cmpMetaX(K=K, nvec=seq(2*m, m*10, 2*m),  m=m, q=q, methods=methods, it=it, type=type) #  (9, 9) or (9, 5)
	m=10; q=m; df = cmpMetaX(K=K, nvec=seq(2*m, m*10, 2*m),  m=m, q=q, methods=methods, it=it, type=type) # (10, 4) or (10, 10)
	
}


run.fig1=function(K=10, it=200)
{ # This is very slow due to mulitple models and using R codes.

	# min run size for 2nd-order model
#	min.N = c(NA, NA, 5, 9, 14, 20, 27, 35, 44, 54) 	# (m-1)*(m+2)/2
    DE.methods=c("DE1","DE2", "DE3", "DE4", "SRS")
    PSO.methods=c("PSO1", "PSO2", "PSO3", "PSO4", "SRS")

	par(mfcol=c(2,3)); 
	# A3opt: type=23, p=(m-1)*(m+2)/2
 	m=6; (nvec=seq((m-1)*(m+2)/2, len=5, by=m)) # 20-44, 50min
	df = cmpMetaX(K=K, nvec=nvec,  m=m,  methods=PSO.methods, it=it, type=23)
	df = cmpMetaX(K=K, nvec=nvec,  m=m, methods=DE.methods, it=it, type=23)
	
	## A2opt: type=22, p=2*m-1
	m=7; (nvec=seq(2*m-1, len=5, by=m)) # 15-43, 9min
	df = cmpMetaX(K=K, nvec=nvec,  m=m, methods=DE.methods, it=it, type=22)
	df = cmpMetaX(K=K, nvec=nvec,  m=m,  methods=PSO.methods, it=it,type=22)

	## A1opt: type=21, p=m
	m=8; (nvec=seq(m, len=5, by=m)) # 8-40, 8min
	df = cmpMetaX(K=K, nvec=nvec,  m=m, methods=DE.methods, it=it, type=21)
	df = cmpMetaX(K=K, nvec=nvec,  m=m,  methods=PSO.methods, it=it,type=21)
	
}

run.fig5=function(K=10, it=1000)
{
	methods=c("GA", "DE1", "PSO1", "SA", "TA")
	m=8; mm= m*(m-1)/2; (nvec=seq(mm+1, len=5, by=2*m));  # 13 min
	m=9; mm= m*(m-1)/2; (nvec=seq(mm+1, len=5, by=2*m));  # 30 min
	m=10; mm= m*(m-1)/2; (nvec=seq(mm+1, len=5, by=2*m)); # 74 min

#	par(mfrow=c(1,3))
	for(m in 8:10){
		mm= m*(m-1)/2; (nvec=seq(mm+1, len=5, by=2*m)); 
		df = cmpMetaX(K=K, nvec=nvec,  m=m,  methods=methods, type=25, it=it, plot=T)
#		plot.Kn(df, max, type=25)	# compare the max
		}
}

## code to generate all fig 1-5.
run.all=function()
{ # use it=200 for fig 1-3 and it=1000 for fig 4-5; 8/30/23
# include fig 1-3 with it=1000 in appendix.
  run.fig2(it=200)	# < 1 min
  run.fig1(it=200)	# about 30m 
  date() # end time
  run.fig3(it=200)	# about 50m 
  date() # end time
 
  par(mfrow=c(2,3))
  run.fig4(it=1000)	 # 6 min for it=1000
  date() # end time
  run.fig5(it=1000) # about 2h
  date() # end time
}

