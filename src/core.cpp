/*
 * core.cpp : core functions called by DE, PSO, GA for constructing LHDs and order-of-addition designs.

Export function: DEx_Rcpp() and TAx_Rcpp() called by MetaSWX.R based on type
  	type <= 0 , a PointNum * FactorNum LHD; 0=maximin, -1=BID, -10=UPD
  	type > 0, a FactorNum * PosNum OAD; 1=maximin, 11=minimax, 21-23=A1opt-A3opt, 25=D-PWO, 
  	default:   // see objFn.Rcpp(x, type, ...) in R

The algorithms return the best design after itermax iterations or maxTime=1800 seconds.

Reference: 
 	Stokes, Z., Wong, W.-K. and Xu, H. (2023). Metaheuristic Solutions for Order-of-Addition Design Problems. Journal of Computational and Graphical Statistics.

Notes:
 	Each agent is a design, stored as a vector of PointNum*FactorNum.
 	Each factor is a permutation of 1 to PointNum
 	PosNum: number of positions (q), which is <= PointNum, if type > 0; 
 		  or number of levels (q) if type <=0.
 	PValue=15: used to compute the phip criterion 
 			or beta used to compute BID value (type=-1)
	pCR=0.5 : crossover probability 
 	pMut=0.1: mutation probability (for crossover factors)
 // For DE strategies, see myDE() in MetaSWX.R
 	pGBest = 1.0 or 0.5: donor selection probability of using GBest as a donor
  	pSelf = 0 or 0.5: donor selection probability of using self as a donor
	(pRandom = 1-pGBest-pSelf: donor selection probability of using a random donor)

 	seed: random seed 
  	vT: a vector of thresholds or temperatures used by TA or SA.
 
* Some codes are borrowed from LaPSO.cpp used in Chen, Hsieh, Hung and Wang (2013)

 Author: Hongquan Xu
 Date: 9/18/2023, add type=-1 for Bayesian-inspired distance criterion
  10/18/2023, add type=-10 for uniform projection criterion/design
  10/23/23, add q=NumPos for BID, type=-1
 */

#include <Rcpp.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace Rcpp;

struct Sort_ref
{
	double fval;
	int    index;
};


// common functions shared by PSO and DE 
void InitializeSwarm(int, int, int, int, int *, struct Sort_ref *);
void InitializeGBest(int, int, int, int, double *, int *, int *, double *);

void UpdatePBest(int, int, int *, double *, int *, double *);
void UpdateGBest(int, int, int *, double *, int *, double *);

// common functions used by PSO, GA and DE
NumericMatrix getFull_R(int, int);
int* Matrix2Vector(NumericMatrix, bool);
NumericMatrix Convert2Design(int, int, int, double, int *, int);

// key function to compute FVal based on type
void ComputeFVal_type(int, int, int, int, int, double, int *, double *, int, int*, int);

// other common functions from LaPSO.cpp
int  MaxMinDist(int, int, int *);
int   compare(const void *, const void *);
void  randperm(int, struct Sort_ref *);
int   minValIdx(int, double *, double *);
int   find(int *, int, int);
int CheckPosNum(int PointNum, int PosNum, int type);
void Swap2Elements(int PosNum,int PointNum, int *Particle);

// functions from yuhao yin

// [[Rcpp::export()]]
double bid_log(NumericMatrix D, double beta)
{ // same as gmp_log 
  int n=D.rows(), m=D.cols();
  int q=max(D)-min(D)+1;		// number of levels

	if(beta==0 && q<n)	return(log(0.0));		// not properly defined
	  
  double obj = 0;
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      double prod = 1.0;
      for(int k=0; k<m; k++){
      	double dijk = (double) (D(i,k) - D(j,k))/q; // scale to [0,1]
  		if(beta==0 && dijk==0){	 // happens if q < n
  			prod *= 1;	// treat beta=0 as epsilon
 // 			print(D);	
 // 			Rcout << "beta=0 with duplicated numbers at i=" << i+1 << "j="<< j+1 << "k=" << k+1 <<"\n";
  		}
        else prod *= 1/(beta + pow(dijk, 2.0));
      }
      obj += prod;
    }
  }
  if(beta != 0)	obj *= pow(beta, m);		// BID
//  printf("beta=%lf, obj=%lf\n", beta, obj);
  
  return (1.0/m)*log(2.0/(n*(n-1))*obj);
}

// [[Rcpp::export()]]
double bid(NumericMatrix D, double beta)
{ // if beta=0, bid is the maxpro criterion
  return exp( bid_log(D, beta) );
}

void ComputeFVal_BID(int PointNum, int FactorNum, int PosNum, int DimensionNum, int ParticleNum, double beta, int *Swarm, double *fSwarm)
{	// q=PosNum is the number of levels, 10/23/23
	int n = PointNum, m = FactorNum, q=PosNum; 
	for( int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
	{
		double sum_phip = 0.0;
		for(int i = 0; i < n-1; i++ )
			for(int j = i+1; j < n; j++ )
			{
				double prod = 1.0;
				for(int k = 0; k < m; k++ )
				{
					int offset = ptcIdx*DimensionNum + k*n;
					double xik   = *(Swarm + offset + i);  
					double xjk   = *(Swarm + offset + j);
					double dijk = (xik-xjk)/q ; // scale to [0,1] 
  					if(beta==0 && dijk==0)	prod *= 1;	// happens if q < n	 
					else prod *= 1/( beta + pow( dijk, 2) );
				}

				sum_phip += prod;
			}
  		if(beta != 0)	sum_phip *= pow(beta, m);		// BID

		fSwarm[ptcIdx] = (1.0/m)*log(2.0/(n*(n-1))* sum_phip); 
	}
}


// -----------------------------------------
// common functions used by PSO, GA and DE
// -----------------------------------------
void Swap2Elements(int PosNum,int PointNum, int *Particle)
{ // Swap 2 elements in Particle; same as SwapQ() in MetaX.R
 // at least one random number is from the first q positions
 
	int randNum = rand() % PosNum;		// PointNum -> PosNum
	int randIdx = rand() % PointNum;
	while( randIdx == randNum ) 	randIdx = rand() % PointNum;

	int temp = *(Particle + randNum);
	*(Particle + randNum) = *(Particle + randIdx);
	*(Particle + randIdx) = temp;	
}

NumericMatrix getFull_R(int FactorNum, int PosNum)
{ // same as R function getFull(FactorNum, PosNum) 
	Environment pkg = Environment::namespace_env("gtools");
    // Picking up permutations() function from the package
    Function getFull = pkg["permutations"];		
    NumericMatrix perms = getFull(FactorNum, PosNum); // same as R
    
	return(perms);
}

int* Matrix2Vector(NumericMatrix mat, bool byrow=true)
{ // convert a R matrix to a C vector for quick computation 
    
    int * xvec = (int*) calloc(mat.nrow() * mat.ncol(), sizeof(int));
    int * x = xvec;		// temp
    
    if(byrow){
 	 	for(int i=0; i<mat.nrow(); i++)
 	 		for(int j=0; j<mat.ncol(); j++)	
 	 			*x++ = mat(i, j);
 	 }
 	 else{
 	 	for(int j=0; j<mat.ncol(); j++)
  	 		for(int i=0; i<mat.nrow(); i++)
 	 			*x++ = mat(i, j);
	 	}  
 	
	return(xvec);
}


NumericMatrix Convert2Design(int PointNum, int FactorNum, int PosNum, double PValue, int *Particle, int type)
{ 	// convert Particle to an LHD if type <=0 or an OAD otherwise
	// PValue is not used
	if(type <= 0){ // LHD
		NumericMatrix design(PointNum, FactorNum);	// R matrix
		for( int ftIdx = 0; ftIdx < FactorNum; ftIdx++ )
		{
			for( int ptIdx = 0; ptIdx < PointNum; ptIdx++ ) 
				design(ptIdx, ftIdx) =  *(Particle + ftIdx*PointNum + ptIdx);
		}
		return(design);  // a PointNum*FactorNum LHD
	}

//	else if(type > 0){ // OAD
		NumericMatrix design(FactorNum, PosNum);	 // PosNum columns
		for( int ftIdx = 0; ftIdx < FactorNum; ftIdx++ )
		{
			for( int ptIdx = 0; ptIdx < PosNum; ptIdx++ ) // PosNum 
				design(ftIdx, ptIdx) =  *(Particle + ftIdx*PointNum + ptIdx);
		}
		return(design);  // a FactorNum*PosNum OAD
//	}
}

void FindOrder(int *x, int n, int *y, struct Sort_ref *sortRef)
{ // y=order(x)
	//--- initialize the sorting reference
	for( int i = 0; i < n; i++ )
	{
		sortRef[i].fval = (double) x[i];
		sortRef[i].index = i;
	}

	//--- sorting
	qsort(sortRef, n, sizeof(struct Sort_ref), compare);

	// y[i] is the order of elements of x[i]
	for( int i = 0; i < n; i++ )		y[i] = sortRef[i].index;
}

void PrintVector(int*x, int n)
{
    for(int i=0; i<n; i++) printf(" %d", x[i]);	printf("\n");
}


NumericMatrix PWO_model(int PointNum, int FactorNum, int PosNum,  int *Particle)
{ // assuming PointNum = PosNum, PosNum is not used.	
		
	struct Sort_ref *sortRef= (struct Sort_ref*)calloc(PointNum, sizeof(struct Sort_ref));
	int *y = (int*)calloc(PointNum, sizeof(int));
	
	int p = 1 + PointNum*(PointNum-1)/2;  // 
	NumericMatrix mat(FactorNum, p);	 // model matrix 
	for( int ftIdx = 0; ftIdx < FactorNum; ftIdx++ )		// each run
	{	
			int *x = Particle + ftIdx*PointNum;		// current row
			
	//		PrintVector(x, PointNum);
			FindOrder(x, PointNum, y, sortRef);	
	//		PrintVector(y, PointNum);
			
			mat(ftIdx, 0) = 1;		// intercept	
			// for each pair of distinct positions
			for( int i = 0, ij=1; i < PointNum -1 ; i++ ){
				int yi = *(y+i);
					for(int j=i+1; j< PointNum; j++){
						int yj = *(y+j);
						mat(ftIdx, ij++)	= (yj - yi) > 0 ? 1 : -1;  // sign(yj-yi)
						}				
				}		 
	//		for(int i=0; i<p; i++) printf(" %.0lf", mat(ftIdx, i));	printf("\n");
	
	}
	free(y);
	free(sortRef);
	return(mat);  // a FactorNum*p matrix
}

NumericMatrix innerproduct_matrix(NumericMatrix X)
{  // information matrix X'X
	int p = X.ncol(), n= X.nrow();
	NumericMatrix XX(p, p);	
  
    for (int i=0; i<p; i++){
          for(int j=0; j<p; j++){
          	int sum = 0;
            for (int k=0; k<n; k++)		sum += X(k, i) * X(k, j) ;
            XX(j, i) = sum; 
            XX(i, j) = sum;
            }
       }
   	return XX;	// p*p matrix
}

void ComputeFVal(int PointNum, int FactorNum, int DimensionNum, int ParticleNum, double PValue, int *Swarm, double *fSwarm)
{
	for( int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
	{
		double sum_phip = 0.0;
		for(int ptIdx1 = 0; ptIdx1 < PointNum-1; ptIdx1++ )
			for(int ptIdx2 = ptIdx1+1; ptIdx2 < PointNum; ptIdx2++ )
			{
				int sum_dist = 0;
				for(int ftIdx = 0; ftIdx < FactorNum; ftIdx++ )
				{
					int offset = ptcIdx*DimensionNum + ftIdx*PointNum;
					int ele1   = *(Swarm + offset + ptIdx1);
					int ele2   = *(Swarm + offset + ptIdx2);

					sum_dist += (ele1 - ele2)*(ele1 - ele2);
				}

				sum_phip += pow( (sqrt(sum_dist)), -PValue);
			}

		fSwarm[ptcIdx] = pow(sum_phip,  1/PValue);
	}

	return;
}

void ComputeFVal_L1(int PointNum, int FactorNum, int DimensionNum, int ParticleNum, double PValue, int *Swarm, double *fSwarm)
{
	for( int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
	{
		double sum_phip = 0.0;
		for(int ptIdx1 = 0; ptIdx1 < PointNum-1; ptIdx1++ )
			for(int ptIdx2 = ptIdx1+1; ptIdx2 < PointNum; ptIdx2++ )
			{
				int sum_dist = 0;
				for(int ftIdx = 0; ftIdx < FactorNum; ftIdx++ )
				{
					int offset = ptcIdx*DimensionNum + ftIdx*PointNum;
					int ele1   = *(Swarm + offset + ptIdx1);
					int ele2   = *(Swarm + offset + ptIdx2);

					sum_dist += abs(ele1 - ele2);
				}

				sum_phip += pow( ((sum_dist)), -PValue);
			}

		fSwarm[ptcIdx] = pow(sum_phip,  1/PValue);
	}

	return;
}

void ComputeFVal_UPD(int PointNum, int FactorNum, int PosNum, int DimensionNum, int ParticleNum, double PValue, int *Swarm, double *fSwarm)
{ // UPD 10/18/23; differ from pcd in Meta-SWX.R by some constant if s<n
 // implement (5.2) in Sun, Wang and Xu (2019, Theorem 2)
 
 	int 	n = PointNum, m= FactorNum;
	int s = PosNum;  // # of levels
	double Cms=(4.0*(5*m-2)*pow(s,4)+30*(3*m-5)*pow(s,2)+15*m+33)/(720*(m-1)*pow(s,4));
	if(s % 2 == 0) Cms += 2.0/(64*pow(s,4));

	for( int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
	{
		double sum_d1ij2 = 0.0, sum_d1i = 0.0;
		for(int i = 0; i < PointNum; i++ ){
			double d1xi = 0.0;
			for(int j = 0; j < PointNum; j++ ){
				int d1ij = 0;	// L1-distance d1(xi,xj)
				for(int k = 0; k < FactorNum; k++ ){
					int offset = ptcIdx*DimensionNum + k*PointNum;
					int ele1   = *(Swarm + offset + i);
					int ele2   = *(Swarm + offset + j);
					d1ij += abs(ele1 - ele2);	// L1-distance
				}	// for k
				d1xi += d1ij;
				sum_d1ij2 += d1ij * d1ij;
			} // for j
			sum_d1i += d1xi * d1xi;
		}	// for i
//		printf("sum_d1ij2=%lf, sum_d1i=%lf\n", sum_d1ij2,  sum_d1i);
		fSwarm[ptcIdx] = (sum_d1ij2 - 2.0 * sum_d1i/PointNum)/(4*m*(m-1)*pow(n*s,2)) + Cms; 
	}

	return;
}

void ComputeFVal_OAD(int PointNum, int FactorNum, int PosNum, int DimensionNum, int ParticleNum, double PValue, int *Swarm, double *fSwarm)
{
	for(int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
	{
		double sum_phip = 0.0;
		int offset0 = ptcIdx*DimensionNum;
		// compute phip for each pair of factors
		for(int ftIdx = 0; ftIdx < FactorNum-1; ftIdx++ ){
			int offset1 = offset0 + ftIdx*PointNum;
			for(int ftIdx2 = ftIdx+1; ftIdx2 < FactorNum; ftIdx2++ )
			{
				int offset2 = offset0 + ftIdx2*PointNum;
				int sum_dist = 0;
				for(int ptIdx1 = 0; ptIdx1 < PosNum; ptIdx1++ )  // first q positions only
				{
					int ele1   = *(Swarm + offset1 + ptIdx1);
					int ele2   = *(Swarm + offset2 + ptIdx1);

					sum_dist += (ele1 - ele2)*(ele1 - ele2);
				}

				sum_phip += pow( (sqrt(sum_dist)), -PValue);
			}
		}
		fSwarm[ptcIdx] = pow(sum_phip,  1/PValue);
	}

	return;
}

int L1dist(int* x, int* y, int n)
{
    int i, val=0;
    for(i=0; i<n; i++) val += abs(x[i]-y[i]); 
    return (val);
}

int L2dist(int* x, int* y, int n)
{
    int i, val=0;
    for(i=0; i<n; i++) val += (x[i]-y[i]) * (x[i]-y[i]); 
    return (val);
}


int dist2Particle(int PointNum, int FactorNum, int PosNum, int *Particle, int *xPt)
{ // return the L2distance of xPt to Particle
	int mindist = 9999999;
	
	// each factor is a run of OAD
	for(int ftIdx = 0; ftIdx < FactorNum; ftIdx++ ){  
		int *curPoint =  Particle + ftIdx*PointNum;
		int dist = L2dist(curPoint, xPt, PosNum );
		if(dist < mindist)	mindist = dist;
		}
	return(mindist);	
}

void ComputeFVal_minimaxOAD(int PointNum, int FactorNum, int PosNum, int DimensionNum, int ParticleNum, int *Swarm, double *fSwarm, int*full, int fullsize)
{ // compute the minimax distance+index/1000 to the full design for each particle
//----------
// see MiniMax(CDesign *design) under work021223/plug/criteria.c
// if CDesign *grid = GridDesign(design->Q, design->D);

	for(int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
	{ // for each particle
		int max_dist = 0, idx=1;
		int* Particle = Swarm + ptcIdx*DimensionNum;  // each particle is an OAD

		// for each point in full, find the distance to the Particle
		for(int i=0; i<fullsize; i++){
			int *xPt = full + i * PosNum ;  // full is a fullRuns * PosNum vector
			int dist = dist2Particle(PointNum, FactorNum, PosNum, Particle, xPt);
			
			// update max_dist and idx
			if(dist > max_dist){
				max_dist = dist;
				idx = 1; // reset idx = 1
			}
			else if(dist == max_dist){
				++idx;  //
			}
		}	// for i
		
		// max_dist is the max distance of points in full to the particle
		fSwarm[ptcIdx] = (double) max_dist + (double) (idx)/ 1000.0;  // L2dist + idx/1000
		if(idx >= 1000)	printf("Warning: max_dist=%d and idx=%d\n", max_dist, idx);
		fSwarm[ptcIdx] = sqrt(fSwarm[ptcIdx]);  // L2dist + idx/1000
	}

	return;
}


void ComputeFVal_R(int PointNum, int FactorNum, int PosNum, int DimensionNum, int ParticleNum, double PValue, int *Swarm, double *fSwarm, int type)
{ // call R function objFn.Rcpp to compute FVal, based on type 
// Retrieving the global environment
//  Environment env = Environment::global_env();
//	Function objFn = env["objFn.Rcpp"];
	Function objFn("objFn.Rcpp");  // objective function in R
	
	for(int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
	{
		int *Particle = Swarm + ptcIdx*DimensionNum; 
		NumericMatrix design = Convert2Design(PointNum, FactorNum, PosNum, PValue, Particle, type); 

		double fval = *REAL(objFn(design, type)); // convert SEXP to double
		fSwarm[ptcIdx] = (double) fval;
	}
	
	return;
}

void ComputeFVal_PWO(int PointNum, int FactorNum, int PosNum, int DimensionNum, int ParticleNum, int *Swarm, double *fSwarm)
{ // same as Deff.PWO()
 	Environment pkg = Environment::namespace_env("base");
    Function det = pkg["det"];		

	for(int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
	{
		int *Particle = Swarm + ptcIdx*DimensionNum; 
		NumericMatrix mat = PWO_model(PointNum, FactorNum, PosNum, Particle); 
		NumericMatrix info = innerproduct_matrix(mat);
		int p = info.ncol();
		int m = PointNum;	// number of components
		
		double fval = 0.0;
		fval = *REAL(det(info)); // convert SEXP to double
		if(isnan(fval))	fval = 0;
//		printf("PWO: ptcIdx=%d p=%d m=%d det=%lf\n", ptcIdx, p, m, fval);

		fval = pow(fval, 1.0/p) / FactorNum; // mat.nrow();	// normalized D-value
		fval /= pow(  pow(m+1, m-1)/pow(3, p-1), 1.0/p);	// D-efficiency
		fSwarm[ptcIdx] = - fval;		// minimize -Deff
//		printf("PWO: ptcIdx=%d fval=%lf\n", ptcIdx, fval);
	}

}

void ComputeFVal_type(int PointNum, int FactorNum, int PosNum, int DimensionNum, int ParticleNum, double PValue, int *Swarm, double *fSwarm, int type, int *full, int fullsize)
{ // key function to compute FVal based on type 
	switch(type){
		case 0: // L2maximin LHD
			ComputeFVal(PointNum, FactorNum, DimensionNum, ParticleNum, PValue, Swarm, fSwarm); 
			break;

		case -5: // L1maximin LHD, 12/19/23
			ComputeFVal_L1(PointNum, FactorNum, DimensionNum, ParticleNum, PValue, Swarm, fSwarm); 
			break;
		
		case -1: // BID, 9/18/23,  PValue is beta
			ComputeFVal_BID(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, PValue, Swarm, fSwarm); 
			break;

		case -10: // UPD, 10/18/23,  PValue is not used
			ComputeFVal_UPD(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, PValue, Swarm, fSwarm); 
			break;
	
		case 1:  // maximin OAD
			ComputeFVal_OAD(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, PValue, Swarm, fSwarm); // OAD
			break;
			
		case 11:  // OAD minimax
			ComputeFVal_minimaxOAD(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, Swarm, fSwarm, full, fullsize); // OAD minimax
			break;
			
		case 25:  // Deff.PWO()
			ComputeFVal_PWO(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, Swarm, fSwarm); // Deff.PWO 
			break; 
			
		default:  // using objFn.Rcpp(x, type, ...) in R 
			ComputeFVal_R(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, PValue, Swarm, fSwarm, type); //  in R	
	}
	return; 
}


int MaxMinDist(int PointNum, int FactorNum, int *LHDarr)
{
	int min_dist = 9999999;

	for( int ptIdx1 = 0; ptIdx1 < PointNum-1; ptIdx1++ )
		for( int ptIdx2 = ptIdx1+1; ptIdx2 < PointNum; ptIdx2++ )
		{
			int sum_dist = 0;
			for( int ftIdx = 0; ftIdx < FactorNum; ftIdx++ )
			{
				int ele1 = *(LHDarr + ftIdx*PointNum + ptIdx1);
				int ele2 = *(LHDarr + ftIdx*PointNum + ptIdx2);

				sum_dist += (ele1 - ele2)*(ele1 - ele2);
			}
//			sum_dist += (ptIdx1 - ptIdx2)*(ptIdx1 - ptIdx2);

			if( sum_dist < min_dist )
				min_dist = sum_dist;
		}

	return -min_dist;
}

int compare(const void *A, const void *B)
{
	Sort_ref *a = (Sort_ref *) A;
	Sort_ref *b = (Sort_ref *) B;

	if( a->fval > b->fval ) return  1;
	if( a->fval < b->fval ) return -1;

	return 0;
}

void randperm(int length, Sort_ref *array)
{
	for( int idx = 0; idx < length; idx++ )
	{
		array[idx].fval  = rand();
		array[idx].index = idx + 1;
	}

	qsort(array, length, sizeof(struct Sort_ref), compare);

	return;
}

int  minValIdx(int length, double *arr, double *minVal)
{
	int minIdx;

	 minIdx = 0;
	*minVal = arr[0];

	for( int idx = 1; idx < length; idx++ )
		if( arr[idx] < *minVal )
		{
			 minIdx = idx;
			*minVal = arr[idx];
		}

	return minIdx;
}

int find(int *arr, int length, int value)
{
	int targetIdx = -1;

	for( int idx = 0; idx < length; idx++ )
		if( arr[idx] == value )
		{
			targetIdx = idx;
			break;
		}

	return targetIdx;
}

