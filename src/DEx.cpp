/*
 * DEx.cpp : DE, TA and SA for constructing LHDs and order-of-addition designs.

Export function: DEx_Rcpp() and TAx_Rcpp() called by MetaSWX.R based on type
  	type <= 0 , a PointNum * FactorNum LHD; 0=maximin, -1=BidLHD
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
 Date: 9/13/2023 
 9/25/23: treat PosNum (q) as number of levels if type <= 0 
 11/4/23: end search if there is no improvement for maxNoImprove=100 consecutive iterations 
 */

#include <Rcpp.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// #define _RCPP_  
// need to include all functions in a single file for sourceCpp("DEx.cpp")
// see common functions in core.cpp

using namespace Rcpp;

struct Sort_ref
{
	double fval;
	int    index;
};

// DE specific function
void UpdatePosition_DE(int, int, int, int, int, double, double, double, double, int *, int *);

// TA, SA specific functions
void Mutate_1Factor(int, int, int, int, int,  int *); // TA and SA
void UpdatePBest_SA(int, int, int *, double *, int *, double *, double);
void UpdatePBest_TA(int, int, int *, double *, int *, double *, double);

// common functions shared by PSO and DE 
void InitializeGBest(int, int, int, int, double *, int *, int *, double *);

void UpdatePBest(int, int, int *, double *, int *, double *);
int UpdateGBest(int, int, int *, double *, int *, double *);

// common functions used by PSO, GA and DE
int CheckPosNum(int PointNum, int PosNum, int type);
void InitializeSwarm(int, int, int, int, int, int *, struct Sort_ref *, int);
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

int CheckPosNum(int PointNum, int PosNum, int type)
{
	if(PosNum > PointNum || PosNum < 0){
		printf("Error: PosNum=%d, but PointNum=%d. Reset PosNum=%d\n", PosNum, PointNum, PointNum);
		PosNum = PointNum;
		}
	if(type > 20) PosNum = PointNum;  //  A- or D-optimality

	return PosNum;
}

// [[Rcpp::export()]]
NumericMatrix DEx_Rcpp(int PointNum, int FactorNum, int PosNum, int ParticleNum, int IterationNum, double PValue,  double pMut,  double pCR, double pGBest, double pSelf, int type, int seed, int maxTime=1800, int maxNoImprove=100)
{ 
/* new variables: 
	type: 0 = maximin LHD; 1 = maximin OAD; 11 = minimax OAD; 	others = see objFn.Rcpp(x, type, ...) in R
	PosNum (q): number of positions (<= PointNum)
	pMut=0.1 : probability to mutate
	pCR=0.5 : probability to crossover (note: at least one crossover factor for DE)
	pGBest=0.5 : probability to use GBest as a donor
	pSelf=0.25 : probability to use self as a donor
	pRandom = 1-pGBest-pSelf: probability to use a random donor
*/

	// reset or check PosNum
	PosNum = CheckPosNum(PointNum, PosNum, type);

	// get full design if type==11 (minimax OAD)
	int *full = NULL, fullsize = 0;
	if(type == 11){
		NumericMatrix mat = getFull_R(PointNum, PosNum);	
		fullsize = mat.nrow();
		full = Matrix2Vector(mat, true);		// alloc memory for full
	}	
	 
	// set seed and DimensionNum
	srand(seed); //	srand(time(NULL));
	int DimensionNum = PointNum * FactorNum; 

	// memory allocation ---------------------------------------------------------
	int* Swarm  = (int    *) calloc(DimensionNum*ParticleNum, sizeof(int   )); 
	int* PBest  = (int    *) calloc(DimensionNum*ParticleNum, sizeof(int   ));
	int* GBest  = (int    *) calloc(DimensionNum            , sizeof(int   ));
	double* fSwarm = (double *) calloc(ParticleNum             , sizeof(double));
	double* fPBest = (double *) calloc(ParticleNum             , sizeof(double));
	struct  Sort_ref *pArr   = (struct Sort_ref*) calloc(PointNum        , sizeof(struct Sort_ref));
	// ---------------------------------------------------------------------------

	//----------------//
	// Initialization //
	//----------------//

	//--- initialize swarm & personal best
	InitializeSwarm(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, Swarm, pArr, type);
	memcpy(PBest, Swarm, sizeof(int)*DimensionNum*ParticleNum);

	//--- compute function values depending type
	ComputeFVal_type(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, PValue, Swarm, fSwarm, type, full, fullsize);
	memcpy(fPBest, fSwarm, sizeof(double)*ParticleNum);

	//--- initialize global best
	double  fGBest;
	InitializeGBest(PointNum, FactorNum, DimensionNum, ParticleNum, fPBest, PBest, GBest, &fGBest);

	//-----------//
	// DE LOOP //
	//-----------//
	time_t start_time = time(NULL);
	int noImprove = 0;

	for( int itIdx = 0; itIdx < IterationNum; itIdx++ )
	{
		//---- DE: set current Swarm as PBest  -----
		memcpy(Swarm, PBest,  sizeof(int)*DimensionNum*ParticleNum);

		//--- DE: update velocity and position
		UpdatePosition_DE(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, pMut, pCR, pGBest, pSelf, Swarm, GBest);

		//--- compute function values
		ComputeFVal_type(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, PValue, Swarm, fSwarm, type, full, fullsize);

		//--- update personal best -- DE selection
		UpdatePBest(DimensionNum, ParticleNum, Swarm, fSwarm, PBest, fPBest);

		//--- update global best
		int suc = UpdateGBest(DimensionNum, ParticleNum, PBest, fPBest, GBest, &fGBest);
		if(suc == 1)		noImprove = 0;
		else ++noImprove;
		if(noImprove >= maxNoImprove)	break;	// end search after maxNoImprove
	
//		printf("iter=%d: fGBest=%lf, type=%d\n", itIdx, fGBest, type);

		time_t cur_time = time(NULL);
		if(difftime(cur_time, start_time) >= maxTime ) 	break; // end search after 30 minutes
	}
 
	//------------------------//
	// convert GBest to xBest:  
	//  LHD: PointNum * FactorNum; OAD: FactorNum *PosNum
	//------------------------//
	NumericMatrix xBest = Convert2Design(PointNum, FactorNum, PosNum, PValue, GBest, type); 

	// free memory ---
	free(Swarm);
	free(PBest);
	free(GBest);
	free(fSwarm);
	free(fPBest);
	free(pArr);
	if(full != NULL)	free(full);
	// ---------------
	
	return xBest;
}

// [[Rcpp::export()]]
NumericMatrix TAx_Rcpp(int PointNum, int FactorNum, int PosNum, int ParticleNum, int IterationNum, double PValue,  NumericVector vT, int SA, int type, int seed, int maxTime=1800, int maxNoImprove=100)
{ 
//  vT is the vector of thresholds from TAoptX in R if SA=0
//  vT is the vector of temperatures if SA=1 (for SA algorithm)
	int nT = vT.length(); // number of thresholds
	int nS = ceil(IterationNum / nT);		// steps per threshold

	// reset or check PosNum
	PosNum = CheckPosNum(PointNum, PosNum, type);

	// get full design if type==11 (minimax OAD)
	int *full = NULL, fullsize = 0;
	if(type == 11){
		NumericMatrix mat = getFull_R(PointNum, PosNum);	
		fullsize = mat.nrow();
		full = Matrix2Vector(mat, true);		// alloc memory for full
	}	
	 
	// set seed and DimensionNum
	srand(seed); //	srand(time(NULL));
	int DimensionNum = PointNum * FactorNum; 

	// memory allocation ---------------------------------------------------------
	int* Swarm  = (int    *) calloc(DimensionNum*ParticleNum, sizeof(int   )); 
	int* PBest  = (int    *) calloc(DimensionNum*ParticleNum, sizeof(int   ));
	int* GBest  = (int    *) calloc(DimensionNum            , sizeof(int   ));
	double* fSwarm = (double *) calloc(ParticleNum             , sizeof(double));
	double* fPBest = (double *) calloc(ParticleNum             , sizeof(double));
	struct  Sort_ref *pArr   = (struct Sort_ref*) calloc(PointNum        , sizeof(struct Sort_ref));
	// ---------------------------------------------------------------------------

	//----------------//
	// Initialization //
	//----------------//

	//--- initialize swarm & personal best
	InitializeSwarm(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, Swarm, pArr, type);
	memcpy(PBest, Swarm, sizeof(int)*DimensionNum*ParticleNum);

	//--- compute function values depending type
	ComputeFVal_type(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, PValue, Swarm, fSwarm, type, full, fullsize);
	memcpy(fPBest, fSwarm, sizeof(double)*ParticleNum);

	//--- initialize global best
	double  fGBest;
	InitializeGBest(PointNum, FactorNum, DimensionNum, ParticleNum, fPBest, PBest, GBest, &fGBest);

	//-----------//
	// TA LOOP //
	//-----------//
	time_t start_time = time(NULL);
	int noImprove = 0;

	for( int itIdx = 0; itIdx < nT; itIdx++ )	// for each threshold or temperature 
	{
		double currThreshold = vT[itIdx];	// current threshold or temperature 
		for(int s=0; s < nS; s++){		// each step
	
			//---- DE: set current Swarm as PBest  -----
			memcpy(Swarm, PBest,  sizeof(int)*DimensionNum*ParticleNum);
	
			//--- neighbor function: swap two elements in a random factor
			Mutate_1Factor(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, Swarm);
	
			//--- compute function values
			ComputeFVal_type(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, PValue, Swarm, fSwarm, type, full, fullsize);
	
			//--- update personal best -- SA or TA selection
			if(SA == 1)	UpdatePBest_SA(DimensionNum, ParticleNum, Swarm, fSwarm, PBest, fPBest, currThreshold);
			else 	UpdatePBest_TA(DimensionNum, ParticleNum, Swarm, fSwarm, PBest, fPBest, currThreshold);
	
			//--- update global best
			int suc = UpdateGBest(DimensionNum, ParticleNum, PBest, fPBest, GBest, &fGBest);
			if(suc == 1)		noImprove = 0;
			else ++noImprove;
			if(noImprove >= maxNoImprove)	break;	// end search after maxNoImprove

			time_t cur_time = time(NULL);
			if(difftime(cur_time, start_time) >= maxTime ) 	break; // end search after 30 minutes
		}
		
	time_t cur_time = time(NULL);
	if(difftime(cur_time, start_time) >= maxTime ) 	break; // end search after 30 minutes

// 	printf("iter=%d with curr Threshold/temperature=%lG: fGBest=%lf\n", itIdx, currThreshold, fGBest);
// 	fflush(stdout);
	}

 
	//------------------------//
	// convert GBest to xBest:  
	//  LHD: PointNum * FactorNum; OAD: FactorNum *PosNum
	//------------------------//
	NumericMatrix xBest = Convert2Design(PointNum, FactorNum, PosNum, PValue, GBest, type); 

	// free memory ---
	free(Swarm);
	free(PBest);
	free(GBest);
	free(fSwarm);
	free(fPBest);
	free(pArr);
	if(full != NULL)	free(full);
	// ---------------
	
	return xBest;
}

void InitializeSwarm(int PointNum, int FactorNum, int PosNum, int DimensionNum, int ParticleNum, int *Swarm, struct Sort_ref *pArr, int type)
{
	int nLevel = PointNum;	// number of levels
	if(type <= 0)	nLevel = PosNum;	
	
	for( int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
		for( int ftIdx = 0; ftIdx < FactorNum; ftIdx++ )
		{
			randperm(PointNum, pArr);

			for( int ptIdx = 0; ptIdx < PointNum; ptIdx++ )
				*(Swarm + ptcIdx*DimensionNum + ftIdx*PointNum + ptIdx) = pArr[ptIdx].index % nLevel +1;		// levels need to be coded as 1:nLevel for A- and D-optimality
		}

	return;
}
	
void InitializeGBest(int PointNum, int FactorNum, int DimensionNum, int ParticleNum, double *fPBest, int *PBest, int *GBest, double *fGBest)
{
	int ptcIdx = minValIdx(ParticleNum, fPBest, fGBest);
	memcpy(GBest, PBest + ptcIdx*DimensionNum, sizeof(int)*DimensionNum);
	return;
}

void Mutate_1Factor(int PointNum, int FactorNum, int PosNum, int DimensionNum, int PopulationNum,  int *Population)
{ // swap two elements in a random column, used by TA/SA

	for( int plIdx = 1; plIdx < PopulationNum; plIdx++ ){
		 // randomly pick a factor to mutate
			int ftIdx = rand() % FactorNum;
			int	offset = plIdx*DimensionNum + ftIdx*PointNum;
			int *Particle = Population + offset;
			Swap2Elements(PosNum, PointNum, Particle);
			}
}

void UpdatePosition_DE(int PointNum, int FactorNum, int PosNum,int DimensionNum, int ParticleNum, double pMut,  double pCR, double pGBest, double pSelf, int *Swarm,  int *GBest)
{

	for( int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
	{
		// pick a random particle other than ptcIdx
		int i0 = rand() % ParticleNum;
		while(i0 == ptcIdx && ParticleNum > 1)		// avoid forever loop when ParticleNum = 1
			i0 = rand() % ParticleNum;
	
		// pick a random factor 
		int j0 = rand() % FactorNum;
 
		for( int ftIdx = 0; ftIdx < FactorNum; ftIdx++ ) // for each factor
		{
			int offset = ptcIdx*DimensionNum + ftIdx*PointNum;
			int offset_i0 = i0*DimensionNum + ftIdx*PointNum;
		
             double prob = ( (double)rand() )/RAND_MAX;
            
            // crossover with probability = pCR (and at least one factor, j0)
           if (ftIdx != j0 && prob > pCR) continue; // no crossover, next factor 
            // This is consistent with the standard DE, keeping the current with 1-pCR.
            // but Pan et al. (2007) always perform the mutation operator (random swap with pMut) if the crossover operator is not performed.
                   
            // DE strategy (selection of donors) based on pGBest and pSelf
            double pRandom = 1- pGBest- pSelf;
			prob = ( (double)rand() )/RAND_MAX;
			if(prob <= pGBest)  // use GBest as a donor with probabilty pGBest=0.5 or 1.0
				memcpy(Swarm+offset, GBest + ftIdx*PointNum, sizeof(int)* PointNum);  
			else if(prob <= pGBest + pRandom) // use random particle (i0) as a donor with  probabilty 1-pGBest-pSelf
				memcpy(Swarm+offset, Swarm + offset_i0, sizeof(int)* PointNum);  
			// else (use self as a donor) with probabilty pSelf
             
			//-------------//
			// random swap //
			//-------------//
			int SwapNumR = 1;  // One swap 
			prob = ( (double)rand() )/RAND_MAX;
			if( prob <= pMut ){
				for(int idx = 0; idx < SwapNumR; idx++ )
					Swap2Elements(PosNum, PointNum, Swarm+offset);
			}
		}
	}
}

void UpdatePBest(int DimensionNum, int ParticleNum, int *Swarm, double *fSwarm, int *PBest, double *fPBest)
{
	for( int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
		if( fSwarm[ptcIdx] < fPBest[ptcIdx] )
		{
			//--- update value
			fPBest[ptcIdx] = fSwarm[ptcIdx];
		
			//--- update LHD
			memcpy(PBest + ptcIdx*DimensionNum, Swarm + ptcIdx*DimensionNum, sizeof(int)*DimensionNum);
		}

	return;
}

void UpdatePBest_SA(int DimensionNum, int ParticleNum, int *Swarm, double *fSwarm, int *PBest, double *fPBest, double currTemp)
{
	for( int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ ){
		// acceptance probability with current temperature
		double prob = ( (double)rand() )/RAND_MAX; 
//	    if (xnF <= xcF || exp((xcF - xnF)/T) > runif(1L))       
		if( fSwarm[ptcIdx] <= fPBest[ptcIdx] || exp((fPBest[ptcIdx] - fSwarm[ptcIdx])/currTemp) > prob )		
		{
			//--- update value
			fPBest[ptcIdx] = fSwarm[ptcIdx];
		
			//--- update LHD
			memcpy(PBest + ptcIdx*DimensionNum, Swarm + ptcIdx*DimensionNum, sizeof(int)*DimensionNum);
		}
	}
	return;
}


void UpdatePBest_TA(int DimensionNum, int ParticleNum, int *Swarm, double *fSwarm, int *PBest, double *fPBest, double currThreshold)
{
	for( int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
		if( fSwarm[ptcIdx] < fPBest[ptcIdx] + currThreshold)
		{
			//--- update value
			fPBest[ptcIdx] = fSwarm[ptcIdx];
		
			//--- update LHD
			memcpy(PBest + ptcIdx*DimensionNum, Swarm + ptcIdx*DimensionNum, sizeof(int)*DimensionNum);
		}

	return;
}

int UpdateGBest(int DimensionNum, int ParticleNum, int *PBest, double *fPBest, int *GBest, double *fGBest)
{
	int    ptcIdx;
	double tempValue;

	ptcIdx = minValIdx(ParticleNum, fPBest, &tempValue);

	if( tempValue < *fGBest )
	{
		//--- update value
		*fGBest = tempValue;

		//--- update LHD
		memcpy(GBest, PBest + ptcIdx*DimensionNum, sizeof(int)*DimensionNum);
		return 1;	// successful
	}

	return 0;	// no improvement
}
