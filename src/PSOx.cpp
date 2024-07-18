/*
 * PSOx.cpp : Discrete PSO for constructing LHDs and order-of-addition designs.
 * This extends the original program (LaPSO.cpp by Hsish) used in Chen, Hsieh, Hung and Wang (2013)

Export function: PSOx_Rcpp() called by MetaSWX.R based on type
  	type <= 0 , a PointNum * FactorNum LHD; 0=maximin
   	type > 0, a FactorNum * PosNum OAD; 1=maximin, 11=minimax, 21-23=A1opt-A3opt, 25=D-PWO, 
  	default:  // see objFn.Rcpp(x, type, ...) in R 

The algorithms return the best design after itermax iterations or maxTime=1800 seconds.

References: 
	Chen, R. B., Hsieh, D. N., Hung, Y. and Wang, W. (2013). Optimizing Latin hypercube designs by particle swarm. Statistics and computing, 23, 663-676.
 	Stokes, Z., Wong, W.-K. and Xu, H. (2023). Metaheuristic Solutions for Order-of-Addition Design Problems. Journal of Computational and Graphical Statistics.

Notes:
 	Each agent is a design, stored as a vector of PointNum*FactorNum.
 	Each factor is a permutation of 1-PointNum
	PosNum: number of positions (q), which is <= PointNum, if type > 0; 
 		  or number of levels (q) if type <=0.
 	PValue=15: used to compute the phip criterion 
 			or beta used to compute BID value (type=-1)
 	pMut=0.1 : mutation probability for each factor
 	SwapNumR=1 : number of random swaps
// for PSO strategies, see myPSO() in MetaSWX.R
	pGBest=0.75 : probability to use GBest 
	pPBest=0.25 : probability to use PBest 
	if SameNumG<0 or SameNumP <0, use the entire factor instead of Move2Best.
 
	seed: random seed 
 
 Author: Hongquan Xu
 Date: 9/13/2023 
 */

#include <Rcpp.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// #define _RCPP_  // include all functions in a single file for sourceCpp("PSOx.cpp")

using namespace Rcpp;

struct Sort_ref
{
	double fval;
	int    index;
};

// PSO specific function
void UpdatePosition(int, int, int, int, int, int, int, double, int, double, double, int *, int *, int *, int *);

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
void Swap2Elements(int PosNum,int PointNum, int *Particle);


// [[Rcpp::export()]]
NumericMatrix PSOx_Rcpp(int PointNum, int FactorNum, int PosNum, int ParticleNum, int IterationNum, double PValue, int SameNumP, int SameNumG, double pMut, int SwapNumR,  double pGBest, double pPBest, int type, int seed, int maxTime=1800, int maxNoImprove=100)
{ 
/* new variables: 
	type: 0 = maximin LHD; 1 = maximin OAD; 11 = minimax OAD; 	others = see objFn.Rcpp(x, type, ...) in R
	PosNum (q): number of positions (<= PointNum)
	pGBest : probability to use GBest 
	pPBest : probability to use PBest 
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
	int* ptChk  = (int    *) calloc(PointNum                , sizeof(int   ));
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
	// DPSO LOOP //
	//-----------//
	time_t start_time = time(NULL);
	int noImprove = 0;

	for( int itIdx = 1; itIdx <= IterationNum; itIdx++ )
	{
		//--- update velocity and position
		UpdatePosition(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, SameNumP, SameNumG, pMut, SwapNumR, pGBest, pPBest, Swarm, PBest, GBest, ptChk);

		//--- compute function values
		ComputeFVal_type(PointNum, FactorNum, PosNum, DimensionNum, ParticleNum, PValue, Swarm, fSwarm, type, full, fullsize);

		//--- update personal best
		UpdatePBest(DimensionNum, ParticleNum, Swarm, fSwarm, PBest, fPBest);

		//--- update global best
		int suc = UpdateGBest(DimensionNum, ParticleNum, PBest, fPBest, GBest, &fGBest);
		if(suc == 1)		noImprove = 0;
		else ++noImprove;
		if(noImprove >= maxNoImprove)	break;	// end search after maxNoImprove

//		printf("iter=%ld: fGBest=%lf\n", itIdx, fGBest);	
		
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
	free(ptChk);
	free(fSwarm);
	free(fPBest);
	free(pArr);
	if(full != NULL)	free(full);
	// ---------------
	
	return xBest;
}


void Move2Best(int PosNum,int PointNum, int SameNum,  int *Particle, int *Best,  int *ptChk)
{ // move Particle towards Best with SameNum positions in the first q=PosNum positions
 // This is the Move() function in Chen et al. (2013)
 
 	if(SameNum > PosNum)		SameNum = PosNum;	// ensure SameNum <= PosNum to avoid forever loop
	
			//--- initialize the check-up table
			for(int idx = 0; idx < PointNum; idx++ ) 	ptChk[idx] = 0;

			for(int idx = 0; idx < SameNum; idx++ )
			{ // randIdx is from the first q positions
				int randIdx = rand() % PosNum;		// PointNum -> PosNum
				while( ptChk[randIdx] == 1 )
					randIdx = rand() % PosNum;		// PointNum -> PosNum
				ptChk[randIdx] = 1;

				int moveToValue = *(Best + randIdx);
				int startIdx    = find(Particle, PointNum, moveToValue);

				*(Particle + startIdx) = *(Particle + randIdx);
				*(Particle + randIdx ) = moveToValue;
			}
//		for(int idx = 0; idx < PointNum; idx++ ) 	printf(" %d", Particle[idx]); printf("\n");
}

// [[Rcpp::export()]]
IntegerVector Move2Best_Rcpp(int PosNum,int PointNum, int SameNum,  IntegerVector vParticle, IntegerVector vBest)
{	// used to verify Move2Best()
	int ptChk[100];
	int Particle[100];
	int Best[100];
	for(int i=0; i<PointNum; i++)	Particle[i] = vParticle(i);
	for(int i=0; i<PointNum; i++)	Best[i] = vBest(i);
	Move2Best(PosNum, PointNum, SameNum,  Particle, Best, ptChk)	;
	
	IntegerVector vParticleR(PointNum);	// copy 
	for(int i=0; i<PointNum; i++)	vParticleR(i) = Particle[i];
	return(vParticleR);
// x=sample(7); y=sample(7);  z=Move2Best_Rcpp(4,7,8,x,y); rbind(x,y,z)
}

void UpdatePosition(int PointNum, int FactorNum, int PosNum,int DimensionNum, int ParticleNum, int SameNumP, int SameNumG, double pMut, int SwapNumR, double pGBest, double pPBest, int *Swarm, int *PBest, int *GBest, int *ptChk)
{
	for( int ptcIdx = 0; ptcIdx < ParticleNum; ptcIdx++ )
		for( int ftIdx = 0; ftIdx < FactorNum; ftIdx++ )
		{
			int offset = ptcIdx*DimensionNum + ftIdx*PointNum;
			int *Particle = Swarm + offset;
			double prob = ( (double)rand() )/RAND_MAX;
			
			//-------------------------------------------------------
			// new PSO feature: move to GBest or PBest with probability pGBest or pPBest
			// copy an entire factor from GBest or PBest if SameNumP=-1 or SameNumG=-1
			//-------------------------------------------------------
			if(prob < pGBest){ // use GBest with 50% or more chance
				//-------------------------//
				// move toward global best //
				//-------------------------//		
				if(SameNumG < 0 || SameNumG > PosNum)
					memcpy(Particle, GBest + ftIdx*PointNum, sizeof(int)* PointNum);  
				else 
					Move2Best(PosNum,  PointNum, SameNumG, Particle, GBest + ftIdx*PointNum,  ptChk);
			}
			else if(prob < pGBest + pPBest){ // use PBest with 25% chance
				//---------------------------//
				// move toward personal best //
				//---------------------------//
				if(SameNumP < 0 || SameNumP > PosNum)
					memcpy(Particle, PBest + offset, sizeof(int)* PointNum);
				else  
					Move2Best(PosNum,  PointNum, SameNumP, Particle, PBest + offset,  ptChk);
			}
			
			// else // keep the current factor with probability 1-pGBest-pPBest

			//-------------//
			// random swap 2 elements SwapNumR times //
			//-------------//
			prob = ( (double)rand() )/RAND_MAX;
			if( prob < pMut ){
				for(int idx = 0; idx < SwapNumR; idx++ )
					Swap2Elements(PosNum, PointNum, Particle);
				}
		}
	
	return;
}
