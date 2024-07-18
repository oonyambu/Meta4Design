/*
 * GAx.cpp : Discrete GA for constructing LHDs and order-of-addition designs.
 * This extends the original program (GA.cpp by Hsish) used in Chen, Hsieh, Hung and Wang (2013)

Export function: GAx_Rcpp() and SRSx_Rcpp() called by MetaSWX.R based on type
  	type <= 0 , a PointNum * FactorNum design; 0=maximin, -1=BID, -10=UPD
  	type > 0, a FactorNum * PosNum OAD; 1=maximin, 11=minimax, 21-23=A1opt-A3opt, 25=D-PWO,
  	default:  // see objFn.Rcpp(x, type, ...) in R

The algorithms return the best design after itermax iterations or maxTime=1800 seconds.

References:
	Liefvendahl, M. and Stocki, R. (2006). A study on algorithms for optimization of Latin hypercubes. Journal of Statistical Planning and Inference, 136, 3231-3247.
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
 	seed: random seed

 Author: Hongquan Xu
 Date: 9/13/2023
 */

#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include "metaSWX.h"

// #define _RCPP_  // include all functions in a single file for sourceCpp("PSOx.cpp")


using namespace Rcpp;

// defined in LaPSOx.cpp
struct Sort_ref
{
	double fval;
	int    index;
};

// GA specific functions
void InitializePopulation_OAD(int, int, int, int, int, int *, int *, int);

void SelectSurvivor(int, int, int, int, int *, double *, int *, int *, double *, struct Sort_ref *);
void CrossOver(int, int, int, int, int *, int *, int *);
void Mutate(int, int, int, int, int, double, int *);

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
NumericMatrix GAx_Rcpp(int PointNum, int FactorNum, int PosNum, int PopulationNum, int IterationNum, int PValue, double pMut,  int type, int seed, int maxTime=1800, int maxNoImprove=100)
{
/* new variables:
	type: 0 = maximin LHD; 1 = maximin OAD; 11 = minimax OAD; 	others = see objFn.Rcpp(x, type, ...) in R
	PosNum (q): number of positions (<= PointNum)
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

	double fBest;
	int    minIdx;

	// memory allocation -----------------------------------------------------------
	int* Population  = (int   *) malloc(sizeof(int  )*DimensionNum* PopulationNum    );
	int* Survivor    = (int   *) malloc(sizeof(int  )*DimensionNum*(PopulationNum>>1));
	int* Best        = (int   *) malloc(sizeof(int  )*DimensionNum);
	int* randIdx     = (int    *) malloc(sizeof(int   )*PopulationNum);
	double *fPopulation = (double *) malloc(sizeof(double)*PopulationNum);
	struct Sort_ref *pArr  = (struct Sort_ref *) malloc(sizeof(struct Sort_ref)*PointNum     );
	struct Sort_ref *sortRef = (struct Sort_ref *) malloc(sizeof(struct Sort_ref)*PopulationNum);

		//--- initialize population
		InitializeSwarm(PointNum, FactorNum, PosNum, DimensionNum, PopulationNum, Population, pArr, type);

		//--- compute function values
		ComputeFVal_type(PointNum, FactorNum, PosNum, DimensionNum, PopulationNum, PValue, Population, fPopulation, type, full, fullsize);

		//---------//
		// GA LOOP //
		//---------//
		time_t start_time = time(NULL);
		double minFval = 99999999.0;
		int noImprove = 0;

		for( int itIdx = 0; itIdx < IterationNum; itIdx++ )
		{
			//--- select survivor
			SelectSurvivor(PointNum, FactorNum, DimensionNum, PopulationNum, Population, fPopulation, Survivor, Best, &fBest, sortRef);
			
			if( itIdx == 0  || fBest < minFval)
			{
				minIdx = minValIdx(PopulationNum, fPopulation, &fBest);
				minFval = fBest; // fBest = fPopulation[minIdx]
				memcpy(Best, Population + minIdx*DimensionNum, sizeof(int)*DimensionNum);
				noImprove = 0;
			}
			else ++noImprove;
			if(noImprove >= maxNoImprove)	break;	// end search after maxNoImprove

			//--- cross-over the survivor
			CrossOver(PointNum, FactorNum, DimensionNum, PopulationNum, Survivor, Population, randIdx);

			//--- mutate the new population (except for the best one)
			Mutate(PointNum, FactorNum, PosNum, DimensionNum, PopulationNum, pMut, Population);

			//--- compute function values
			ComputeFVal_type(PointNum, FactorNum, PosNum, DimensionNum, PopulationNum, PValue, Population, fPopulation, type, full, fullsize);

			time_t cur_time = time(NULL);
			if(difftime(cur_time, start_time) >= maxTime ) 	break; // end search after 30 minutes

		} //--- end of GA loop

		minIdx = minValIdx(PopulationNum, fPopulation, &fBest);
		memcpy(Best, Population + minIdx*DimensionNum, sizeof(int)*DimensionNum);
//		int mmDist = MaxMinDist(PointNum, FactorNum, Best);

	//------------------------//
	// convert Best to xBest:
	//  LHD: PointNum * FactorNum; OAD: FactorNum *PosNum
	//------------------------//
	NumericMatrix xBest = Convert2Design(PointNum, FactorNum, PosNum, PValue, Best, type);

	// memory deallocation ---
	free(Population);
	free(Survivor);
	free(Best);
	free(randIdx);
	free(fPopulation);
	free(pArr);
	free(sortRef);
	if(full != NULL)	free(full);
	// -----------------------

	return xBest;
}

// [[Rcpp::export()]]
NumericMatrix SRSx_Rcpp(int PointNum, int FactorNum, int PosNum, int PopulationNum, int IterationNum, int PValue,  int type, int seed, int maxTime=1800, int maxNoImprove=100)
{
	// reset or check PosNum
	PosNum = CheckPosNum(PointNum, PosNum, type);

	// get full design if type > 0 (OAD)
	int *full = NULL, fullsize = 0;
	if(type > 0){
		NumericMatrix mat = getFull_R(PointNum, PosNum);
		fullsize = mat.nrow();
		full = Matrix2Vector(mat, true);		// alloc memory for full
	}

	// set seed and DimensionNum
	srand(seed); //	srand(time(NULL));
	int DimensionNum = PointNum * FactorNum;

	// memory allocation -----------------------------------------------------------
	int* Population  = (int   *) malloc(sizeof(int  )*DimensionNum* PopulationNum    );
	int* Best        = (int   *) malloc(sizeof(int  )*DimensionNum);
	double *fPopulation = (double *) malloc(sizeof(double)*PopulationNum);
	struct Sort_ref *pArr  = (struct Sort_ref *) malloc(sizeof(struct Sort_ref)*PointNum     );

		time_t start_time = time(NULL);
		double minFval = 99999999.0;
		int noImprove = 0;
		
		for( int itIdx = 0; itIdx < IterationNum; itIdx++ )
		{
			//--- initialize population
			if(type > 0)		// initial OAD using SRS of full design
				InitializePopulation_OAD(PointNum, FactorNum, PosNum, DimensionNum, PopulationNum, Population, full, fullsize);
			else  // LHD
				InitializeSwarm(PointNum, FactorNum, PosNum, DimensionNum, PopulationNum, Population, pArr, type);

			//--- compute function values
			ComputeFVal_type(PointNum, FactorNum, PosNum, DimensionNum, PopulationNum, PValue, Population, fPopulation, type, full, fullsize);

			double fBest;
			int minIdx = minValIdx(PopulationNum, fPopulation, &fBest);
//			printf("iter=%d: fBest=%lf minFval=%lf\n", itIdx, fBest, minFval);

			// save the best design
			if( itIdx == 0  || fBest < minFval)
			{
				minFval = fBest; // fBest = fPopulation[minIdx]
				memcpy(Best, Population + minIdx*DimensionNum, sizeof(int)*DimensionNum);
				noImprove = 0;
			}
			else ++noImprove;
			if(noImprove >= maxNoImprove)	break;	// end search after maxNoImprove

			time_t cur_time = time(NULL);
			if(difftime(cur_time, start_time) >= maxTime ) 	break; // end search after 30 minutes

		} //--- end of loop

	//------------------------//
	// convert Best to xBest:
	//  LHD: PointNum * FactorNum; OAD: FactorNum *PosNum
	//------------------------//
	NumericMatrix xBest = Convert2Design(PointNum, FactorNum, PosNum, PValue, Best, type);

	// memory deallocation ---
	free(Population);
	free(Best);
	free(fPopulation);
	free(pArr);
	if(full != NULL)	free(full);
	// -----------------------

	return xBest;
}


void InitializePopulation_OAD(int PointNum, int FactorNum, int PosNum, int DimensionNum, int PopulationNum, int *Population, int *full, int fullsize)
{	// initial OAD with SRS of the full design

	for( int plIdx = 0; plIdx < PopulationNum; plIdx++ ){
		// take an SRS
		IntegerVector srs = sample(fullsize, FactorNum, false) - 1;	// R style

		for( int ftIdx = 0; ftIdx < FactorNum; ftIdx++ )
		{
			// only first PosNum positions are copied and used if PosNum < PointNum
			memcpy(Population + plIdx*DimensionNum + ftIdx*PointNum,  full +  srs[ftIdx] * PosNum, sizeof(int)*PosNum );		// full is a vector of fullsize * PosNum
		}

	}

	return;
}


void SelectSurvivor(int PointNum, int FactorNum, int DimensionNum, int PopulationNum, int *Population, double *fPopulation, int *Survivor, int *Best, double *fBest, struct Sort_ref *sortRef)
{
	//--- initialize the sorting reference
	for( int plIdx = 0; plIdx < PopulationNum; plIdx++ )
	{
		sortRef[plIdx].fval = fPopulation[plIdx];
		sortRef[plIdx].index = plIdx;
	}

	//--- sorting
	qsort(sortRef, PopulationNum, sizeof(struct Sort_ref), compare);

	//--- select survivor
	for( int plIdx = 0; plIdx < (PopulationNum>>1); plIdx++ )
		memcpy(Survivor + plIdx*DimensionNum, Population + (sortRef[plIdx].index)*DimensionNum, sizeof(int)*DimensionNum);

	return;
}


void CrossOver(int PointNum, int FactorNum, int DimensionNum, int PopulationNum, int *Survivor, int *Population, int *randIdx)
{ // as Liefvendahl, M. and Stocki, R. (2006).

	int plIdx;
	int ftIdx;
	int halfPopulationNum = PopulationNum>>1;

	//--- first half children -- copy the best for each child
	for( plIdx = 0; plIdx < halfPopulationNum; plIdx++ )
		memcpy(Population + plIdx*DimensionNum, Survivor, sizeof(int)*DimensionNum);

	//--- second half children -- copy top half population
	memcpy(Population + halfPopulationNum*DimensionNum, Survivor, sizeof(int)*halfPopulationNum*DimensionNum);

	//--- cross over a random factor
	for( plIdx = 0; plIdx < halfPopulationNum; plIdx++ )
	{
		//--- first half -- crossover with each of top half population
		ftIdx = rand() % (FactorNum);	// same random factor for all the half and the best	
		memcpy(Population + plIdx*DimensionNum + ftIdx*PointNum, Survivor + plIdx*DimensionNum + ftIdx*PointNum, sizeof(int)*PointNum);
		#ifdef DEBUG
		randIdx[plIdx] = ftIdx;
		#endif

		//--- second half -- crossover with the best
		ftIdx = rand() % (FactorNum);	// same random factor for all the half and the best
		memcpy(Population + (plIdx + halfPopulationNum)*DimensionNum + ftIdx*PointNum, Survivor + ftIdx*PointNum, sizeof(int)*PointNum);
		#ifdef DEBUG
		randIdx[plIdx + halfPopulationNum] = ftIdx;
		#endif
	}

	return;
}


void Mutate(int PointNum, int FactorNum, int PosNum, int DimensionNum, int PopulationNum, double pMut, int *Population)
{ // mutate all (except for the best one, plIdx=0)
	for( int plIdx = 1; plIdx < PopulationNum; plIdx++ )
		for( int ftIdx = 0; ftIdx < FactorNum; ftIdx++ )
		{
			int offset = plIdx*DimensionNum + ftIdx*PointNum;
			int * Particle = Population + offset;

			double prob = (double)rand()/RAND_MAX;
			if( prob < pMut ) 	Swap2Elements(PosNum, PointNum, Particle);
		}
}

