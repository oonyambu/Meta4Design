/* SAxy.cpp: modified from simulated-annealing.cpp, Simulated annealing from Yuhao Yin 
		for constructing BID LHD designs 
	called by mySAxy() in MetaX.R
	9/19/23: removed redunant functions
	removed SimulatedAnnealing, a special case of MultiSimulatedAnnealing with s=n
	renamed as SA_XY() and SA_MM(), 9/21/23

Note: normalize design D01 = D/q; 10/24/23
 see also bid_log and bid in core.cpp
*/

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export()]]
double gmp_log(NumericMatrix D, double beta)
{ // same as bid_log in core.cpp
  int n=D.rows(), m=D.cols();
  int q=max(D)-min(D)+1;		// number of levels
  NumericMatrix D01 = D/q;			// normalized to [0,1]
  
  double obj = 0;
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      double prod = 1.0;
      for(int k=0; k<m; k++){
        prod *= 1/(beta + pow(D01(i, k) - D01(j, k), 2.0));
      }
      obj += prod;
    }
  }
  
  if(beta != 0)	obj *= pow(beta, m);		// BID
 
  return (1.0/m)*log(2.0/(n*(n-1))*obj);
}

// [[Rcpp::export()]]
double gmp(NumericMatrix D, double beta)
{ // if beta=0, bid is the maxpro criterion

  return exp( gmp_log(D, beta) );
}

// [[Rcpp::export()]]
double gmp_diff(NumericMatrix D, int i0, int j0, int k0, double beta=1.0)
{ 
// compute the difference in gmp_normalized values after switching elements of row i0 and j0

  int n = D.rows(), m = D.cols();
  int q=max(D)-min(D)+1;		// number of levels
  NumericMatrix D01 = D/q;			// normalized to [0,1]
  NumericMatrix D01_prop = clone(D01);
  D01_prop(i0, k0) = D01(j0, k0);
  D01_prop(j0, k0) = D01(i0, k0);
  
  double diff = 0;
  for(int i=0; i<n; i++){
    double prod_i0_new = 1.0, prod_i0_old = 1.0, prod_j0_new = 1.0, prod_j0_old = 1.0;
    for(int k=0; k<m; k++){
      if(i != i0)
      {
        prod_i0_new *= beta/(beta + pow(D01_prop(i, k) - D01_prop(i0, k), 2.0));
        prod_i0_old *= beta/(beta + pow(D01(i, k) - D01(i0, k), 2.0));
      }
      if(i != j0)
      {
        prod_j0_new *= beta/(beta + pow(D01_prop(i, k) - D01_prop(j0, k), 2.0));
        prod_j0_old *= beta/(beta + pow(D01(i, k) - D01(j0, k), 2.0));
      }
    }
    diff += (prod_i0_new - prod_i0_old + prod_j0_new - prod_j0_old);
  }

  return diff;
}

//// phip measure for the maximin criterion, 9/23/23

// [[Rcpp::export()]]
double phip_log(NumericMatrix D, double p=15){
  int n=D.rows(), m=D.cols();
  int q=max(D)-min(D)+1;		// number of levels
  NumericMatrix D01 = D/q;			// normalized to [0,1]
  
  double sum_phip = 0;
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      double dist = 0.0;
      for(int k=0; k<m; k++){
        dist += pow(D01(i, k) - D01(j, k), 2.0);
      }
      sum_phip += pow( sqrt(dist), -p);
    }
  }
  
  return (1.0/p)*log(sum_phip);
}

// [[Rcpp::export()]]
double phip(NumericMatrix D, double p=15){
   
  return exp( phip_log(D, p) );
}


// [[Rcpp::export()]]
double phip_diff(NumericMatrix D, int i0, int j0, int k0, double p=15)
{ 
// compute the difference in phip values after switching elements of row i0 and j0
	if( D(i0, k0) == D(j0, k0) )		return 0;	// the same element, which happens for multi-level designs
	
  int n = D.rows(), m = D.cols();
  int q=max(D)-min(D)+1;		// number of levels
  NumericMatrix D01 = D/q;			// normalized to [0,1]
  NumericMatrix D01_prop = clone(D01);
  D01_prop(i0, k0) = D01(j0, k0);
  D01_prop(j0, k0) = D01(i0, k0);
  
  double diff = 0;
  for(int i=0; i<n; i++){
    double sum_i0_new = 0.0, sum_i0_old = 0.0, sum_j0_new = 0.0, sum_j0_old = 0.0;
    if(i != i0){
        for(int k=0; k<m; k++){
	        sum_i0_new += pow(D01_prop(i, k) - D01_prop(i0, k), 2.0);
	        sum_i0_old += pow(D01(i, k) - D01(i0, k), 2.0);
        } // for k
     diff += pow(sum_i0_new, -p/2.0) - pow(sum_i0_old, -p/2.0);
     }
    if(i != j0){
   	  for(int k=0; k<m; k++){
	        sum_j0_new += pow(D01_prop(i, k) - D01_prop(j0, k), 2.0);
	        sum_j0_old += pow(D01(i, k) - D01(j0, k), 2.0);
	      }  // for k
	 diff += pow(sum_j0_new, -p/2.0) - pow(sum_j0_old, -p/2.0);
    }  
 
 }

  return diff;
}

// [[Rcpp::export()]]
double crt_log(NumericMatrix D, double beta, int type)
{
  	if(type == 0)   return phip_log(D, beta);
  	
  	return gmp_log(D, beta);
}

// [[Rcpp::export()]]
double crt_diff(NumericMatrix D, int i0, int j0, int k0, double beta, int type)
{ 
  	if(type == 0)   return phip_diff(D, i0, j0, k0, beta);
  	
  	return gmp_diff(D, i0, j0, k0, beta);
}

// [[Rcpp::export()]]
NumericMatrix generate_rbdesign(int n, int m, int s)
{ // return a nearly balanced n*m design with s levels
  NumericMatrix random_design(n, m);

  for(int j=0; j<m; j++)
  {
    random_design(_, j) = sample(n, n, false)-1;
    
     /* change to level s */
  	if(s < n) for (int i=0; i<n; i++) random_design(i, j) = (int)random_design(i, j) % s;
  }

  return random_design;
}

// NumericMatrix generate_rdesign(int n, int m){
//	return generate_rbdesign(n, m, n);
// }

void swap_elements(int i, int j, int k, NumericMatrix &D)
{
  double tmp = D(i, k);
  D(i, k) = D(j, k);
  D(j, k) = tmp;
}


// [[Rcpp::export()]]
NumericMatrix SA_XY(int n, int m, int s, double beta, int Ntrials, int itermax, int maxTime, double T0, int Npairs, int nFail, int type, bool verbose)
{ // Modified from MultiSimulatedAnnealing() from Yuhao Yin
// ends at max_total iteration or maxTime or no improvement (accept_flag==false) for 1-3 times
// Npairs is similar to Imax in SA_MM
// type = 0(maximin), -1 (bid)

  	NumericMatrix D_optim = generate_rbdesign(n, m, s);
   	double	min_gmp = crt_log(D_optim, beta, type);
 
	// max total iterations allowed for all trials, consistent with SLHD and DE/PSO/GA
  	int max_total = Ntrials*itermax; 	
 	int itotal = 0;		// total number of iterations
 	time_t end = time(NULL) + maxTime;

	bool b_stop = false;		// whether to stop early
  for(int trial=0; trial<Ntrials && !b_stop; trial++){
    NumericMatrix D = generate_rbdesign(n, m, s);
    double T = T0;
    int accept_count = 0;
    
//	  bool accept_flag = true;	    
	int nFailed = 0;
	while(nFailed <= nFail && !b_stop)	 // run up to itermax or maxTime or failed to improve for 1-3 times.
	{
	bool accept_flag = false;
      
      // randomly choose column k
      int k = sample(m, 1)[0]-1;
      
      double cost_diff_min = R_PosInf;
      int i_min, j_min;
      
      // choose the optimal pair that achieves the lowest cost_diff
      for(int p=0; p<Npairs; p++){
        IntegerVector v = sample(n, 2, false)-1;
        int i = v[0], j = v[1];
        
        double cost_diff_curr = crt_diff(D, i, j, k, beta, type);
        if(cost_diff_curr < cost_diff_min){
          i_min = i;
          j_min = j;
          cost_diff_min = cost_diff_curr;
        }
      }
      
      double u = runif(1)[0];
      if(exp(-(cost_diff_min)/T) >= u){
        swap_elements(i_min, j_min, k, D);
        accept_flag = true;
      }
  
		double gmp_value = crt_log(D, beta, type);
	    if( gmp_value < min_gmp ){	
	      D_optim = clone(D);
	      min_gmp = gmp_value;
	      } 
          
      itotal += Npairs; 	// each step has Npairs swaps
      if(itotal >= max_total || time(NULL) >= end) b_stop=true; // reached itermax or maxTime
 
      // decrease T if a new design is accepted. 
      // This differs from MM, as there is no while(ipert < Imax)
      if(accept_flag){
        accept_count += 1;
        if(T > pow(10, -15))       T *= 0.95;
      }
      	else ++nFailed;    
      
    }  // while (accept_flag) -> while (nFailed <= 3)

    if(verbose){
      	Rcout << "Trial=" << trial << " | " 
      		<< "crt_log=" << min_gmp 
	        << "; accept_count=" << accept_count
	        << "; itotal=" << itotal     
	        << "; T_final: " << T << ";\n";
    }
    
  }
 
  return(D_optim);
}


// [[Rcpp::export()]]
NumericMatrix SA_MM(int n, int m, int s, double beta, int Ntrials, int itermax, int maxTime, double T0, int Imax, int nFail, int type, bool verbose)
{ // ends at itermax or maxTime or no improvement (accept_flag==false -> nFailed <=3)
// This is adopted from Morris and Mitchell (1995) with Imax to control the temperature and stop the search. 9/20/23.

  	NumericMatrix D_optim = generate_rbdesign(n, m, s);
  	double	min_gmp = crt_log(D_optim, beta, type);

	// max total iterations allowed for all trials, consistent with SLHD and DE/PSO/GA
  	int max_total = Ntrials*itermax; 	
 	int itotal = 0;		// total number of iterations
 	time_t end = time(NULL) + maxTime;
 	
 	bool	 b_stop = false; 	// whether to stop early
  	for(int trial=0; trial<Ntrials && !b_stop; trial++){
 
	  NumericMatrix D = generate_rbdesign(n, m, s);
	  double T = T0;
	  int accept_count = 0;
	  	  
//	  bool accept_flag = true;
	int nFailed = 0;
	  while(nFailed <= nFail && !b_stop)	// run up to itermax or maxTime or failed to improve for 1-3 times.
	  {	    
	    bool accept_flag = false;
	    
	    //// constant temperature loop ////
		int ipert=1;
		while (ipert < Imax && !b_stop)
		{						 
		    // randomly choose column k
		    int k = sample(m, 1)[0]-1;
		    
		    // randomly choose row pair (i,j)
		    IntegerVector v = sample(n, 2, false)-1;
		    int i = v[0], j = v[1];
		    
		    double cost_diff = crt_diff(D, i, j, k, beta, type);
		    
		    if(cost_diff<0){
			      swap_elements(i, j, k, D);
			      accept_flag = true;
		    } else{
			      double u = runif(1)[0];
			      if(exp(-(cost_diff)/T) >= u){
			        swap_elements(i, j, k, D);
			        accept_flag = true;
			      }
		    }
		    		    
			double gmp_value = crt_log(D, beta, type);
		    if( gmp_value < min_gmp ){
		      D_optim = clone(D);
		      min_gmp = gmp_value; 
		      ipert = 1; 		// reset ipert=1 
		    }
		    else ++ipert;	// 

			itotal += 1;
	   		if(itotal >= max_total || time(NULL) >= end) b_stop=true; // reached itermax or maxTime

	    }	// while(ipert < Imax), end of constant temperature loop, 

		// decrease T if a new design is accepted
		if(accept_flag){	
	      accept_count += 1;
	      if(T > pow(10, -15))	     T *= 0.95;
	    }
	    else ++nFailed;		// try
	
    	}  // while (accept_flag) -> while (nFailed <= 3)
	   
	  if(verbose){
      	Rcout << "Trial=" << trial << " | " 
      		<< "crt_log=" << min_gmp 
	        << "; accept_count=" << accept_count
	        << "; itotal=" << itotal   << "; T_final: " << T << ";\n";
	    }
	   
	} // for(trial)

  return(D_optim);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R

*/
