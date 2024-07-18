/* GAxy.cpp: genetic-algorithm.cpp, Genetic Algorithm from Yuhao Yin 
		for constructing BID LHD designs 
	called by myGAxy() in MetaX.R
	9/19/23: removed redunant functions
	9/24/24: add type=-1 (BID) or 0 (maximin)
	
References:
	Liefvendahl, M. and Stocki, R. (2006). A study on algorithms for optimization of Latin hypercubes. Journal of Statistical Planning and Inference, 136, 3231-3247.
	(which uses a common column (k2=k1) and is the same as GA in GAx.cpp.)
 
*/

#include <Rcpp.h>
#include <time.h>
using namespace Rcpp;

#define MAX(a, b) (a)>(b) ? (a) : (b)
#define MIN(a, b) (a)<(b) ? (a) : (b)

// common functions used by SA and GA, defined in mySAxy.cpp
double crt_log(NumericMatrix D, double beta, int type);
NumericMatrix generate_rbdesign(int n, int m, int s);

class Individual
{
public:
  NumericMatrix design;
  double measure;
  double beta;
  int design_type;		// design type
  Individual(NumericMatrix design, double gmp_beta, int type);
  Individual mate(Individual par2, double pMut);
  Individual mate_left(Individual par2, double pMut);
  Individual mate_right(Individual par2, double pMut);
  double compute_criterion();
};

Individual::Individual(NumericMatrix design, double gmp_beta, int type)
{
  this->design = design;
  beta = gmp_beta;
  design_type = type;
  measure = compute_criterion();
};

// Perform mating and produce new offspring
Individual Individual::mate(Individual par2, double pMut)
{
  int n=design.rows(), m=design.cols();
  int s = max( design ) + 1;	// number of levels
  NumericMatrix child(n, m);
  
  for(int j=0; j<m; j++){
    double u = runif(1)[0];
    // with 45% prob inherit from parent 1
    if(u<(1-pMut)/2){
      child(_, j) = design(_, j);
    }
    // with another 45% prob inherit from parent 2
    else if(u<1-pMut){
      child(_, j) = par2.design(_, j);
    }
    // with 10% prob mutation occurs
    else{
      child(_, j) = sample(n, n, false)-1;
      /* change to level s */
  	 if(s < n) for (int i=0; i<n; i++) child(i, j) = (int)child(i, j) % s;

    }
  }
  
  return Individual(child, beta, design_type);
};

double Individual::compute_criterion()
{
  return crt_log(design, beta, design_type);
};


Individual Individual::mate_left(Individual par2, double pMut)
{
  NumericMatrix child = clone(design);
  
  // replace a randomly chosen column, k1, in par1
  // with a randomly chosen column, k2, in par2
  int n=design.rows(), m=design.cols();
  int k1 = sample(m, 1)[0]-1;
  int k2 = k1;  // k2 = sample(m, 1)[0]-1;   // use same k 9/20/23
  
  child(_, k1) = par2.design(_, k2);
  
  for(int k=0; k<m; k++)
  {
    double u = runif(1)[0];
    // mutation: randomly switching two elements in a column
    if(u < pMut)
    {
      IntegerVector v = sample(n, 2, false)-1;
      int i = v[0], j = v[1];
      
      double tmp = child(i, k);
      child(i, k) = child(j, k);
      child(j, k) = tmp;
    }
  }
  
  return Individual(child, beta, design_type);
};


Individual Individual::mate_right(Individual par2, double pMut)
{
  NumericMatrix child = clone(par2.design);
  
  // replace a randomly chosen column, k1, in par2
  // with a randomly chosen column, k2, in par1
  int n=design.rows(), m=design.cols();
  int k1 = sample(m, 1)[0]-1;
  int k2 = k1;  // k2 = sample(m, 1)[0]-1;   // use same k 9/20/23
  
  child(_, k1) = design(_, k2);
  
  for(int k=0; k<m; k++)
  {
    double u = runif(1)[0];
    // mutation: randomly switching two elements in a column
    if(u < pMut)
    {
      IntegerVector v = sample(n, 2, false)-1;
      int i = v[0], j = v[1];
      
      double tmp = child(i, k);
      child(i, k) = child(j, k);
      child(j, k) = tmp;
    }
  }
  
  return Individual(child, beta, design_type);
}

// Overloading < operator
bool operator<(const Individual &ind1, const Individual &ind2)
{
  return ind1.measure < ind2.measure;
};
  


// [[Rcpp::export()]]
NumericMatrix GeneticAlg_YY(int n, int m, int s, double gmp_beta, int popSize, double elitism_ratio, 
	double crossover_ratio, double mutation_prob, int maxiter, int type, bool verbose=false)
{ // GA by Yuhao Yin
	
  std::vector<Individual> curr_gen;
   
  // create initial population, i.e. sets of random LHD
  for(int i=0; i<popSize; i++)
  {
    NumericMatrix design = generate_rbdesign(n, m, s);
    curr_gen.push_back(Individual(design, gmp_beta, type));
  }
  
  for(int g=0; g<maxiter; g++)
  {
    sort(curr_gen.begin(), curr_gen.end());
    
    std::vector<Individual> next_gen;
    
    // keep the top elitism_ratio% individuals (design)
    int elitist_size = MAX(popSize * elitism_ratio, 1);		// at least 1
    for(int i=0; i<elitist_size; i++)
      next_gen.push_back(curr_gen[i]);
    
    // From 50% of fittest population, Individuals
    // will mate to produce offspring
    int mate_size = popSize - elitist_size;
    for(int i=0; i<mate_size; i++)
    {
    	int crossover_size = MAX(int(popSize*crossover_ratio), 2);	// at least 2
      IntegerVector v = sample(crossover_size, 2, false)-1;
      Individual par1 = curr_gen[v[0]];
      Individual par2 = curr_gen[v[1]];
      
      Individual offspring = par1.mate(par2, mutation_prob);
      next_gen.push_back(offspring);
    }
    
    curr_gen = next_gen;
    
    if(verbose && (g % 1000 == 0))
    {
      Rcout<< "Generation=" << g << "; Measure: "<< curr_gen[0].measure << "\n";
    }
  }
  
  return curr_gen[0].design;
};


// [[Rcpp::export()]]
NumericMatrix GeneticAlg_base(int n, int m, int s, double beta, int popSize, double pMut, int maxiter, int type, bool verbose)
{ // GAx1: modified according to Liefvendahl, M. and Stocki, R. (2006).
// same as GA in GAx.cpp 
	
  std::vector<Individual> curr_gen;
  
  // create initial population, i.e. sets of random LHD
  for(int i=0; i<popSize; i++)
  {
    NumericMatrix design = generate_rbdesign(n, m, s);
    curr_gen.push_back(Individual(design, beta, type));
  }
  
//  time_t end = time(NULL) + secs;
  
//  while(time(NULL) < end)
  for(int g=0; g<maxiter; g++)
  {
    sort(curr_gen.begin(), curr_gen.end());
    
    if(verbose){
      Rcout << "crt_log=" << crt_log(curr_gen[0].design, beta, type) << ";\n";
    }
    
    std::vector<Individual> next_gen;
    
    // From 50% of fittest population with the best, Individuals
    // will mate to produce offspring
    Individual par1 = curr_gen[0];
    
    next_gen.push_back(curr_gen[0]);	// copy the best, 9/20/23
    next_gen.push_back(par1.mate_left(curr_gen[0], pMut));	// mutate the best, 9/20/23
    
    int mate_size = popSize / 2;
    for(int i=1; i<mate_size; i++)	// one pair less 9/20/23
    {
      Individual par2 = curr_gen[i];	
      next_gen.push_back(par1.mate_left(par2, pMut));
      next_gen.push_back(par1.mate_right(par2, pMut));
    }
    
    curr_gen = next_gen;
  }
  
  return curr_gen[0].design;
};



// [[Rcpp::export()]]
double test(int a)
{
  return 1.0/a;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
test(2)
*/
