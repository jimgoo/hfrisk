/*
  Author: Jimmie Goode
  Created: 2012-09-01
*/

#include "myskewt.hpp"

//-------------------------------------------------------------------------------


int fxcalls = 0;


//-------------------------------------------------------------------------------


// string arr2str(const unsigned n, const double* x) {
//   std::ostringstream ss;
//   ss << std::fixed << std::setprecision(5) << "[";  
//   for (int i = 0; i < n; i++) {
// 	ss << x[i];
// 	if (i != n-1) ss << ", ";
//   }
//   ss << "]";
//   return ss.str();
// }

//-------------------------------------------------------------------------------
//
/*
  void myskewt::setTolerances(START_TYPE startType) {
  if (startType == COLD) {

  TOL_DIFF  = 0.000001;
  TOL_GRAD  = 0.0001;
  STEP_SIZE = 0.01;
  TOL_MIN   = 0.0001;
  MAX_ITER  = 200;
	
  } else if (startType == HOT) {

  TOL_DIFF  = 0.000001;
  TOL_GRAD  = 0.0001;
  STEP_SIZE = 0.01;
  TOL_MIN   = 0.001;
  MAX_ITER  = 100;
	
  } else {
  assert(0 && "Unknown start type\n");
  }
  }*/

//-------------------------------------------------------------------------------

double sign(double x) {
  if (x < 0.0)
	return -1.0;
  else if (x > 0.0)
	return 1.0;
  else if (x == 0.0)
	return 0.0;
  else {
	string msg = "Bad value for sign, x = " + boost::lexical_cast<string>(x);
	cout << msg << endl;
	return 0.0;
  }
}


//-------------------------------------------------------------------------------
//
  
void printPars(DistPars p) {
  printf("asymTpars = [%g, %g, %g, %g]\n", p.gamma, p.mu, p.df, p.sigma);
}

//-------------------------------------------------------------------------------
//
  
void myskewt::DistPars_print(DistPars p) {
  printf("asymTpars = [%g, %g, %g, %g], algoPars = [%g, %u, %g, %g]\n",
		 p.gamma, p.mu, p.df, p.sigma,
		 p.iters, p.fxcalls, p.LLF, p.status);
}

//-------------------------------------------------------------------------------
//
  
void printGradient(double dGamma, double dMu, double dDf, double dSigma) {
  printf("Gradient = [%g, %g, %g, %g]\n", dGamma, dMu, dDf, dSigma);
}

//-------------------------------------------------------------------------------
//

DistPars copyDistPars(DistPars dpIn) {
  DistPars dpOut;

  dpOut.gamma  = dpIn.gamma;
  dpOut.mu     = dpIn.mu;
  dpOut.df     = dpIn.df;
  dpOut.sigma  = dpIn.sigma;
  dpOut.iters  = dpIn.iters;
  dpOut.LLF    = dpIn.LLF;
  dpOut.status = dpIn.status;

  return dpOut;
}

//-------------------------------------------------------------------------------
// Pre estimation tranform of paramaeters

DistPars preEstTransform(START_TYPE startType, gsl_vector* vnLogrets) {

  DistPars pout;

  if (startType == COLD) {
	
	Moments moments = Stats::xMomentSummary(vnLogrets);
	double df = 5.0;
	pout.gamma = sign(moments.skewness) * 0.0001;
	pout.mu    = moments.mean - pout.gamma*df/(df - 2.0);
	pout.df    = log(df - 4.0);
	pout.sigma = log(moments.std);

  } else if (startType == HOT) {

	// pout.gamma = initialParsPrev.gamma;
	// pout.mu    = initialParsPrev.mu;
	// pout.df    = log(initialParsPrev.df - 4.0);
	// pout.sigma = log(initialParsPrev.sigma);

	assert(0 && "needs implementation");
	
  } else {
	assert(0 && "Unknown START_TYPE\n");
  }

  return pout;
}

//-------------------------------------------------------------------------------
// Transform of parameters by Likelihood function during estimation

DistPars duringEstTransform(DistPars pin) {
  DistPars pout;

  pout.gamma = pin.gamma;
  pout.mu = pin.mu;
  pout.df = exp(pin.df) + 4.0;
  pout.sigma = exp(pin.sigma) + 1e-10;

  return pout;
}

//-------------------------------------------------------------------------------
//  Post estimation tranform of paramaeters

DistPars postEstTransform(DistPars pin) {
  DistPars pout;
  pout.gamma  = pin.gamma;
  pout.mu     = pin.mu;
  pout.df     = std::max(1.0, boost::math::round(exp(pin.df))) + 4.0;
  pout.sigma  = exp(pin.sigma);
  pout.iters  = pin.iters;
  pout.LLF    = pin.LLF;
  pout.status = pin.status;

  //delete pin;
  return pout;
}

//-------------------------------------------------------------------------------
// Set initial guesses for MLE

DistPars getInitialGuess(START_TYPE startType, gsl_vector *vnLogrets) {
  return preEstTransform(startType, vnLogrets);
}

//******************************************************************************
// BFGS - GSL
//******************************************************************************

//-------------------------------------------------------------------------------
// Called to get negative log likelihood
double myskewt::negLLF(const unsigned n, const double* x, const gsl_vector* y) {

  fxcalls++;
  
  // during estimation transform
  double gamma = x[0];
  double mu = x[1];
  double df = exp(x[2]) + 4.0;
  double sigma = exp(x[3]) + 1e-10;

  // Cap the degrees of freedom
  if (df > 200.0) {
	df = 200.0;
  }

  double LL;
  int nobs = y->size;
  
  gsl_vector* logPDF = gsl_vector_alloc(nobs);
  gsl_vector* negOne = gsl_vector_alloc(nobs);
  gsl_vector_set_all(negOne, -1.0);

  try {
	skeweduvstudenttpdf(logPDF, y, gamma, mu, df, sigma, 1);
	gsl_blas_ddot(negOne, logPDF, &LL);	
  } catch (...) {
	// Catch any and all exceptions by setting LLF to very large value.
	// This includes cases where DOF is very large, causing Boost Bessel functions to throw errors.
    LL = 1e+20;
  }

  if (gsl_finite(LL) != 1) {
	LL = 1e+20;
  }

  gsl_vector_free(negOne);
  gsl_vector_free(logPDF);
  
  return LL;
}

//-------------------------------------------------------------------------------
//
  
DistPars pars2struct(double gamma, double mu, double df, double sigma) {
  DistPars p;
  p.gamma = gamma;
  p.mu = mu;
  p.df = df;
  p.sigma = sigma;
  return p;
}

//-------------------------------------------------------------------------------
// GSL version of "negLLF"

double my_f(const gsl_vector *x, void *p) {
  gsl_vector* y = (gsl_vector *) p;
  return myskewt::negLLF(x->size, x->data, y); 
}

//-------------------------------------------------------------------------------
// 

typedef struct {
  const gsl_vector *x;
  gsl_vector *y;
  int idx;
} diff_struct;

//-------------------------------------------------------------------------------
//

double myskewt::negLLF_par(double x, void *params) {
  diff_struct* s = (diff_struct *) params;

  int n = s->x->size;
  
  // copy existing x[] array into xx[] so that it can be changed.
  double xx[n];
  for (int i = 0; i < n; i++)
	xx[i] = gsl_vector_get(s->x, i);

  // set the idx element to double x.
  xx[s->idx] = x;
  
  return myskewt::negLLF(n, xx, s->y);
}


//-------------------------------------------------------------------------------
// GSL: The gradient of f, df = (df/dx_1, df/dx_2, ...).

void my_df(const gsl_vector *x, void *params, gsl_vector *grad) {

  gsl_vector* y = (gsl_vector *) params;

  diff_struct ds;
  ds.y = y;
  ds.x = x;

  double dx, dxerr;

  gsl_function F;  
  double tol =  0.000001; //1.0e-6;

  for (int i = 0; i < x->size; i++) {
	ds.idx = i;
	F.function = &myskewt::negLLF_par;
	F.params = &ds;
  
    gsl_deriv_central(&F, gsl_vector_get(x,i), tol, &dx, &dxerr);
    gsl_vector_set(grad, i, dx);
  }

  //double llf = negLLF(x->size, x->data, y);
  //cout << " grad = " << arr2str(grad->size, grad->data) << ", -llf = " << llf << endl;
}

//-------------------------------------------------------------------------------
// GSL: Compute both f and df together.

void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *grad) {
  *f = my_f(x, params);
  my_df(x, params, grad);
}


//-------------------------------------------------------------------------------
// BFGS method with finite differencing of gradient

DistPars myskewt::skewedstudenttfit_bfgs(START_TYPE startType, gsl_vector *y) {

  fxcalls = 0;
  
  double TOL_GRAD  = 0.0001;
  double STEP_SIZE = 0.01;
  double TOL_MIN   = 0.0001;
  double MAX_ITER  = 200;
  
  DistPars distPars; // initialize output structure

  // Create function to be minimized (with gradient)
  int dimension = 4;
  gsl_multimin_function_fdf my_func;
  my_func.n = dimension;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = (void *) y;

  // Set initial starting points
  DistPars x0;
  x0 = getInitialGuess(startType, y);
  
  gsl_vector* x = gsl_vector_alloc(dimension); // should be free;
  gsl_vector_set(x, 0, x0.gamma);
  gsl_vector_set(x, 1, x0.mu);
  gsl_vector_set(x, 2, x0.df);
  gsl_vector_set(x, 3, x0.sigma);

  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  
  T = gsl_multimin_fdfminimizer_vector_bfgs2;
  s = gsl_multimin_fdfminimizer_alloc(T, dimension);
  gsl_multimin_fdfminimizer_set(s, &my_func, x, STEP_SIZE, TOL_MIN);

  int status;
  int iter = 0;
  
  do {
	iter++;
	status = gsl_multimin_fdfminimizer_iterate(s);

	if (status) {
	  break;
	}

	status = gsl_multimin_test_gradient(s->gradient, TOL_GRAD);

  } while (status == GSL_CONTINUE && iter < MAX_ITER);
  
  distPars.gamma  = gsl_vector_get(s->x, 0);
  distPars.mu     = gsl_vector_get(s->x, 1);
  distPars.df     = gsl_vector_get(s->x, 2);
  distPars.sigma  = gsl_vector_get(s->x, 3);
  distPars.iters  = iter;
  distPars.fxcalls = fxcalls;
  distPars.LLF    = myskewt::negLLF(x->size, s->x->data, y);
  distPars.status = status;
	
  // final transformed parameters
  distPars = postEstTransform(distPars);

  distPars.df = std::min(distPars.df, 200.0);
  
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);

  return distPars;
}


//******************************************************************************
// Nelder-Mead Simplex (Gradient-free) - GSL
//******************************************************************************

// Nelder-Mead simplex method
DistPars myskewt::skewedstudenttfit_nmsimplex(START_TYPE startType, gsl_vector *y) {

  fxcalls = 0;
  
  double MAX_ITER  = 2000;
  //double STEP_SIZE = 1.0e-6;

  DistPars distPars; // initialize output structure

  // Create function to be minimized (with gradient)
  int dimension = 4; // number of parameters;
  
  gsl_multimin_function functionInfo;
  functionInfo.f = &my_f;
  functionInfo.n = dimension;
  functionInfo.params = (void *) y;
  
  // Set initial starting points
  DistPars x0;
  x0 = getInitialGuess(startType, y);
  
  gsl_vector* x = gsl_vector_alloc(dimension); // should be free;
  gsl_vector_set(x, 0, x0.gamma);
  gsl_vector_set(x, 1, x0.mu);
  gsl_vector_set(x, 2, x0.df);
  gsl_vector_set(x, 3, x0.sigma);

  const gsl_multimin_fminimizer_type *T;
  gsl_multimin_fminimizer *s;
  
  T = gsl_multimin_fminimizer_nmsimplex;
  s = gsl_multimin_fminimizer_alloc(T, dimension);

  // set step sizes for each variable
  gsl_vector* stepSizes = gsl_vector_alloc(dimension);

  gsl_vector_set(stepSizes, 0, 1.0e-6);
  gsl_vector_set(stepSizes, 1, 1.0e-6);
  gsl_vector_set(stepSizes, 2, 1.0e-2);
  gsl_vector_set(stepSizes, 3, 1.0e-6);
  
  //gsl_vector_set_all(stepSizes, STEP_SIZE); // <TODO> what are the best simplex starting steps?
  gsl_multimin_fminimizer_set(s, &functionInfo, x, stepSizes);

  gsl_vector_free(stepSizes);
  
  //bool isSuccess = false;
  int status;
  int iter = 0;
  
  do {
	iter++;
	status = gsl_multimin_fminimizer_iterate(s);

	if (status) {
	  break;
	}
	
  } while (status == GSL_CONTINUE && iter < MAX_ITER);
  
  distPars.gamma   = gsl_vector_get(s->x, 0);
  distPars.mu      = gsl_vector_get(s->x, 1);
  distPars.df      = gsl_vector_get(s->x, 2);
  distPars.sigma   = gsl_vector_get(s->x, 3);
  distPars.iters   = iter;
  distPars.fxcalls = fxcalls;
  distPars.LLF     = myskewt::negLLF(x->size, s->x->data, y);
  distPars.status  = status;

  //printf("beforePostEst: \n");
  //printPars(distPars);
  
  // final transformed parameters
  distPars = postEstTransform(distPars);
  distPars.df = min(distPars.df, 200.0);
  
  gsl_multimin_fminimizer_free(s);
  gsl_vector_free(x);
  //  gsl_vector_free(stepSizes); // These are already freed by _free(s).

  return distPars;
}


//-------------------------------------------------------------------------------
// Called to get negative log likelihood

double myskewt::negLLF_nlopt(unsigned n, const double *x, double *grad, void *params) {

  fxcalls++;
  
  double gamma = x[0];
  double mu    = x[1];
  double df    = x[2];
  double sigma = x[3];

  gsl_vector* y = (gsl_vector *) params;
  int nobs = y->size;
  
  gsl_vector* logPDF = gsl_vector_alloc(nobs);
  gsl_vector* negOne = gsl_vector_alloc(nobs);
  gsl_vector_set_all(negOne, -1.0);

  double negLL;

  // Catch any and all exceptions by setting LLF to very large value.
  try {
	skeweduvstudenttpdf(logPDF, y, gamma, mu, df, sigma, 1);
	gsl_blas_ddot(negOne, logPDF, &negLL);	
  } catch (...) {
    negLL = 1.0e+20;
  }

  if (gsl_finite(negLL) != 1) {
	negLL = 1.0e+20;
  }

  gsl_vector_free(negOne);
  gsl_vector_free(logPDF);

  return negLL;
}


//-------------------------------------------------------------------------------

double myskewt::negLLF_nlopt_std(unsigned n, const double *x, double *grad, void *params) {

  fxcalls++;
  
  double gamma = x[0];
  double df    = x[1];

  double mu, sigma;
  myskewt::std_pars(gamma,df,mu,sigma);

  gsl_vector* y = (gsl_vector *) params;
  int nobs = y->size;
  
  gsl_vector* logPDF = gsl_vector_alloc(nobs);
  gsl_vector* negOne = gsl_vector_alloc(nobs);
  gsl_vector_set_all(negOne, -1.0);

  double negLL;

  // Catch any and all exceptions by setting LLF to very large value.
  try {
	skeweduvstudenttpdf(logPDF, y, gamma, mu, df, sigma, 1);
	gsl_blas_ddot(negOne, logPDF, &negLL);	
  } catch (...) {
    negLL = 1.0e+20;
  }

  if (gsl_finite(negLL) != 1) {
	negLL = 1.0e+20;
  }

  gsl_vector_free(negOne);
  gsl_vector_free(logPDF);

  return negLL;
}


//-------------------------------------------------------------------------------
// Fit with NLOPT

DistPars myskewt::fit_nlopt(gsl_vector *y) {

  int iPars = 4;
  fxcalls = 0;
  
  // <TODO> Free these?
  vector<double>  x(iPars); // parameter vector
  vector<double> lb(iPars); // lower param bounds
  vector<double> ub(iPars); // upper param bounds
  
  Moments mom = Stats::xMomentSummary(y);

  double tol = 1.0e-6;
  double BIG = 1.0e6;
  
  // gamma
  x[0]  = sign(mom.skewness) * 0.0001;
  lb[0] = -BIG;
  ub[0] =  BIG;
  
  // mu
  x[1]  = mom.mean;
  lb[1] = -BIG;
  ub[1] =  BIG;

  // df
  x[2]  = 5.0;
  lb[2] = 4.0 + tol;
  ub[2] = 200.0;

  // sigma
  x[3]  = mom.std;
  lb[3] = tol;
  ub[3] = BIG;

  nlopt::algorithm algo = nlopt::LN_NELDERMEAD;
  nlopt::opt opt(algo, iPars);
    
  opt.set_min_objective(myskewt::negLLF_nlopt, y);
  opt.set_xtol_rel(1e-6);
  opt.set_ftol_abs(1e-6);
  opt.set_maxeval(5000);
  //opt.set_stopval(-9000.0);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

  double minf;
  double time = MPI::Wtime();
  nlopt::result nloptResult;
		
  nloptResult = opt.optimize(x, minf);

  DistPars dp;
  dp.gamma   = x[0];
  dp.mu      = x[1];
  dp.df      = x[2];
  dp.sigma   = x[3];
  dp.iters   = -1;
  dp.fxcalls = fxcalls;
  dp.LLF     = minf;
  dp.status  = (int) nloptResult;
  dp.runTime = MPI::Wtime() - time;
  dp.iters   = 0;

  return dp;
}


//-------------------------------------------------------------------------------

void myskewt::std_pars(double gamma, double df, double &mu, double &sigma) {

  mu = -1.0*(df/(df - 2.0))*gamma;
  sigma = ((df - 2.0)/df)*(1. - 2.*df*df/(pow(df - 2., 2.)*(df - 4.))*gamma*gamma);
  
}


//-------------------------------------------------------------------------------

DistPars myskewt::fit_nlopt_std(gsl_vector *y) {

  int iPars = 2;
  fxcalls = 0;
  
  // <TODO> Free these?
  vector<double>  x(iPars); // parameter vector
  vector<double> lb(iPars); // lower param bounds
  vector<double> ub(iPars); // upper param bounds
  
  Moments mom = Stats::xMomentSummary(y);

  double tol = 1.0e-6;
  double BIG = 1.0e6;
  
  // gamma
  x[0]  = sign(mom.skewness) * 0.0001;
  lb[0] = -BIG;
  ub[0] =  BIG;

  // df
  x[1]  = 5.0;
  lb[1] = 4.0 + tol;
  ub[1] = 200.0;

  nlopt::algorithm algo = nlopt::LN_NELDERMEAD;
  nlopt::opt opt(algo, iPars);
    
  opt.set_min_objective(myskewt::negLLF_nlopt_std, y);
  opt.set_xtol_rel(1e-6);
  opt.set_ftol_abs(1e-6);
  opt.set_maxeval(5000);
  //opt.set_stopval(-9000.0);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

  double minf;
  double time = MPI::Wtime();
  nlopt::result nloptResult;
		
  nloptResult = opt.optimize(x, minf);

  DistPars dp;
  dp.gamma   = x[0];
  dp.df      = x[1];
  dp.iters   = -1;
  dp.fxcalls = fxcalls;
  dp.LLF     = minf;
  dp.status  = (int) nloptResult;
  dp.runTime = MPI::Wtime() - time;
  dp.iters   = 0;

  return dp;
}


//-------------------------------------------------------------------------------
// Mean(X) when X ~ t_d(gamma, mu, df, C)

vec myskewt::mean(vec gamma, vec mu, double df,  mat C) {

  return mu + (df / (df - 2.0)) * gamma;
}

//-------------------------------------------------------------------------------
// Covariance(X) when X ~ t_d(gamma, mu, df, C)

mat myskewt::cov(vec gamma, vec mu, double df,  mat C) {

  double df2 = df - 2.0;
  return (df/df2)*C + (2*df*df/(df2*df2*(df-4.0)))*gamma*gamma.t();
}
	

//-------------------------------------------------------------------------------
// Generate samples from X when X ~ t_d(gamma, mu, df, C) for d > 1

mat myskewt::mvskewtrnd_1(vec gamma, vec mu, double df,  mat C, int nSim) {

  int nVars = gamma.n_rows;

  double t_mvn = MPI::Wtime();
  mat Z = Stats::mvnrnd(arma::zeros(nVars), C, nSim);
  t_mvn = MPI::Wtime() - t_mvn;

  double t_gam = MPI::Wtime();
  mat W = df / Stats::gamrnd(df/2.0, 2, nSim, nVars);
  t_gam = MPI::Wtime() - t_gam;
  
  double t_mult = MPI::Wtime();
  mat sims = repmat(mu.t(),nSim,1) + repmat(gamma.t(),nSim,1) % W + arma::sqrt(W) % Z;
  t_mult = MPI::Wtime() - t_mult;
  
  //cout << "---> Runtimes for mvskewtrnd_1: t_mvn = " << t_mvn
  // << ", t_gam = " << t_gam << ", t_mult = " << t_mult << endl;
  
  return sims;
}

//-------------------------------------------------------------------------------
// Generate samples from X when X ~ t_d(gamma, mu, df, C) for d > 1
// Uses Cholesky factor A.

/*
mat myskewt::mvskewtrnd_2(vec gamma, vec mu, double df, int nSim, mat A) {

  int nVars = gamma.n_rows;

  double t_mvn = MPI::Wtime();
  mat Z = Stats::mvnrnd_2(arma::zeros(nVars), nSim, A);
  t_mvn = MPI::Wtime() - t_mvn;

  cout << "t_mvn = " << t_mvn << endl;
  
  double t_gam = MPI::Wtime();
  mat W = df / Stats::gamrnd(df/2.0, 2, nSim, nVars);
  t_gam = MPI::Wtime() - t_gam;

  cout << "t_gam = " << t_gam << endl;
  
  double t_mult = MPI::Wtime();
  mat sims = repmat(mu.t(),nSim,1) + repmat(gamma.t(),nSim,1) % W + arma::sqrt(W) % Z;
  t_mult = MPI::Wtime() - t_mult;

  cout << "t_mult = " << t_mult << endl;
  
  //cout << "---> Runtimes for mvskewtrnd_2: t_mvn = " << t_mvn
  //   << ", t_gam = " << t_gam << ", t_mult = " << t_mult << endl;
  
  return sims;
}
*/


//-------------------------------------------------------------------------------

mat myskewt::mvskewtrnd_2(vec gamma, vec mu, double df, int nSim, mat A) {

  int nVars = gamma.n_rows;

  mat W = df / Stats::gamrnd(df/2.0, 2, nSim, nVars);
  
  mat sims = repmat(mu.t(),nSim,1)
	+ (repmat(gamma.t(),nSim,1) % W)
	+ (arma::sqrt(W) % Stats::mvnrnd_2(arma::zeros(nVars), nSim, A));

  return sims;
}

//-------------------------------------------------------------------------------
// Generate samples from X when X ~ t_1(gamma, mu, df, C)

mat myskewt::skewtrnd(int nRows, int nCols, double gamma, double mu, double df, double sigma) {

  mat Z = Stats::randn(0, 1, nRows, nCols);
  mat W = df / Stats::gamrnd(df/2.0, 2, nRows, nCols);
  return mu + gamma*W + sigma*(arma::sqrt(W) % Z);
}

//-------------------------------------------------------------------------------
//

void myskewt::mvskewtrnd_test() {

  vec gamma(2);
  vec mu(2);
  double df = 5.0;
  mat C;
	
  gamma(0) = 0.1;
  gamma(1) = -0.1;

  mu(0) = 0.01;
  mu(1) = 0.02;

  C << 1.5 << 0.1 << endr
	<< 0.1 << 2.0 << endr;
  
  mat samp = myskewt::mvskewtrnd_1(gamma, mu, df, C, 10000);
  //samp.print("samp = ");

  // analytic moments
  vec mu_ana = myskewt::mean(gamma, mu, df, C);
  mat C_ana  = myskewt::cov(gamma, mu, df, C);
  mu_ana.print("mu_ana = ");
  C_ana.print("C_ana = ");
 
  // empirical moments
  mat mu_hat = arma::mean(samp, 0).t();
  mat C_hat = arma::cov(samp);
  mu_hat.print("mu_hat = ");
  C_hat.print("C_hat = ");

}


//-------------------------------------------------------------------------------
// CDF for skewed t

vec myskewt::cdf(vec x, int iSamp, double gamma, double mu, double df, double sigma) {
  
  mat sample = myskewt::skewtrnd(iSamp, 1, gamma, mu, df, sigma);
  
  return Stats::empCDF(sample.col(0), x);
  
}

//-------------------------------------------------------------------------------
// Inverse CDF for skewed t

// THIS ISN'T WORKING - I think emp_CDF call needs to sort cdfi.

vec myskewt::inv_cdf(vec pvalues, int iN, double gamma, double mu, double df, double sigma) {

  mat sample = myskewt::skewtrnd(iN, 1, gamma, mu, df, sigma);  
  vec q(pvalues.n_rows);
  vec x(iN);
  vec cdfi(iN);
  Stats::empCDF_fast(sample.col(0), x, cdfi); // problem

  // tail cropping tolerance
  double tol = 1.0e-6;

  cdfi.save("/Users/jimmiegoode/Documents/Glimm/Toolbox/Copulas/e_cdf_pre.csv", csv_ascii);
  x.save("/Users/jimmiegoode/Documents/Glimm/Toolbox/Copulas/e_x_pre.csv", csv_ascii);
  
   // force CDF to be monotonic
  Stats::cleanupCDF(x, cdfi, tol);
  
  double xMax = max(x);
  double xMin = min(x);
  double cdfMin = min(cdfi);
  double cdfMax = max(cdfi);
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, iN);
  gsl_interp_init(interp, cdfi.memptr(), x.memptr(), iN);

  /*
  for (int i = 0; i < pvalues.n_rows; i++) {
	if (pvalues(i) >= cdfMax)
	  q(i) = xMax;
	else if (pvalues(i) <= cdfMin)
	  q(i) = xMin;
	else
	  q(i) = gsl_interp_eval(interp, cdfi.memptr(), x.memptr(), pvalues(i), acc);
  }
  */
  gsl_interp_accel_free(acc);
  gsl_interp_free(interp);

  return q;
}

//-------------------------------------------------------------------------------
// Estimates the matrix Sigma in the stochastic form for t_d.

void myskewt::estimateSigma(const mat sample,
							const vec gamma,
							const vec mu,
							const double df,
							const bool doShrink,
							double &shrinkage,
							mat &C) {

  // Shrink the covariance
  shrinkage = -1;

  if (doShrink) {
	C = Stats::cov2para(sample, shrinkage);
  }

  C = C*((df - 2.0)/df) - gamma*gamma.t() * (2.0*df)/((df - 2.0)*(df - 4.0));
  
}

//-------------------------------------------------------------------------------
// NOTE: This updates resid in place with GARCH filtered returns.

void myskewt::portSample(const int nSim,
							 const vec gamma,
							 const vec mu,
							 const double df,
							 const mat mnGarchPars,
							 const vec wts,
							 const mat A,
							 mat &resid,
							 vec &nextMeans,
							 vec &nextSigmas,
							 vec &portRets) {

  // number of stocks
  int iS = gamma.n_rows;
  
  resid = myskewt::mvskewtrnd_2(gamma, mu, df, nSim, A);

  for (int i = 0; i < iS; i++) {

	vec nxtRet(nSim);

	mygarch::forecast_fromVec(mnGarchPars.col(i), resid.col(i),
							  nxtRet, nextMeans(i), nextSigmas(i));

	// Update resid in place with returns.
	resid.col(i) = nxtRet;
  }

  // portfolio returns for each path
  portRets = resid*wts;

}
  

//-------------------------------------------------------------------------------

void myskewt::tester(int c) {

  cout << "\n---> myskewt tester: case " << c << " <---\n\n";

  switch (c) {
  case 1:
	{
  
	  //string fname = "/Users/jimmiegoode/Dropbox/Projects/Glimm/data/scalingSeries/2.csv";
	  //string fname = "/Users/jimmiegoode/Dropbox/Projects/Glimm/data/NASDAQ_Returns_z.csv";
	  string fname = "/Users/jimmiegoode/Documents/Glimm/svn/GLIMM/trunk/data/exports/u_-1.csv";
  
	  gsl_matrix* mnRet = IO::importCSVmatrix(fname);
	  gsl_vector* vnRet = gsl_vector_alloc(mnRet->size1);
	  gsl_matrix_get_col(vnRet, mnRet, 0);

	  DistPars pars;
	  double MPI_Wtime(); 
	  double t;

	  // t = MPI::Wtime();
	  // pars = myskewt::skewedstudenttfit_bfgs(COLD, vnRet);
	  // myskewt::DistPars_print(pars);
	  // cout << "DT_test_1 = " << MPI::Wtime() - t << endl;

	  t = MPI::Wtime();

	  //pars = myskewt::skewedstudenttfit_nmsimplex(COLD, vnRet);
	  pars = myskewt::fit_nlopt(vnRet);
  
	  myskewt::DistPars_print(pars);

	  cout << "DT_test_1 = " << MPI::Wtime() - t << endl;
	}
	break;

  case 2:
	{
	  cout << "\n CDF \n";
	  
	  vec x = linspace(-20., 20., 1000);
	  vec cdf = myskewt::cdf(x, 1.e6, -1., 0., 5.0, 1.);
	  cdf.save("/Users/jimmiegoode/Documents/Glimm/Toolbox/Copulas/e_cdf.csv", csv_ascii);

	  //vec myskewt::inv_cdf(vec pvalues, int iN, double gamma, double mu, double df, double sigma) {

	  cout << "\n INV \n";
	  
	  vec pvalues = linspace(0., 1., 1000);
	  vec inv = myskewt::inv_cdf(pvalues, 1.e6, -1., 0., 5.0, 1.);
	  inv.save("/Users/jimmiegoode/Documents/Glimm/Toolbox/Copulas/e_inv.csv", csv_ascii);
	}
	break;
	
  } //switch
}


//-------------------------------------------------------------------------------
// <main>

// int main(int argc, char** argv) {
  
//   if (argc < 2)
// 	assert(0 && "Invalid Arguments");
  
//   int testCase = atoi(argv[1]);

//   myskewt::tester(testCase);
  
//   return 0;
// }

