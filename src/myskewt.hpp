/*
  Author: Jimmie Goode
  Created: 2012-09-01
*/

#ifndef __MYSKEWT_HPP__
#define __MYSKEWT_HPP__

// STD
#include <iostream>
#include <math.h>
#include <time.h>
#include <sstream>
#include <iomanip>

// GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_deriv.h> // finite differencing
#include <gsl/gsl_multimin.h> // optimization
#include <gsl/gsl_spline.h> // interp for pdfs
#include <gsl/gsl_interp.h> // interp for pdfs

//BOOST
#include <boost/math/special_functions/round.hpp>
#include <boost/math/special_functions/bessel.hpp>

// OTHER EXT
#include <armadillo>
#include "mpi.h"
#include <nlopt.hpp>

// LOCAL
#include "Constants.h"
#include "gamma.h"
#include "skeweduvstudenttpdf.h"
#include "IO.hpp"
#include "Stats.hpp"

using namespace std;

// structure to hold univariate skewed t parameters
struct DistPars {
  // distribution parameters
  double gamma;
  double mu;
  double df;
  double sigma;
  // estimation parameters
  double iters;
  int fxcalls;
  double LLF;
  double status;
  double runTime;
  //gsl_vector *y;
};


class myskewt {

public:

  static DistPars skewedstudenttfit_bfgs(START_TYPE startType, gsl_vector *y);
  static DistPars skewedstudenttfit_nmsimplex(START_TYPE startType, gsl_vector *y);

  static double negLLF(const unsigned n, const double* x, const gsl_vector* y);
  static double negLLF_par(double x, void *params);
  
  static gsl_matrix* computeSkewtCopulaCovariance(gsl_vector* gamma,
														   gsl_vector* mu,
														   gsl_matrix* mnReturns);

  static mat mvskewtrnd_1(vec gamma, vec mu, double df, mat C, int nSim);
  static mat mvskewtrnd_2(vec gamma, vec mu, double df, int nSim, mat A);

  static void mvskewtrnd_test();

  static vec mean(vec gamma, vec mu, double df, mat C);
  static mat cov(vec gamma, vec mu, double df, mat C);

  static void DistPars_print(DistPars p);

  static mat skewtrnd(int nRows, int nCols, double gamma, double mu, double df, double sigma);

  static double negLLF_nlopt(unsigned n, const double *x, double *grad, void *params);
  static double negLLF_nlopt_std(unsigned n, const double *x, double *grad, void *params);
  
  static DistPars fit_nlopt(gsl_vector *y);
  static DistPars fit_nlopt_std(gsl_vector *y);

  static void std_pars(double gamma, double df, double &mu, double &sigma);

  static vec cdf(vec X, double gamma, double mu, double df, double sigma);

  static void tester(int c);

  static vec inv_cdf(vec u, double gamma, double mu, double df, double sigma);
};


#endif
