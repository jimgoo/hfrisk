#ifndef __MYGARCH_HPP__
#define __MYGARCH_HPP__

// STD
#include <iostream>
#include <math.h>
#include <time.h>
#include <sstream>
#include <iomanip>

// GSL
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_deriv.h> // finite differencing for gradient
#include <gsl/gsl_sf.h> // special function gamma for t-pdf
#include <gsl/gsl_poly.h> // root finding for arma0
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sort.h>

// EXT
#include "mpi.h"
#include <nlopt.hpp>
#include <armadillo>

// LOCAL
#include "IO.hpp"
#include "Stats.hpp"
#include "myskewt.hpp"


typedef struct {
  unsigned n; //
  int R; //
  int M; //
  int P; //
  int Q; //
  int maxRMPQ;//
  string algo;
  double negLL;   //
  int exitStatus; //
  int fxCalls;    //
  int iters;      //
  double runTime;
  
  gsl_vector *yy; // returns
  gsl_vector *ee; // innovations
  gsl_vector *hh; // variances
  gsl_vector *u;  // residuals
  
  double *x; // only used in gradient
  int idx;   // only used in gradient

  gsl_vector *x0; // initial parameters
  gsl_vector *x1; // final parameters

} garch_struct;


using namespace std;
using namespace arma;

class mygarch {

  
  
public:

  static const int iNonCoeffs = 11;
  
  static void garch_struct_print(garch_struct s);
  
  static void getInitialPars(const int n, const gsl_vector* y,
							 const int R, const int M, const int P, const int Q,
							 double x[], double lb[], double ub[]);

  static gsl_vector* negLLF_grad(const unsigned n, const double* x,  garch_struct *s, const double *steps);

  static double negLLF_struct(const double* x, garch_struct *s);

  static garch_struct fit_nlopt(const gsl_vector* y, const int R, const int M, const int P, const int Q);

  static vec forecast(garch_struct s, vec resid);

  //static vec forecast_fromVec(vec x, vec resid);

  static void forecast_fromVec(vec x, vec resid,
								 vec &nextRet, double &nextMean, double &nextSigma);

  static int getMessageSize(garch_struct s);

  static vec garch_struct_vec(garch_struct s);

  static void garch_struct_free(garch_struct s);
};


#endif
