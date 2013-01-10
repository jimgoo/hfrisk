#ifndef __MYSTABLE_HPP__
#define __MYSTABLE_HPP__

// STD
#include <iostream>
#include <math.h>
#include <sstream>
#include <iomanip>

// GSL
#include <gsl/gsl_sf.h> // gamma fx
#include <gsl/gsl_math.h> // pi constant
#include <gsl/gsl_spline.h> // interp for pdfs
#include <gsl/gsl_interp.h> // interp for pdfs

// EXT
#include "mpi.h" // timing
#include <nlopt.hpp> // MLE optimization
#include <armadillo>
#include <complex> 
#include <fftw3.h> // FFT for pdfs
#include "Stats.hpp" // randn, rand


using namespace std;
using namespace arma;


class mystable {
public:
  
  enum dist {stdCTS, stdNTS, stdAS, symAS, AS};

  const static double gridSize = 8192; //2^13 = 8192, 2^15 = 32768
  const static double gridStep = 0.01;

  typedef struct {
	double* pars;
	double iters;
	double fxcalls;
	double LLF;
	double status;
	double runTime;
	mystable::dist dist;
  
  } ts_struct;
  
  static string arr2str(const unsigned n, const double* x);
  static string vec2str(const unsigned n, const vector<double> x);
  
  static ts_struct mle_nlopt(vec y, mystable::dist dist);
  static double assg_alphaEst(vec alpha);
  static mat assg_dispersionEst(mat X, double alpha, vec sigma, vec mu);
  
  static vec cdf_FFT(vec arg, const int nparam, const double param[], mystable::dist dist);
  static vec inv_FFT(vec pvalues, const int nparam, const double param[], mystable::dist dist);

  static int getParCount(mystable::dist dist);

  static string dist2str(mystable::dist dist);

  static double assg_scaleEst(mat X, double alpha, vec mu, vec a);

  static void ts_struct_print(mystable::ts_struct s);
};


#endif

