#ifndef __STATS_H__
#define __STATS_H__

// STD
#include <string>
#include <assert.h>
#include <time.h>

// GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>

// EXT
#include <armadillo>
#include "mpi.h"

using namespace std;
using namespace arma;

struct Moments {
  double mean;
  double std;
  double skewness;
};

class Stats {

public:
  
  Stats();
  ~Stats();

  static gsl_rng* rng;
  static void gsl_rng_init();
  static void gsl_rng_init(unsigned long int seed);
  
  static void xCovariance(gsl_matrix *r, gsl_matrix *m);
  static void xCovarianceTest(string inFile, string outFile);
  static Moments xMomentSummary(gsl_vector* data);
  static void printMoments(Moments moments);
  static double sum_vector(gsl_vector* d);
  static double cov(const gsl_vector *a, const gsl_vector *b);
  static double corr(const gsl_vector *a, const gsl_vector *b);
  static double autocorr(const gsl_vector *a, const int lag);
  static double sum(const gsl_vector *a);
  static double mean(const gsl_vector *a);

  //static double normpdf(const double x, const double sigma);
  static vec normpdf(const vec x, const double mu, const double sigma);
  static double tpdf(const double x, const double nu);
  static double tinv(const double p, const double nu);
  
  static mat randn(double mu, double sigma, int nRows, int nCols);
  static mat mvnrnd(vec mu, mat sigma, int nObs);
  static mat mvnrnd_2(vec mu, int nObs, mat A);
  static void mvnrnd_test();
  
  static mat gamrnd(double a, double b, int nRows, int nCols);

  static mat cholcov(mat Sigma);
  static mat cov2para(mat X, double &shrinkage);

  static vec empCDF(vec data, vec y);
  static void empCDF_fast(vec x0, vec &x1, vec &Fx);

  static vec diff(vec x);
  
  static void tester(int c);

  static void cleanupCDF(vec &x, vec &cdf, double tol);

  static vec ksdensity(vec xi, vec x);
  
};

#endif
