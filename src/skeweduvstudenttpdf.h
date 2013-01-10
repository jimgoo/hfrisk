/*
skeweduvstudenttpdf.cpp:
input:
		X: a vector of univariate distribution.
		gamma,mu,df,sigma: all scalars;
		logIndicator: 1 for log of pdf; 0 for pdf.
output: 
		a vector of PDFs.
usage:
		gsl_vector *logPDF=gsl_vector_alloc(X-size);
Dec 2011, Yan Huang. 
*/

#ifndef SKEWEDUVSTUDENTTPDF_H
#define SKEWEDUVSTUDENTTPDF_H

#include <iostream>
#include <vector>
#include <math.h>
//#include <omp.h>
#include <time.h>
// path might need to change
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include "gamma.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define pi 3.1415926535897
using namespace std;

/*
void skeweduvstudenttpdf(gsl_vector* logPDF,const gsl_vector* X,
						 const F<double>& gamma,
						 const F<double>& mu,
						 const F<double>& df,
						 const F<double>& sigma,
						 const int logIndicator);
*/

void skeweduvstudenttpdf(gsl_vector* logPDF,	const gsl_vector* X, 
											const double gamma, 
											const double mu,
											const double df, 
											const double sigma,
											const int logIndicator);

#endif
