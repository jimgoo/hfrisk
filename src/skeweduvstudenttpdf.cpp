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
yhuang0824@gmail.com
*/

#include "skeweduvstudenttpdf.h"


void skeweduvstudenttpdf(gsl_vector* logPDF,const gsl_vector* X, 
											const double gamma, 
											const double mu,
											const double df, 
											const double sigma,
											const int logIndicator)
{
  
	int nPoints = X->size;
	gsl_vector *Rho=gsl_vector_alloc(nPoints);
	gsl_vector *bessel_arg = gsl_vector_alloc(nPoints);
	for ( int i =0; i < nPoints; i++)
	{
            //
		gsl_vector_set(Rho, i, pow( (gsl_vector_get(X,i) - mu )/sigma, 2  ) );
	}
	
	if (gamma ==0) 
	{  
		for ( int i =0; i < nPoints;  i++)
		{ 
			gsl_vector_set(logPDF, i, log (gammafun(  (df+1)/2.0  )  ) - log (gammafun(  df/2.0  )  )
									 - (1/2) *log(pi*df ) - log(sigma) - ( (df + 1)/2 ) 
									 * log(1+gsl_vector_get(Rho,i)/df ) );
		}
	}
	else 
	{
	for (int jj =0; jj<nPoints;jj++) 
	{
		gsl_vector_set(bessel_arg, jj, pow( (df + gsl_vector_get(Rho, jj) )* (pow(gamma,2))/(pow(sigma,2))  ,0.5) );
		gsl_vector_set(logPDF, jj, (log(2) + ( df / 2)*log(pi* df) - log(gammafun( df/2 ))-log(sigma)
								  + log(exp(  gsl_vector_get(bessel_arg, jj) )
								  *boost::math::cyl_bessel_k( ((df+1)/2),gsl_vector_get(bessel_arg, jj)  )  ) 
								  -gsl_vector_get(bessel_arg, jj) + ( ( gsl_vector_get(X, jj) - mu)*gamma / pow(sigma,2))
								  - ((df+1)/2) * log((( 1+ gsl_vector_get(Rho, jj) / df) * (2*pi*df) ) / gsl_vector_get(bessel_arg, jj))));
	}
  }
  
  if (logIndicator == 0)
  {
	for ( int i =0; i < nPoints; i++)
	{
		gsl_vector_set(logPDF, i, exp( gsl_vector_get(logPDF,i))  );
	}
  }
  gsl_vector_free(Rho);
  gsl_vector_free(bessel_arg);
}
