/*
  Author: Jimmie Goode
  Created: 2012-09-01
*/

#include "Stats.hpp"

gsl_rng* Stats::rng = gsl_rng_alloc(gsl_rng_taus);

//-------------------------------------------------------------------------------

Stats::Stats() {
}

//-------------------------------------------------------------------------------

Stats::~Stats() {
}


//-------------------------------------------------------------------------------
// initialize the GSL Random number generator

void Stats::gsl_rng_init() {
  
  unsigned long int seed = (long) time (NULL);
  gsl_rng_set(Stats::rng, seed);

  cout << "---> GSL RNG SEEDED AUTOMATICALLY TO " << seed << " <---\n";
}

//-------------------------------------------------------------------------------
// initialize the GSL Random number generator

void Stats::gsl_rng_init(unsigned long int seed) {
  
  gsl_rng_set(Stats::rng, seed);

  cout << "---> GSL RNG SEEDED MANUALLY TO " << seed << " <---\n";
}


//-------------------------------------------------------------------------------

// Estimate covariance
void Stats::xCovariance(gsl_matrix *r, gsl_matrix *m) {
  gsl_vector_view a, b;
  size_t i, j;
  
  for (i = 0; i < m->size2; i++) {
	for (j = i; j < m->size2; j++) {
	  a = gsl_matrix_column (m, i);
	  b = gsl_matrix_column (m, j);
	  double cov = gsl_stats_covariance(a.vector.data, a.vector.stride,
										b.vector.data, b.vector.stride, a.vector.size);
	  gsl_matrix_set (r, i, j, cov);
	}
  }
  // cov is symmetric so we copy the lower diagonal to the upper
  for (i = 1; i < m->size2; i++) {
	for (j = 0; j < i; j++) {
	  gsl_matrix_set(r, i, j, gsl_matrix_get(r, j, i));
	}
  }
}


//-------------------------------------------------------------------------------

// // test the covariance function
// void Stats::xCovarianceTest(string inFile, string outFile) {
//   gsl_matrix *mnData = IO::importCSVmatrix(inFile);
//   gsl_matrix *mnCov  = gsl_matrix_alloc(mnData->size2, mnData->size2);

//   gsl_matrix_set_all(mnCov, 0.0);
  
//   xCovariance(mnCov, mnData);
//   printf("\n------- Data --------\n");
//   IO::printMatrix(mnData);
  
//   printf("\n---- Cov Matrix -----\n");
//   IO::printMatrix(mnCov);

//   IO::exportGslMatrix(mnCov, outFile);

//   gsl_matrix_free(mnData);
//   gsl_matrix_free(mnCov);
// }

//-------------------------------------------------------------------------------

Moments Stats::xMomentSummary(gsl_vector* data) {

  // Sum the data
  int n = data->size;
  double sum = Stats::sum_vector(data);
  double mu = sum / (1.0*n);

  sum = 0.0;
  double sumSkew = 0.0;
  double diff;

  for (int i = 0; i < n; i++) {
	diff = gsl_vector_get(data,i) - mu;
	sum += pow(diff, 2);
	sumSkew += pow(diff, 3);
  }

  double sigma = sqrt(sum/(1.0 * n));
  double skew = (sumSkew/(1.0 * n)) / pow(sigma, 3);

  Moments mom;
  mom.mean     = mu;
  mom.std      = sigma;
  mom.skewness = skew;
  return mom;
}

//-------------------------------------------------------------------------------

void Stats::printMoments(Moments moments) {
  printf("Moments = [%g, %g, %g]\n", moments.mean, moments.std, moments.skewness);
}

//-------------------------------------------------------------------------------

// <TODO> Remove this and replace with sum(const gsl_vector*) below
double Stats::sum_vector(gsl_vector* d) {
  double s = 0.0;
  for (int i = 0; i < d->size; i++) {
	s += gsl_vector_get(d,i);
  }
  return s;
}

//-------------------------------------------------------------------------------

// covariance of two vectors
double Stats::cov(const gsl_vector *a, const gsl_vector *b) {
  return gsl_stats_covariance(a->data, a->stride, b->data, b->stride, a->size);
}

//-------------------------------------------------------------------------------

// correlation of two vectors
double Stats::corr(const gsl_vector *a, const gsl_vector *b) {
  return gsl_stats_correlation(a->data, a->stride, b->data, b->stride, a->size);
}

// sum of one vector
double Stats::sum(const gsl_vector *a) {
  double s = 0.0;
  for (int i = 0; i < a->size; i++) {
	s += gsl_vector_get(a, i);
  }
  return s;
}

//-------------------------------------------------------------------------------

// mean of vector
double Stats::mean(const gsl_vector *a) {
  return gsl_stats_mean(a->data, a->stride, a->size); //sum(a)/(1.0*(a->size));
}

//-------------------------------------------------------------------------------

// autocorrelation of a vector with a single lag
double Stats::autocorr(const gsl_vector *x, const int lag) {
  
  assert(lag >= 0);
  
  if (lag == 0)
	return 1.0;
  
  int n = x->size;

  gsl_vector *xx = gsl_vector_alloc(n-lag);
  gsl_vector *xl = gsl_vector_alloc(n-lag);  

  // 1 : end-lag
  for (int i = 0; i < n-lag; i++) {
	gsl_vector_set(xx, i, gsl_vector_get(x, i));
  }

  // lag : end
  for (int i = lag; i < n; i++) {
	gsl_vector_set(xl, i-lag, gsl_vector_get(x, i));
  }
  
  double nCorr = corr(xx, xl);
  gsl_vector_free(xx);
  gsl_vector_free(xl);

  return nCorr;
}

//-------------------------------------------------------------------------------

// double Stats::normpdf(const double x, const double sigma) {
//   return gsl_ran_gaussian_pdf(x, sigma);
// }

//-------------------------------------------------------------------------------

double Stats::tpdf(const double x, const double nu) {
  return gsl_ran_tdist_pdf(x, nu);
}

//-------------------------------------------------------------------------------

double Stats::tinv(const double p, const double nu) {
  return gsl_cdf_tdist_Pinv(p, nu);
}

//-------------------------------------------------------------------------------
// Gaussian random variate with mean mu and standard deviation sigma.

mat Stats::randn(double mu, double sigma, int nRows, int nCols) {

  mat samp(nRows, nCols);
  
  for (int i = 0; i < nRows; i++)
	for (int j = 0; j < nCols; j++)
	  samp(i,j) = mu + gsl_ran_gaussian(Stats::rng, sigma);

  return samp; 
}

//-------------------------------------------------------------------------------
// Multivariate Gaussian random variates with mean vector mu and covariance Sigma

mat Stats::mvnrnd(vec mu, mat sigma, int nObs) {

  if (mu.n_rows != sigma.n_cols)
	assert(0 && "Dimension Mismatch");

  return Stats::randn(0, 1, nObs, mu.n_rows)*chol(sigma) + repmat(mu.t(), nObs, 1);
}

//-------------------------------------------------------------------------------
// With Cholesky factor A = chol(sigma)

mat Stats::mvnrnd_2(vec mu, int nObs, mat A) {

  if (mu.n_rows != A.n_cols)
	assert(0 && "Dimension Mismatch");
 
  return Stats::randn(0, 1, nObs, mu.n_rows)*A + repmat(mu.t(), nObs, 1);
}


//-------------------------------------------------------------------------------

void Stats::mvnrnd_test() {

  vec mu(2);
  mat Sigma;

  mu(0) = 1.0;
  mu(1) = -1.0;
  
  Sigma << 0.9 << 0.4 << endr
        << 0.4 << 0.3 << endr;

  mat samp = Stats::mvnrnd(mu, Sigma, 10000);

  mu.print("mu = ");
  Sigma.print("Sigma = ");
  arma::mean(samp).print("mu_hat = ");
  arma::cov(samp).print("Sigma_hat = ");

}



//-------------------------------------------------------------------------------
// Gamma random deviates with shape parameter a and scale parameter b.

mat Stats::gamrnd(double a, double b, int nRows, int nCols) {

  mat samp(nRows,nCols);

  for (int i = 0; i < nRows; i++)
	for (int j = 0; j < nCols; j++)
	  samp(i,j) = gsl_ran_gamma(Stats::rng, a, b);
	
  return samp;
}

//------------------------------------------------------------------------------- 
// Cholesky-like decomposition for covariance matrix.
// T = CHOLCOV(SIGMA) computes T such that SIGMA = T'*T.  SIGMA must be
//   square, symmetric, and positive semi-definite.  If SIGMA is positive
//   definite, then T is the square, upper triangular Cholesky factor.
//
//   If SIGMA is not positive definite, T is computed from an eigenvalue
//   decomposition of SIGMA.  T is not necessarily triangular or square in
//   this case.  Any eigenvectors whose corresponding eigenvalue is close to
//   zero (within a small tolerance) are omitted.  If any remaining
//   eigenvalues are negative, T is empty.
//
// <TODO> Remove debugging and exporting of data.

mat Stats::cholcov(mat Sigma) {

  mat T;
  
  // test for being square and symmetric
  int n = Sigma.n_rows;
  int m = Sigma.n_cols;

  double tol = as_scalar(10*eps(arma::max(arma::abs(diagvec(Sigma)))));
  uvec idxBad = find((Sigma - Sigma.t()) >= tol);
  
  if (n == m && idxBad.n_rows == 0) {
	cout << "here\n";

	// Can get factors of the form Sigma==T'*T using the eigenvalue
	// decomposition of a symmetric matrix, so long as the matrix
	// is positive semi-definite.
	// [U,D] = eig(full((Sigma+Sigma')/2));
	
	vec vnD;
	mat U;

	// eigensystem: divide and conquer for larger matrices
	eig_sym(vnD, U, (Sigma + Sigma.t())/2.0, "standard");

	mat D = diagmat(vnD);

	D.save("/Users/jimmiegoode/Documents/Glimm/Verification/D0.csv", csv_ascii);
	U.save("/Users/jimmiegoode/Documents/Glimm/Verification/U0.csv", csv_ascii);

    // Pick eigenvector direction so max abs coordinate is positive
	// [ignore,maxind] = max(abs(U),[],1);

	//U.print("U = ");
	//D.print("D = ");

	cout << "here2\n";

	// find the locations of max(abs(U))
	mat absU = abs(U);
	vec maxind = zeros(m);

	cout << "here3\n";
	
	for (int c = 0; c < m; c++) {
	  double mx = absU(0,c);
	  int mx_idx = 0;
	  
	  for (int r = 1; r < m; r++) {
		if (absU(r,c) > mx) {
		  mx = absU(r,c);
		  mx_idx = r;
		}
	  }
	  maxind(c) = mx_idx;
	}

	cout << "here4\n";

	// set each eigenvector's direction so that max abs coordinate is positive
	for (int c = 0; c < m; c++) {
	  if (U(maxind(c),c) < 0.0) {
		for (int r = 0; r < m; r++)
		  U(r,c) = -U(r,c);
	  }
	}

	U.save("/Users/jimmiegoode/Documents/Glimm/Verification/U.csv", csv_ascii);
	
	tol = eps(max(vnD)) * vnD.n_rows;
	uvec t = find(abs(vnD) > tol);

	//t.print("t = ");

	cout << "here5, t.n_rows = " << t.n_rows << ", tol = " << tol << "\n";
	
	vec vnD2 = zeros(t.n_rows);
	mat U2 = zeros(m, t.n_rows);
    int p = 0;

	cout << "here6\n";
	
	for (int i = 0; i < t.n_rows; i++) {
	  vnD2(i) = vnD(t(i));
	  U2.col(i) = U.col(t(i));
	  
	  if (vnD2(i) < 0.0)
		p++;
	}

	cout << "here7, p = " << p << ", U2.n_rows = " << U2.n_rows << ", U2.n_cols = " << U2.n_cols << "\n"; 
	
	if (p == 0) {
	  T = diagmat(sqrt(vnD2)) * U2.t();
	  T.save("/Users/jimmiegoode/Documents/Glimm/Verification/T.csv", csv_ascii);
	}

	cout << "here8\n";
	
	//vnD2.print("vnD2 = ");
	//T.print("T = ");	
  }

  return T;
}

//-------------------------------------------------------------------------------
// Ledoit & Wolf two parameter shrinkage

mat Stats::cov2para(mat X, double &shrinkage) {

  int t = X.n_rows;
  int n = X.n_cols;

  vec meanx = arma::mean(X,0).t();  
  X = X - repmat(meanx.t(), t, 1);

  mat sample = (1.0/t) * X.t() * X;
  
  double meanvar = as_scalar(arma::mean(diagvec(sample)));
  double meancov = 0.0;
  
  for (int r = 0; r < n; r++)
	for (int c = 0; c < n; c++)
	  if (r != c)
		meancov += sample(r,c);
  
  meancov = meancov/n/(n-1);

  mat mnEye(n,n);
  mnEye.eye();
  
  mat prior = meanvar*mnEye + meancov*(ones(n,n) - diagmat(ones(n)));

  double c = std::pow(arma::norm(sample - prior, "fro"), 2.0);
  mat Y = pow(X, 2.0);
  
  double p = (1.0/t)*arma::sum(arma::sum(Y.t()*Y)) - arma::sum(arma::sum(arma::pow(sample, 2.0)));
  double r = 0.0;
  
  double rdiag = (1.0/n)*arma::sum(arma::sum(arma::cov(Y)));
  double roff = 0.0;

  //cout << "----> Stats:: arma::sum(arma::sum(Y.t()*Y)) = " << arma::sum(arma::sum(Y.t()*Y)) << endl;
  //cout << "----> Stats:: p = " << p << endl;

  //cout << "----> Stats:: arma::sum(arma::sum(arma::cov(Y))) = "
  //	   << arma::sum(arma::sum(arma::cov(Y))) << endl;

  //cout << "----> Stats:: rDiag = " << rdiag << endl;

  Y.print("----> Stats:: Y = ");
  arma::cov(Y).print("----> Stats:: cov(Y) = ");
  arma::sum(Y, 0).t().print("-----> Stats:: sum(Y) = ");

  // <TODO> This as missing before, test now that it's here!
  r = rdiag + roff;
  
  double k = (p - r)/c;
  shrinkage = std::max(0.0, std::min(1.0, k/t));
  
  return shrinkage*prior + (1-shrinkage)*sample;
  
}

/*
//-------------------------------------------------------------------------------
// Ledoit & Wolf diagonal Shrinkage

void Stats::shrink_diag(mat x, mat &sigma, double &shrinkage) {

shrinkage = 0.0;
  
int t = x.n_rows;
int n = x.n_cols;
double td = (double) t;

// de-mean returns
vec meanx = mean(x);
x = x - repmat(meanx.t(), t, 1);

// compute sample covariance matrix
mat sample = (1/t)*(x.t()*x);

// compute prior
mat prior = diagmat(diagvec(sample));

// compute shrinkage parameters
  
// what we call p 
mat y = arma::pow(x,2);
  
mat phiMat = y.t() * (y/td) - 2*(x.t()*x) % sample/td + arma::pow(sample,2);
double phi = sum(sum(phiMat,0));
  
// what we call r
double rho = sum(diagvec(phiMat));
  
// what we call c
double froNorm = norm(sample-prior,"fro");
double gamma = froNorm*froNorm;

// compute shrinkage constant
double kappa = (phi-rho)/gamma;
shrinkage = std::max(0.0, min(1.0, kappa/td));
  
// compute shrinkage estimator
sigma = shrinkage*prior + (1-shrinkage)*sample;

cout << "shrinkage = " << shrinkage << endl;
sigma.print("sigma = ");
  
}
*/

//-------------------------------------------------------------------------------
// Linear interpolation given known points (x0,y0) and (x1,y1).

double interp_linear(double x0, double y0, double x1, double y1, double x) {
  
  return y0 + (x - x0)*((y1 - y0)/(x1 - x0));
}

//-------------------------------------------------------------------------------
// Empirical CDF of data with interpolation for points in y.

vec Stats::empCDF(vec data, vec y) {

  int n = data.n_rows;

  // x values
  vec x(n+1);
  x(0) = min(x);
  x(span(1,n)) = sort(data);
  
  // empirical CDF
  vec F(n+1);
  F(0) = 0.0;
  F(span(1,n)) = cumsum(ones(n,1) * (1./((double) n)));

  vec Fx(y.n_rows);

  for (int i = 0; i < y.n_rows; i++) {

	// expensive find
	uvec idx = find(x < y(i), 1, "last");

	if (idx.n_rows == 0) {
	  Fx(i) = F(0);
	} else if (idx(idx.n_rows-1) == x.n_rows-1) {
	  Fx(i) = F(F.n_rows-1);
	} else {
	  int ind = idx(idx.n_rows-1);
	  Fx(i) = interp_linear(x(ind), F(ind), x(ind+1), F(ind+1), y(i));
	}
  }
  
  return Fx;
}

//-------------------------------------------------------------------------------
// Much faster than empCDF if you need the CDF at every input point

void Stats::empCDF_fast(vec x0, vec &x1, vec &Fx) {

  assert(x1.n_rows == Fx.n_rows);
  assert(x0.n_rows == x1.n_rows);
  
  int n = x0.n_rows;

  // sorted indices for x0
  uvec idxs = sort_index(x0);

  // assign sorted data
  x1 = x0.elem(idxs);
  
  // sort the sorted indices to get indices of F at orignal points in x0
  uvec idxss = sort_index(idxs);

  // empirical CDF
  vec F = cumsum(ones(n,1) * (1./((double) n)));

  for (int i = 0; i < n; i++) {
	Fx(i) = F(idxss(i));
  }
}

//-------------------------------------------------------------------------------
// First differences

vec Stats::diff(vec x) {
  vec out(x.n_rows - 1);

  for (int i = 1; i < x.n_rows; i++)
	out(i-1) = x(i) - x(i-1);

  return out;
}


//-------------------------------------------------------------------------------
// Put CDF into (0,1) domain and force it to be strictly monotonic.
// Make a total of maxIters attempts before quitting.

void Stats::cleanupCDF(vec &x, vec &cdf, double tol) {

  assert(x.n_rows == cdf.n_rows);

  if (cdf.n_rows <= 1) {
	cout << "Stats.cpp >> length(cdf) = 1" << endl;
	return;
  }
  
  uvec idx = find(cdf > tol);
  x = x.elem(idx);
  cdf = cdf.elem(idx);
  if (cdf.n_rows == 0) return;
  
  idx = find(cdf < 1.-tol);
  x = x.elem(idx);
  cdf = cdf.elem(idx);
  if (cdf.n_rows == 0) return;
  
  idx = find(Stats::diff(cdf) != 0.);
  x = x.elem(idx);
  cdf = cdf.elem(idx);
  if (cdf.n_rows == 0) return;
 
  int maxIters = 5;
  int i;
 
  for(i = 0; i < maxIters; i++) {
	uvec idxBad = find(Stats::diff(cdf) <= 0.);

	if (idxBad.n_rows == 0)
	  break;

	idx = find(Stats::diff(cdf) > 0.);
	x = x.elem(idx);
	cdf = cdf.elem(idx);
  }

  if (i == maxIters) {
	cout << "FAILURE: to cleanupCDF" << endl;
  }//  else {
  // 	cout << "cleanupCDF took " << i << " iters." << endl;
  // }
  
  assert(x.n_rows == cdf.n_rows);
}


//-------------------------------------------------------------------------------
// Gaussian PDF

vec Stats::normpdf(const vec x, const double mu, const double sigma) {
  
  return (1.0/(sigma*sqrt(2*M_PI))) * arma::exp(-arma::pow(x-mu, 2.0)/(2.0*sigma*sigma));
  
}

//-------------------------------------------------------------------------------
// Kernel smoothed PDF with Gaussian kernel.
// Bandwidth h is chosen via heuristic and is only applicable for Gaussian xi.
// INPUT:
//     xi - vector of observations
//     x  - vector of points to evaluate kernel density at
// OUPUT:
//     f  - vector of PDF evaluated at each x

vec Stats::ksdensity(vec xi, vec x) {

  int iN = xi.n_rows;
  double dN = (double) iN;

  double med = arma::median(xi);
  double sig = arma::median(arma::abs(xi-med)) / 0.6745;
  double h   = sig * pow(4.0/(3.0*dN), 0.2);

  vec f(x.n_rows);

  for (int i = 0; i < x.n_rows; i++) 
	f(i) = (1.0/(dN*h)) * arma::sum(Stats::normpdf((x(i)-xi)/h, 0, 1));
  
  return f; 
}


//-------------------------------------------------------------------------------

void Stats::tester(int c) {

  cout << "\n---> Stat Tester: Case " << c << " <---\n";
	
  switch (c) {
  case 1:
	{
	  mat Sigma;
  
	  Sigma << 1 << 2 << endr
			<< 2 << 4 << endr;
  
	  Sigma << 2 << 1 << 1 << 2 << endr
			<< 1 << 2 << 1 << 2 << endr
			<< 1 << 1 << 2 << 2 << endr
			<< 2 << 2 << 2 << 3 << endr;

	  Sigma.print("Sigma = ");

	  Sigma.load("/Users/jimmiegoode/Documents/Glimm/Verification/mnCov.csv", csv_ascii);
	  Sigma.load("/Users/jimmiegoode/Documents/Glimm/Verification/mnU_ML.csv", csv_ascii);
  
	  mat T = Stats::cholcov(Sigma);

	  double shrink = 0.0;
	  Sigma = Stats::cov2para(Sigma, shrink);

	  cout << "shrinkage = " << shrink << endl;
	  Sigma.save("/Users/jimmiegoode/Documents/Glimm/Verification/mnCov2.csv", csv_ascii);
  
	  T.print("T = ");
	}
	break;

  case 2:
	{
	  vec data;
	  data.load("/Users/jimmiegoode/Documents/Glimm/data/scalingSeries/1.csv", csv_ascii);

	  vec x(3);

	  x(0) = -1.;
	  x(1) =  0.;
	  x(2) =  1.;

	  x.print("x = ");
	  
	  vec Fx = Stats::empCDF(data, data);

	  Fx.print("Fx = ");
	  Fx.save("/Users/jimmiegoode/Documents/Glimm/skewt_vs_ASSG/e_Fx.csv", csv_ascii);
 
	}
	break;
	
  case 3:
	{
	  vec data;
	  data.load("/Users/jimmiegoode/Documents/Glimm/data/scalingSeries/1.csv", csv_ascii);

	  int n = data.n_rows;
	  
	  vec x(n);
	  vec F(n);
	  Stats::empCDF_fast(data, x, F);

	  x.save("/Users/jimmiegoode/Documents/Glimm/Toolbox/empCDF/testing/e_x.csv", csv_ascii);	  
	  F.save("/Users/jimmiegoode/Documents/Glimm/Toolbox/empCDF/testing/e_Fx.csv", csv_ascii);
 
	}
	break;

  case 4:
	{
	  vec cdfi, x;
	  cdfi.load("/Users/jimmiegoode/Documents/Glimm/src/e_cdfi.csv", csv_ascii);
	  x.load("/Users/jimmiegoode/Documents/Glimm/src/e_x.csv", csv_ascii);

	  cout << x.n_rows << ", " << cdfi.n_rows << endl;
	  cout << arma::mean(cdfi) << endl;
	  
	  Stats::cleanupCDF(x, cdfi, 0.0);

	  cout << x.n_rows << ", " << cdfi.n_rows << endl;

	  cdfi.save("/Users/jimmiegoode/Documents/Glimm/src/e_cdfi3.csv", csv_ascii);
	  x.save("/Users/jimmiegoode/Documents/Glimm/src/e_x3.csv", csv_ascii);
	  
	}
	break;

  case 5:
	{
	  cout << "Testing KS Density" << endl;
	  
	  vec xi, x;

	  xi.load("/Users/jimmiegoode/Documents/Glimm/LookupTables/xi.csv", csv_ascii);
	  x.load("/Users/jimmiegoode/Documents/Glimm/LookupTables/x.csv", csv_ascii);

	  vec f = Stats::ksdensity(xi, x);

	  f.save("/Users/jimmiegoode/Documents/Glimm/LookupTables/e_f.csv", csv_ascii);

	}
	break;
	
  } // end of switch
  
} // end of method




//-------------------------------------------------------------------------------

// int main(int argc, char** argv) {
  
//   if (argc < 2)
// 	assert(0 && "Invalid Arguments");
  
//   int testCase = atoi(argv[1]);

//   Stats::tester(testCase);
  
//   return 0;
// }
