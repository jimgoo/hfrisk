/*
  Author: Jimmie Goode
  Created: 2012-10-10
*/


#include "mystable.hpp"

//-------------------------------------------------------------------------------
// globals

// int iters;
// int fxcalls;
int doExport = 0;
int doPrintBounds = 0;
int doPrintLLF = 0;

cx_double dt(1.0, 0.0);
cx_double onei(0.0, 1.0);

//-------------------------------------------------------------------------------

string mystable::arr2str(const unsigned n, const double* x) {
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(5) << "[";  
  for (int i = 0; i < n; i++) {
	ss << x[i];
	if (i != n-1) ss << ", ";
  }
  ss << "]";
  return ss.str();
}

//-------------------------------------------------------------------------------
//

string mystable::vec2str(const unsigned n, const vector<double> x) {
  double xarr[n];
  std::copy(x.begin(),x.end(),xarr);
  return arr2str(n, xarr);
}

//-------------------------------------------------------------------------------
// Gamma function

double xGamma(double x) {
  return gsl_sf_gamma(x);
}

//-------------------------------------------------------------------------------

double xSign(double x) {
  if (x < 0.0)
	return -1.0;
  else if (x > 0.0)
	return 1.0;
  else if (x == 0.0)
	return 0.0;
  else {
	return 0.0;
  }
}

//-------------------------------------------------------------------------------

vec xSignVec(vec x) {
  vec out(x.n_rows);

  for (int i = 0; i < x.n_rows; i++)
	out(i) = xSign(x(i));

  return out;
}

//-------------------------------------------------------------------------------
// CTS CHF

cx_vec chf_CTS(vec u, double alpha, double C, double lmp, double lmm, double mm) {

  double m = mm - xGamma(1.0 - alpha) * C * (std::pow(lmp, alpha-1.0)
													   - std::pow(lmm, alpha-1.0));
  
  cx_vec phi = arma::exp(dt * (onei * u * m + C * xGamma(-alpha) * (
					  pow(lmp - onei * u, alpha) - pow(lmp,alpha)
					  + pow(lmm + onei * u, alpha) - pow(lmm, alpha))));
  return phi;
}


//-------------------------------------------------------------------------------
// stdCTS CHF

cx_vec chf_stdCTS(vec u, double alpha, double lmp, double lmm) {

  double C = 1.0/xGamma(2.0-alpha)/(pow(lmp, alpha-2.0) + pow(lmm, alpha-2.0));

  cx_vec phi = chf_CTS(u, alpha, C, lmp, lmm, 0.0);

  return phi;
}


//-------------------------------------------------------------------------------
// NTS CHF

cx_vec chf_NTS(vec u, double alpha, double C, double lambda, double beta, double m) {

  double CC = C*sqrt(M_PI)*xGamma(-alpha/2.0)*pow(2, -alpha/2.0 - 0.5);

  cx_vec phi = arma::exp(dt*onei*u*m
						 + dt*CC*(
								  pow(pow(lambda, 2.0) - pow(beta + onei*u, 2.0), alpha/2.0)
								  -pow(lambda*lambda - beta*beta, alpha/2.0)
								  + onei*u*alpha*beta*pow(lambda*lambda - beta*beta, alpha/2.0 - 1.0)));
  return phi;
}

//-------------------------------------------------------------------------------
// stdNTS CHF

cx_vec chf_stdNTS(vec u, double alpha, double lambda, double beta) {

  double al = alpha;
  double lsq = lambda*lambda;
  double bsq = beta*beta;

  double c = pow(2.0, (al+1.0)/2.0)/(sqrt(M_PI)*al*xGamma(-al/2.)*pow(lsq-bsq, al/2.-2.)*(al*bsq-lsq-bsq));

  cx_vec phi = chf_NTS(u, alpha, c, lambda, beta, 0);

  return phi;
}


//-------------------------------------------------------------------------------
// alpha-stable CHF

cx_vec chf_AS(vec u, double alpha, double beta, double sigma, double mu) {

  if (alpha == 1.0) {
	alpha = 1.0 + 1.0e-5;
  }
  
  cx_vec phi = arma::exp(onei*u*mu - pow(sigma,alpha)*pow(abs(u),alpha)
						 %(1. - onei*beta*xSignVec(u)*tan(M_PI*alpha/2.)));
 
  return phi;
}

//-------------------------------------------------------------------------------
// standardized alpha stable CHF

cx_vec chf_stdAS(vec u, double alpha, double beta) {
  return chf_AS(u, alpha, beta, 1.0, 0.0);  
}


//-------------------------------------------------------------------------------
// Generalized Hyperbolic CHF
// SOURCE: Page 125, Financial Modelling With Jump Processes Cont and Tankov, 2004, 1st. ed)

//
// THE PROBLEM: complex bessel function argument not supported by Boost.
//

// cx_vec besselk(double order, cx_vec u) {

//   cx_vec out(u.n_rows);

//   complex<double> val(2.,2.);
//   //boost::math::complex<double> val(2.,2.);
//   //boost::math::cyl_bessel_k(order, val);
//   //boost::complex z3(1., 2.);
  
//   // for (int i = 0; i < u.n_rows; i++) {
//   //   complex<double> val(arma::real(u(i)), arma::imag(u(i)));
//   // 	out(i) = boost::math::cyl_bessel_k(order, val);
//   // }
//   return out;
// }


// cx_vec chf_GH(vec u, double lambda, double alpha, double beta, double delta, double mu) {

//   cx_vec tmp = arma::pow(beta + onei*u, 2.);
//   //(besselk(lambda, delta * arma::pow(lambda*lambda - tmp, 0.5)) /
  
//   cx_vec phi = arma::exp(onei*u*mu) *
// 	arma::pow((alpha*alpha - beta*beta)/(alpha*alpha - tmp), lambda/2.) *
// 	(besselk(lambda, tmp))/
// 	  boost::math::cyl_bessel_k(lambda, delta * std::pow(alpha*alpha - beta*beta, 0.5)));

//   return phi;
// }


//-------------------------------------------------------------------------------
//

void pdf_FFT(const int narg, const int nparam,
			double arg[], const double param[],
			 mystable::dist dist, vec &f, vec &x, vec &pdf) {

  double h = mystable::gridStep; // 0.01
  double N = mystable::gridSize; //pow(2.0, 13);
  const int iN = (int) N;

  // check that sizes are correct
  assert(iN == x.n_rows && iN == pdf.n_rows);
  assert(narg == f.n_rows);
  
  double s = 1.0/(h*N);
  vec t1 = linspace(1,N,N);
  vec t2 = 2.0*M_PI*(t1 - 1.0 - N/2.0)*s;

  cx_vec cfvalues(iN);

  // set named parameters and get CHF values
  switch (dist) {
  case mystable::stdCTS:
	{
	  double alpha = param[0];
	  double lmp   = param[1];
	  double lmm   = param[2];
	  cfvalues = chf_stdCTS(t2, alpha, lmp, lmm);
	}
	break;
  case mystable::stdNTS:
	{
	  double alpha  = param[0];
	  double lambda = param[1];
	  double beta   = param[2];
	  cfvalues = chf_stdNTS(t2, alpha, lambda, beta);
	}
	break;
  case mystable::stdAS:
	{
	  double alpha = param[0];
	  double beta  = param[1];
	  cfvalues = chf_stdAS(t2, alpha, beta);
	}
	break;
  case mystable::symAS:
	{
	  double alpha = param[0];
	  double sigma = param[1];
	  double mu    = param[2];
	  cfvalues = chf_AS(t2, alpha, 0.0, sigma, mu);
	}
	break;
  case mystable::AS:
	{
	  double alpha = param[0];
	  double beta  = param[1];
	  double sigma = param[2];
	  double mu    = param[3];

	  /*
	  cout << "----> alpha = " << alpha << ", " 
		   << "beta = " << beta << ", "
		   << "sigma = " << sigma << ", "
		   << "mu = " << mu << endl;
	  */
	  
	  cfvalues = chf_AS(t2, alpha, beta, sigma, mu);
	}
	break;
  default:
	assert(0 && "not implemented");
  }

  //
  // FFT
  //
  
  fftw_plan p;
  fftw_complex in[iN];
  fftw_complex out[iN];
  
  for (int i = 0; i < iN; i++) {
	complex<double> val = pow(-1.0, t1(i) - 1.0) * cfvalues(i);
	in[i][0] = real(val);
	in[i][1] = imag(val);
  }
  
  //p = fftw_plan_dft_c2r_1d(N, reinterpret_cast<fftw_complex*>(x1), out, FFTW_ESTIMATE);  
  //p = fftw_plan_dft_1d(iN, reinterpret_cast<fftw_complex*>(in), out, FFTW_FORWARD, FFTW_ESTIMATE);
  p = fftw_plan_dft_1d(iN, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);

  for (int i = 0; i < iN; i++) {
	double tmp0 = t1(i) - 1.0 - N/2.0;
	double tmp1 = out[i][0]*s*pow(-1.0, tmp0);
	pdf(i) = max(0.0, tmp1);
	x(i) = tmp0*h;
  }

  //
  // Interpolation
  //

  //gsl_interp_accel *acc = gsl_interp_accel_alloc();
  //gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, iN);
  //gsl_spline_init(spline, x.memptr(), pdf.memptr(), iN);  

  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, iN);
  gsl_interp_init(interp, x.memptr(), pdf.memptr(), iN);

  //vec f(narg);
  double xmin = x[0];
  double xmax = x[iN-1];

  //cout << "extrema: " << xmin << ", " << xmax << endl;
  
  for (int i = 0; i < narg; i++) {
	// set the probability of anything outside the grid to a tiny value (so that log PDF is still finite)
	if (arg[i] >= xmax || arg[i] <= xmin) {
	  f(i) = 1.0e-20;
	} else {
	  f(i) = gsl_interp_eval(interp, x.memptr(), pdf.memptr(), arg[i], acc);
	}
  }

  gsl_interp_accel_free(acc);
  gsl_interp_free(interp);
  fftw_destroy_plan(p);
  //fftw_free(in);
  //fftw_free(out);

  if (doExport) {
	cout << "saving data..." << endl;
	string path = "/Users/jimmiegoode/Documents/Papers/paper2/Toolbox_Project/";
	pdf.save(path + "pdf", csv_ascii);
	x.save(path + "x", csv_ascii);
	f.save(path + "f", csv_ascii);
  }
  
}

//-------------------------------------------------------------------------------

vec mystable::cdf_FFT(vec arg, const int nparam, const double param[],
			mystable::dist dist) {

  int iN = (int) mystable::gridSize;
  vec f(1);
  vec x(iN);
  vec pdf(iN);

  double xpdf[1];
  xpdf[0] = 0.;

  // call PDF to get grid
  pdf_FFT(1, nparam, xpdf, param, dist, f, x, pdf);

  // multipy by the step size in x
  vec cdfi = cumsum(pdf) * (x(1) - x(0));

  // check for values outside of grid
  double xmax = max(x);
  double xmin = min(x);
  
  if (min(arg) < xmin || max(arg) > xmax) {
  	uvec in = find(arg > xmax);
  	arg.elem(in) = xmax * ones(in.n_rows);
  	in = find(arg < xmin);
  	arg.elem(in) = xmin * ones(in.n_rows);
  }
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, iN);
  gsl_interp_init(interp, x.memptr(), cdfi.memptr(), iN);

  int narg = arg.n_rows;
  vec cdf(narg);
  
  for (int i = 0; i < narg; i++) {
	cdf(i) = gsl_interp_eval(interp, x.memptr(), cdfi.memptr(), arg(i), acc);
  }
  gsl_interp_accel_free(acc);
  gsl_interp_free(interp);

  return cdf;
}

//-------------------------------------------------------------------------------

vec mystable::inv_FFT(vec pvalues, const int nparam, const double param[],
			mystable::dist dist) {
  
  int iN = (int) mystable::gridSize;
  vec f(1);
  vec x(iN);
  vec pdf(iN);

  double xpdf[1];
  xpdf[0] = 0.;

  // call PDF to get grid

  pdf_FFT(1, nparam, xpdf, param, dist, f, x, pdf);
  
  /*
  if (dist == mystable::stdAS) {

	cout << "here\n";
	
	double paramAS[4];
	paramAS[0] = param[0];
	paramAS[1] = param[1];
	paramAS[2] = 1.0;
	paramAS[3] = param[1]*std::tan(M_PI); // mu + sigma^2*beta*tan(pi)
	
	pdf_FFT(1, 4, xpdf, paramAS, mystable::AS, f, x, pdf);
	
  } else {
	  pdf_FFT(1, nparam, xpdf, param, dist, f, x, pdf);
  }
  */

  // multipy by the step size in x
  vec cdfi = cumsum(pdf) * (x(1) - x(0));

  /*
  cout << "Saving...\n";
  cdfi.save("e_cdf_b.csv", csv_ascii);
  x.save("e_x_b.csv", csv_ascii);
  cout << "size: " <<  x.n_rows << ", " << cdfi.n_rows << endl;
  */
  
   // tail cropping tolerance
  double tol = 1.0e-6;
  
   // force CDF to be monotonic
  Stats::cleanupCDF(x, cdfi, tol);

  /*
  cout << "Saving after..\n";
  cdfi.save("e_cdf_a.csv", csv_ascii);
  x.save("e_x_a.csv", csv_ascii);
  cout << "size: " <<  x.n_rows << ", " << cdfi.n_rows << endl;
  */
  
  iN = cdfi.n_rows;
  int narg = pvalues.n_rows;
  vec q(narg);

  // If CDF was all bad (reduced to empty vector), set quantiles to 0 and exit.
  if (iN == 0) {
	cout << "----> length(iN) == 0\n";
	q.zeros();
	return q;
  }
  
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, iN);
  gsl_interp_init(interp, cdfi.memptr(), x.memptr(), iN);

  double xMax = max(x);
  double xMin = min(x);
  double cdfMin = min(cdfi);
  double cdfMax = max(cdfi);
  
  for (int i = 0; i < narg; i++) {
   
	if (pvalues(i) >= cdfMax)
	  q(i) = xMax;
	else if (pvalues(i) <= cdfMin)
	  q(i) = xMin;
	else
	  q(i) = gsl_interp_eval(interp, cdfi.memptr(), x.memptr(), pvalues(i), acc);

	if (gsl_finite(q(i)) != 1)
	  cout << "F_inv(" << pvalues(i) << ") = " << q(i) 
		   << ", max(cdf) = " << cdfMin << ", min(cdf) = " << cdfMax << endl;

  }
  
  gsl_interp_accel_free(acc);
  gsl_interp_free(interp);

  return q;
}


//-------------------------------------------------------------------------------

double checkNegLLF(double negLLF) {
  if (gsl_finite(negLLF) == 0) {
	negLLF = 1.0e20;
  }
  return negLLF;
}

//-------------------------------------------------------------------------------
//

typedef struct {
  int ny;
  int np;
  double *y;
  mystable::dist dist;
  int fxcalls;
  int iters;
} fit_struct;

//-------------------------------------------------------------------------------
// Called to get negative log likelihood

double negLLF_stable(unsigned n, const double *x, double *grad, void *data) {
  
  fit_struct *s = (fit_struct *) data;
  s->fxcalls++;

  
  //vec pdf = pdf_FFT(s->ny, s->np, s->y, x, s->dist);
  vec pdf(s->ny);
  vec nil1((int) mystable::gridSize);
  vec nil2((int) mystable::gridSize);
  pdf_FFT(s->ny, s->np, s->y, x, s->dist, pdf, nil1, nil2);
  
  double negLLF = -1.0*as_scalar(sum(log(pdf)));

  if (doPrintLLF) {
	cout << mystable::arr2str(n,x) << ", negLLF = " << negLLF << " --> " << checkNegLLF(negLLF) << endl;
  }
  
  return checkNegLLF(negLLF);
}



//-------------------------------------------------------------------------------
// NTS abs(beta) < lambda constraint

double constraint_NTS(unsigned n, const double *x, double *grad, void *data) {
  
    if (grad) {
	  assert(0 && "shouldn't be a gradient");
    }

	//double alpha  = x[0];
	double lambda = x[1];
	double beta   = x[2];
	
    return abs(beta) - lambda;
 }


//-------------------------------------------------------------------------------
// returns the number of parameters for each distribution

int mystable::getParCount(mystable::dist dist) {

  switch (dist) {
  case mystable::stdCTS:
	return 3;
  case mystable::stdNTS:
	return 3;
  case mystable::stdAS:
	return 2;
  case mystable::symAS:
	return 3;
  case mystable::AS:
	return 4;
  default:
	assert(0 && "bad dist");
  }
}


//-------------------------------------------------------------------------------

string mystable::dist2str(mystable::dist dist) {

  switch (dist) {
  case mystable::stdCTS:
	return "stdCTS";
  case mystable::stdNTS:
	return "stdNTS";
  case mystable::stdAS:
	return "stdAS";
  case mystable::symAS:
	return "symAS";
  case mystable::AS:
	return "AS";
  default:
	assert(0 && "bad dist");
  }
}


//-------------------------------------------------------------------------------
// Initial estimates and box constraints for MLE

void getInitialPars(const int n, vector<double> &x,
					vector<double> &lb, vector<double> &ub, mystable::dist dist) {

  double tol = 1.0e-6;
  
  switch (dist) {
  case mystable::stdCTS:
	{	  
	  x[0] = 1.5;  // alpha
	  x[1] = 0.05; // lam_pos
	  x[2] = 0.05; // lam_neg

	  lb[0] = tol;
	  lb[1] = tol;
	  lb[2] = tol;

	  ub[0] = 2.0 - tol;
	  ub[1] = 10.0;
	  ub[2] = 10.0;
	}
	break;
  case mystable::stdNTS:
	{ 
	  x[0] = 1.8;  // alpha
	  x[1] = 1.0;  // lambda
	  x[2] = -0.2; // beta

	  lb[0] = tol;
	  lb[1] = tol;
	  lb[2] = -1.0 + tol;

	  ub[0] = 2.0 - tol;
	  ub[1] = 10.0;
	  ub[2] = 1.0 - tol;
	}
	break;
  case mystable::stdAS:
	{
	  // This sets alpha > 1.0 for use with GARCH innovations.
	  
	  x[0] =  1.5; // alpha
	  x[1] = -0.1; // beta

	  lb[0] =  1.0 + tol;
	  lb[1] = -1.0 + tol;

	  ub[0] = 2.0 - tol;
	  ub[1] = 1.0 - tol;	  
	}
	break;
  case mystable::symAS:
	{
	  x[0] = 1.5; // alpha
	  x[1] = 1.0; // sigma
	  x[2] = 0.0; // mu

	  lb[0] = 1.0 + tol;
	  lb[1] = tol;
	  lb[2] = -HUGE_VAL;

	  ub[0] = 2.0 - tol;
	  ub[1] = HUGE_VAL;
	  ub[2] = HUGE_VAL;
	}
	break;
	case mystable::AS:
	  {
		x[0] = 1.5; // alpha
		x[1] = 0.0; // beta
		x[2] = 1.0; // sigma
		x[3] = 0.0; // mu

		lb[0] =  0.0 + tol;
		lb[1] = -1.0 + tol;
		lb[2] = -HUGE_VAL;
		lb[3] = -HUGE_VAL;

		ub[0] = 2.0 - tol;
		ub[1] = 1.0 + tol;
		ub[2] = HUGE_VAL;
		ub[3] = HUGE_VAL;
	  }
	  break;
  default:
	assert(0 && "bad dist");
  }

  if (doPrintBounds) {
	cout << "x0 = " << mystable::vec2str(n,x) << endl;
	cout << "lb = " << mystable::vec2str(n,lb) << endl;
	cout << "ub = " << mystable::vec2str(n,ub) << endl;
  }

}


//-------------------------------------------------------------------------------
// MLE of distribution paramters

mystable::ts_struct mystable::mle_nlopt(vec y, mystable::dist dist) {

  // sort input data
  //y = sort(y);
  
  int iPars = mystable::getParCount(dist);
  int iObs = y.n_rows;

  vector<double>  x(iPars); // parameter vector
  vector<double> lb(iPars); // lower param bounds
  vector<double> ub(iPars); // upper param bounds

  getInitialPars(iPars, x, lb, ub, dist);

  fit_struct fs;
  fs.ny = iObs;
  fs.np = iPars;
  fs.y = y.memptr();
  fs.dist = dist;
  fs.fxcalls = 0;
  fs.iters = 0;
  
  nlopt::algorithm algo = nlopt::LN_NELDERMEAD;
  //nlopt::algorithm algo = nlopt::LN_COBYLA;
  nlopt::opt opt(algo, iPars);
  opt.set_xtol_rel(1e-6);
  opt.set_ftol_abs(1e-6);
  opt.set_maxeval(5000);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);
  opt.set_min_objective(negLLF_stable, &fs);
 
  double minf;
  double time = MPI::Wtime();
  nlopt::result nloptResult;
		
  nloptResult = opt.optimize(x, minf);
  double runTime = MPI::Wtime() - time;
  
  ts_struct s;
  s.pars = new double[iPars];
  for (int i = 0; i < iPars; i++)
	s.pars[i] = x[i];
  s.iters = fs.iters;
  s.fxcalls = fs.fxcalls;
  s.LLF = minf;
  s.status = (int) nloptResult;
  s.runTime = runTime;
  s.dist = dist;

  //<TODO> Free fit_struct and ts_struct later
		 	   
  return s;
}


//-------------------------------------------------------------------------------

void mystable::ts_struct_print(mystable::ts_struct s) {
  cout << "-----------------------------------------------------------\n";
  cout << "xF      = " << mystable::arr2str(mystable::getParCount(s.dist), s.pars) << endl;
  cout << "iters   = " << s.iters << endl;
  cout << "fxcalls = " << s.fxcalls << endl;
  cout << "negLLF  = " << s.LLF << endl;
  cout << "status  = " << s.status << endl;
  cout << "runTime = " << s.runTime << endl;
  cout << "-----------------------------------------------------------\n";
}


//-------------------------------------------------------------------------------
void ts_struct_free(mystable::ts_struct s) {
  
}

//-------------------------------------------------------------------------------
// Alpha-stable sub-Gaussian (ASSG) global alpha parameter

double mystable::assg_alphaEst(vec alpha) {
  return mean(alpha);
}


//-------------------------------------------------------------------------------
// ASSG estimate of sigma for linear projetion a'*X
// This is from Kring Dissertation, Theorem 11.

double mystable::assg_scaleEst(mat X, double alpha, vec mu, vec a) {

  int N = X.n_rows;
  double n = (double) N;
  double c = (2./M_PI)*sin(M_PI/(2.*n)) * xGamma(1./n) * xGamma(1. - 1./(n*alpha));
  double aMu = sum(a % mu); //dot(a, mu);
  vec inner(N);
  
  for (int i = 0; i < N; i++) {
	double aX = sum(a % X.row(i).t()); //dot(a, X.row(i).t());
	inner(i) = pow(abs(aX - aMu), 1./n);
  }
  return pow(c,-1.0*n)*prod(inner);
}


//-------------------------------------------------------------------------------
// ASSG dispersion matrix estimator

void mystable::assg_dispersionEst(mat X, double alpha, vec sigma, vec mu, mat &Sigma) {

  assert(sigma.n_rows == mu.n_rows);
  const int n = mu.n_rows;
  
  //mat Sigma(n,n);

  // set main and upper diagonal
  for (int i = 0; i < n; i++) {
	for (int j = i; j < n; j++) {
	  if (i == j) {
		Sigma(i,i) = 2.0*sigma(i)*sigma(i);
	  } else {
		vec a = zeros(n);
		a(i) = 1.0;
		a(j) = 1.0;
		double sigmaPlus = assg_scaleEst(X, alpha, mu, a);
		a.zeros();
		a(i) = 1.0;
		a(j) = -1.0;
		double sigmaMinus = assg_scaleEst(X, alpha, mu, a);
		Sigma(i,j) = (sigmaPlus*sigmaPlus - sigmaMinus*sigmaMinus)/2.;		
	  }
	}
  }

  // copy upper diagonal to lower
  for (int i = 1; i < n; i++)
	for (int j = 0; j < i; j++)
	  Sigma(i,j) = Sigma(j,i);
}


//-------------------------------------------------------------------------------
//  RNG of stable and TS variables.

mat mystable::stablernd(const int nparam, const double param[], mystable::dist dist,
						int m, int n) {

  /*
  // random uniforms
  vec u(iRows*iCols);
  u.randu();

  u.save("/Users/jimmiegoode/Documents/Glimm/Toolbox/alphastable/testing/e_u.csv", csv_ascii);

  vec inv = mystable::inv_FFT(u, nparam, param, dist);

  inv.save("/Users/jimmiegoode/Documents/Glimm/Toolbox/alphastable/testing/e_inv.csv", csv_ascii);
 
  mat out = reshape(inv, iRows, iCols);

  return out; 
  */

  mat x(m, n);

  //double pi = M_PI;
  
  switch (dist) {
  case mystable::AS:
	{
	  double alpha = param[0];
	  double beta = param[1];
	  double c = param[2];
	  double delta = param[3];
	  
	  mat w = -1.0 * arma::log(Stats::rand(m,n));
	  mat phi = (Stats::rand(m,n)-.5)*M_PI;

	  if (alpha == 2.) {
		x = (2.*arma::sqrt(w) % arma::sin(phi));
		x = delta + c*x;
		return x;
	  }

	  // Symmetric cases:
	  if (beta == 0.) {
		if (alpha == 1.) { // Cauchy case 
		  x = arma::tan(phi);
		} else {
		  x = ( pow((cos((1.-alpha)*phi) / w), (1./alpha - 1.)) 
				% sin(alpha * phi) / pow(cos(phi), (1./alpha)));
		}

	  // General cases:
	  } else {

		mat cosphi = cos(phi);

		if (abs(alpha-1) > 1.e-8) {

		  double zeta = beta * tan(M_PI*alpha/2.);
		  mat aphi = alpha * phi;
		  mat a1phi = (1. - alpha) * phi;

		  x = ((sin(aphi) + zeta * cos(aphi)) / cosphi) 
            % arma::pow(((cos(a1phi) + zeta * sin(a1phi)) / (w % cosphi)), ((1.-alpha)/alpha));
		  
		} else {
		  
		  mat bphi = (M_PI/2.) + beta * phi;

		  x = (2./M_PI) * (bphi % tan(phi) - beta * log((M_PI/2.) * w % cosphi / bphi));

		  if (alpha != 1.) {
            x = x + beta * tan(M_PI * alpha/2.);
		  }
		}			
	  }

	  // Finale:
	  x = delta + c * x;
	  return x;
	  
	}
	break;
  default:
	assert(0 && "RNG for dist not implemented");
  }
}

//-------------------------------------------------------------------------------



//-------------------------------------------------------------------------------
// ASSG sample

mat mystable::assg_rnd(double alpha, mat Sigma, vec mu, int nSim) {

  // check dimensions
  assert(Sigma.n_rows == Sigma.n_cols);
  assert(Sigma.n_rows == mu.n_rows);

  int iCols = mu.n_rows;

  double param[] = {alpha/2., 1., pow(cos(M_PI*alpha/4.), 2./alpha), 0.};

  mat W = mystable::stablernd(4, param, mystable::AS, nSim, iCols);

  return repmat(mu.t(), nSim, 1) +
	arma::sqrt(W) % Stats::mvnrnd(arma::zeros(iCols), Sigma, nSim);
}

//-------------------------------------------------------------------------------
// ASSG sample - Cholesky factor A supplied

mat mystable::assg_rnd_2(double alpha, mat A, vec mu, int nSim) {

  // check dimensions
  assert(A.n_rows == A.n_cols);
  assert(A.n_rows == mu.n_rows);

  int iCols = mu.n_rows;

  double param[] = {alpha/2., 1., pow(cos(M_PI*alpha/4.), 2./alpha), 0.};

  mat W = mystable::stablernd(4, param, mystable::AS, nSim, iCols);

  return repmat(mu.t(), nSim, 1) +
	arma::sqrt(W) % Stats::mvnrnd_2(arma::zeros(iCols), nSim, A);
}



//-------------------------------------------------------------------------------
// NOTE: This updates resid in place with GARCH filtered returns.

void mystable::portSample(const int nSim,
						  const double alpha,
						  const vec mu,
						  const mat mnGarchPars,
						  const vec wts,
						  const mat A, // A is chol(Sigma)
						  mat &resid,
						  vec &nextMeans,
						  vec &nextSigmas,
						  vec &portRets) {

  // number of stocks
  int iS = mu.n_rows;

  resid = mystable::assg_rnd_2(alpha, A, mu, nSim);

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

void mystable::tester(int c) {

  cout << ">> Testing mystable.cpp (Case " << c << ")" << endl;
  
  string fname = "/Users/jimmiegoode/Documents/Glimm/data/scalingSeries/1.csv";
  //string fname = "/Users/jimmiegoode/Documents/Papers/paper2/Toolbox_Project/rand_AS.csv";
  //string fname = "/Users/jimmiegoode/Documents/Papers/paper2/Toolbox_Project/rand_AS2.csv";
  //string fname = "/Users/jimmiegoode/Documents/Papers/paper2/Toolbox_Project/rand_AS3.csv";
  //string fname = "/Users/jimmiegoode/Documents/Papers/paper2/Toolbox_Project/rand_symAS.csv";
  
  vec logret;

  //cout << "mom =    [" << mean(logret) << ", " << var(logret) << "]" << endl;
  //cout << "minmax = [" << min(logret)  << ", " << max(logret) << "]" << endl;
  
  double x[3];
  
  switch (c) {
  case 1:
	{
	  vec u = linspace(-2,2,10);
	  double alpha = 1.5;
	  double lmp = 0.8;
	  double lmm = 1.2;
  
	  cx_vec phi = chf_stdCTS(u, alpha, lmp, lmm);

	  phi.print("phi_CTS = ");
	}
	break;
	
  case 2:
	{

	  //
	  // Test MLE
	  //

	  logret.load(fname, csv_ascii);
	  
	  // [1.5 0.05 0.05], [STOL, STOL, STOL], [2-STOL 10 10]

	  vec x0(3), xL(3), xU(3);
	  double TOL = 0.01;
	
	  x0(0) = 1.5;
	  xL(0) = TOL;
	  xU(0) = 2-TOL;

	  x0(1) = 0.05;
	  xL(1) = TOL;
	  xU(1) = 10.0;

	  x0(2) = 0.05;
	  xL(2) = TOL;
	  xU(2) = 10.0;

	  mystable::ts_struct s = mystable::mle_nlopt(logret, mystable::stdCTS);
	  mystable::ts_struct_print(s);
	  
	}
	break;
	
  case 3:
	{
	  logret.load(fname, csv_ascii);

	  x[0] = 0.77987;
	  x[1] = 0.010001;
	  x[2] = 0.010193;

	  logret = sort(logret);

	  //vec pdf = pdf_FFT(logret.n_rows, 3, logret.memptr(), x, mystable::stdCTS);
	  //cout << "negLLF = " << -1.0*sum(log(pdf)) << endl;
	  
	  // double xx[1];
	  // xx[0] = 0.014927;
	  // cout << "f(0.014927) = " << pdf_FFT(1, 3, xx, x, mystable::CTS);
	}
	break;

  case 4:
	{
	  vec u = linspace(-2,2,5);

	  cout << chf_stdNTS(u, 1.8, 1.0, -0.2) << endl;
	  cout << chf_stdCTS(u, 1.5, 0.05, 0.1) << endl;

	  cout << "--------------------" << endl;

	  cx_double dt(1.0, 0.0);
	  cx_double onei(0.0, 1.0);

	  cout << dt << endl;
	  cout << onei << endl;
	}
	break;

  case 5:
	{
	  logret.load(fname, csv_ascii);
	  
	  mystable::ts_struct s = mystable::mle_nlopt(logret, mystable::stdNTS);
	  mystable::ts_struct_print(s);
	  
	}
	break;

  case 6:
	{
	  //chf_AS(vec u, double alpha, double beta, double sigma, double mu)
	  vec u = linspace(-2,2,5);

	  cx_vec phi = chf_AS(u, 1.5, -0.5, 1.3, 0.1);

	  phi.print("phi = ");
	}
	break;

  case 7:
	{
	  logret.load(fname, csv_ascii);
	  
	  doPrintBounds = 1;
	  doPrintLLF = 1;
	  mystable::ts_struct s = mystable::mle_nlopt(logret, mystable::stdAS);
	  mystable::ts_struct_print(s);
	}
	break;

  case 8:
	{
	  logret.load(fname, csv_ascii);
	  
	  doExport = 1;
	  
	  double xa[2];
	  xa[0] = 1.5;
	  xa[1] = -0.5;
	  
	  logret = sort(logret);

	  //vec pdf = pdf_FFT(logret.n_rows, 2, logret.memptr(), xa, mystable::stdAS);
	  //cout << "negLLF = " << -1.0*sum(log(pdf)) << endl;
	  
	}
	break;

  case 9:
	{
	  logret.load(fname, csv_ascii);
	  
	  double t0 = MPI::Wtime();
	  
	  doPrintBounds = 1;
	  doPrintLLF = 0;
	  mystable::ts_struct s = mystable::mle_nlopt(logret, mystable::symAS);

	  cout << "t_fit = " << MPI::Wtime()-t0 << endl;
	  
	  mystable::ts_struct_print(s);
	}
	break;
  case 10:
	{
	  mat X;
	  X.load("/Users/jimmiegoode/Documents/PhD/Matlab/ASSG/rand_ASSG.csv", csv_ascii);

	  int n = X.n_cols;

	  vec alpha(n);
	  vec sigma(n);
	  vec mu(n);

	  for (int i = 0; i < n; i++) {
		mystable::ts_struct s = mystable::mle_nlopt(X.col(i), mystable::symAS);
		mystable::ts_struct_print(s);

		alpha(i) = s.pars[0];
		sigma(i) = s.pars[1];
		mu(i)    = s.pars[2];
	  }
	  
	  double nalpha = mean(alpha);
	  mat SigmaHat(n,n);
	  mystable::assg_dispersionEst(X, nalpha, sigma, mu, SigmaHat);

	  cout << "alpha = " << nalpha << endl;
	  mu.print("mu = ");
	  SigmaHat.print("Sigma = ");
	  
	}
	break;

  case 11:
	{
	  /*
	  vec X;
	  X.load("/Users/jimmiegoode/Documents/Glimm/data/scalingSeries/7.csv", csv_ascii);

	  // double t0 = MPI::Wtime();
	  // vec Fx = Stats::empCDF(X);
	  // cout << MPI::Wtime() - t0 << " sec \n";

	  double t0 = MPI::Wtime();
	  vec Fx = Stats::empCDF_fast(X);
	  cout << MPI::Wtime() - t0 << " sec \n";
	  
	  Fx.save("/Users/jimmiegoode/Documents/Glimm/skewt_vs_ASSG/Fx.csv", csv_ascii);
	  
	  //Fx.print("Fx = ");
	  */
	}
	break;
	
  case 12:
	{
	  doPrintBounds = 0;

	  vec X;
	  X.load("/Users/jimmiegoode/Documents/Glimm/Toolbox/alphastable/testing/rnd_stdAS.csv", csv_ascii);
	  string fout = "/Users/jimmiegoode/Documents/Glimm/Toolbox/alphastable/testing/";

	  vec u = linspace(-10., 10., 1000);
	  vec p = linspace(0.01, 0.99, 1000);

	  int iN = mystable::gridSize;
	  mystable::dist dist;
      int npars = 0;
	  mystable::ts_struct s;
	  
	  vector<mystable::dist> dists;
	  dists.push_back(mystable::stdAS);
	  dists.push_back(mystable::stdCTS);
	  dists.push_back(mystable::stdNTS);
	  dists.push_back(mystable::AS);
	  dists.push_back(mystable::symAS);

	  for (int i = 0; i < dists.size(); i++)
	  {	
		dist = dists[i];
		string sDist = mystable::dist2str(dist);
		npars = mystable::getParCount(dist);
	  
		s = mystable::mle_nlopt(X, dist);
		mystable::ts_struct_print(s);

		// CDF
		vec cdf = mystable::cdf_FFT(u, npars, s.pars, dist);
		cdf.save(fout + "e_cdf_" + sDist + ".csv", csv_ascii);

		// INV
        vec q = mystable::inv_FFT(p, npars, s.pars, dist);
		q.save(fout + "e_inv_" + sDist + ".csv", csv_ascii);

		// PDF
		vec f(u.n_rows);
		vec x(iN);
		vec pdf(iN);
		pdf_FFT(u.n_rows, npars, u.memptr(), s.pars, dist, f, x, pdf);
		f.save(fout + "e_f_" + sDist + ".csv", csv_ascii);	
	  }
	}
	break;
	
  case 13:
	{
	  vec arg;
	  arg.load("e_res.csv", csv_ascii);

	  cout << "mean(arg) = " << mean(arg) << ", " << max(arg) << ", " << min(arg) << endl;
	  
	  vec param;
	  param << 1.0e-6 << 1.0e-6 << 1.0e-6 << endr; 
	  
	  vec icdf = mystable::inv_FFT(arg, 3, param.memptr(), mystable::stdCTS);

	  icdf.save("e_icdf.csv", csv_ascii);
	}
	break;

  case 14:
	{

	  mat X;
	  int n;
	  
	  if (0) {

		cout << "----> Importing\n";
		
		X.load("/Users/jimmiegoode/Documents/PhD/Matlab/ASSG/rand_ASSG.csv", csv_ascii);
		n = X.n_cols;
		
	  } else {

		cout << "----> Simulating\n";
	  
		n = 2;

		double nAlpha = 1.5;

		mat mnSigma(n,n);  
		mnSigma << 5 << 11 << endr
				<< 11 << 25 << endr;

		vec vnMu(n);
		vnMu.zeros();

		X = mystable::assg_rnd(nAlpha, mnSigma, vnMu, 10000);
	  }


	  mat sampleMean = mean(X);
	  sampleMean.print("sampleMean = ");
	  
	  vec alpha(n);
	  vec sigma(n);
	  vec mu(n);

	  for (int i = 0; i < n; i++) {
	  	mystable::ts_struct s = mystable::mle_nlopt(X.col(i), mystable::symAS);
	  	mystable::ts_struct_print(s);

	  	alpha(i) = s.pars[0];
	  	sigma(i) = s.pars[1];
	  	mu(i)    = s.pars[2];
	  }
	  
	  double nalpha = mean(alpha);
	  mat SigmaHat(n,n);
	  mystable::assg_dispersionEst(X, nalpha, sigma, mu, SigmaHat);

	  cout << "alpha = " << nalpha << endl;
	  mu.print("mu = ");
	  SigmaHat.print("Sigma = ");
	}
	break;

  case 15:
	{
	  double alpha = 1.5;
	  
	  //double param[] = {alpha/2., 1., pow(cos(M_PI*alpha/4.), 2./alpha), 0.};
	  //double param[] = {1.5, 1., 0.75, 0.};
	  double param[] = {alpha/2., 1., pow(cos(M_PI*alpha/4.), 2./alpha), 0.};
 
	  mat W = mystable::stablernd(4, param, mystable::AS, 5000, 1);
	  
	  //W.save("/Users/jimmiegoode/Documents/Glimm/Toolbox/alphastable/testing/e_W.csv", csv_ascii);

	  mystable::ts_struct s = mystable::mle_nlopt(W.col(0), mystable::AS);
	  mystable::ts_struct_print(s);
	  
	}
	break;
	
  } // end of swtich
}

//-------------------------------------------------------------------------------


// int main(int argc, char** argv) {
  
//   if (argc < 2)
// 	assert(0 && "Invalid Arguments");
  
//   int testCase = atoi(argv[1]);

//   tester(testCase);
  
//   return 0;
// }
