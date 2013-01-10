/*
  Author: Jimmie Goode
  Created: 2012-08-21
*/

#include "mygarch.hpp"

//-------------------------------------------------------------------------------
//

bool debug = false;

double MPI_Wtime();

//-------------------------------------------------------------------------------
// Print a progress bar

void printProgressBar( int percent ){
  std::string bar;

  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }

  std::cout<< "\r" "[" << bar << "] ";
  std::cout.width( 3 );
  std::cout<< percent << "%     " << std::flush;
}

//-------------------------------------------------------------------------------
//
  
string arr2str(const unsigned n, const double* x) {
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

string vec2str(const unsigned n, const vector<double> x) {
  double xarr[n];
  std::copy(x.begin(),x.end(),xarr);
  return arr2str(n, xarr);
}

//-------------------------------------------------------------------------------
//

void printStruct(garch_struct s) {
  cout << "n = " << s.n;
  cout << "iters = " << s.iters;
  cout << "idx = " << s.idx;
}


//-------------------------------------------------------------------------------
// Log-likelihood function

double mygarch::negLLF_struct(const double* x, garch_struct *s) {

  // Allocate named parameter values
  double C;             // Conditional mean constant
  gsl_vector *AR;       // Conditional mean AR coefficients
  gsl_vector *MA;       // Conditional mean MA coefficients
  double K;             // Conditional variance constant  
  gsl_vector *GARCH;    // Conditional variance coefficients for lagged variances
  gsl_vector *ARCH;     // Conditional variance coefficients for lagged squared innovations
  double DoF;           // Degrees of freedom for Student's t innovations.

  if (s->R > 0)  AR    = gsl_vector_alloc(s->R);
  if (s->M > 0)  MA    = gsl_vector_alloc(s->M);
  if (s->P > 0)  GARCH = gsl_vector_alloc(s->P);
  if (s->Q > 0)  ARCH  = gsl_vector_alloc(s->Q);

  // Set parameters from x[] to named values.
  int R = s->R;
  int M = s->M;
  int P = s->P;
  int Q = s->Q;
  
  C = x[0];

  for (int i = 1; i < 1+R; i++) 
	gsl_vector_set(AR, i-1, x[i]);

  for (int i = 1+R; i < 1+R+M; i++)
	gsl_vector_set(MA, i-(1+R), x[i]);

  K = x[1+R+M];

  int i0 = 1+R+M+1;
  
  for (int i = i0; i < i0+P; i++)
	gsl_vector_set(GARCH, i-i0, x[i]);

  i0 = 1+R+M+1+P;
  
  for (int i = i0; i < i0+Q; i++)
	gsl_vector_set(ARCH, i-i0, x[i]);

  DoF = x[(s->n)-1];
  
  // allocate padded vectors for returns (yy), innovations (ee), and variances (hh)
  int iPad = s->yy->size;
  
  // This is the unconditional variance used to initialize GARCH recursion.
  double sigma0 = sqrt(K/(1.0 - Stats::sum(ARCH) - Stats::sum(GARCH)));

  // This is the unconditional innovation used to initialize ARMA-GARCH recursion.
  double e0 = sigma0;

  // initialize first maxRMPQ elements of e and h
  for (int i = 0; i < s->maxRMPQ; i++) {
	gsl_vector_set(s->ee, i, e0);
	gsl_vector_set(s->hh, i, sigma0*sigma0);
  }

  double negLL = 0.0;
  double myPI = 3.14159265359;
  double cons = gsl_sf_gamma((DoF + 1.0)/2.0) / (sqrt(DoF*myPI) * gsl_sf_gamma(DoF/2.0));
  
  // filter out innovations and conditional variances
  for (int t = s->maxRMPQ; t < iPad; t++) {

	double sumAR = 0.0;
	for (int j = 0; j < s->R; j++) {
	  sumAR += gsl_vector_get(AR, j) * gsl_vector_get(s->yy, t-j-1);
	}
	double sumMA = 0.0;
	for (int j = 0; j < s->M; j++) {
	  sumMA += gsl_vector_get(MA, j) * gsl_vector_get(s->ee, t-j-1);
	}

	// e(t) = y(t) - C - sumAR - sumMA;
	double e_t = gsl_vector_get(s->yy, t) - C - sumAR - sumMA;
	gsl_vector_set(s->ee, t, e_t);

	double h_t = K;
	for (int j = 0; j < s->P; j++) {
	  h_t += gsl_vector_get(GARCH, j) * gsl_vector_get(s->hh, t-j-1);
	}
	for (int j = 0; j < s->Q; j++) {
	  double et = gsl_vector_get(s->ee, t-j-1);
	  h_t += gsl_vector_get(ARCH, j) * et*et;
	}
	gsl_vector_set(s->hh, t, h_t);

	h_t = sqrt(abs(h_t));
	double u_t = e_t / h_t;

	// student's t
	negLL -= log( cons*pow(1 + (u_t*u_t)/DoF, -1.0*(DoF + 1.0)/2.0)/h_t );
	
  } // for t

  if (s->R > 0)  gsl_vector_free(AR);
  if (s->M > 0)  gsl_vector_free(MA);
  if (s->P > 0)  gsl_vector_free(GARCH);
  if (s->Q > 0)  gsl_vector_free(ARCH);

  return negLL;
}



//-------------------------------------------------------------------------------
//

double negLLF_par(double x, void *params) {
  
  // cast void pointer to parameter structure
  garch_struct* s = (garch_struct *) params;
  
  // copy existing x[] array into xx[] so that it can be changed.
  double xx[s->n];
  for (int i = 0; i < s->n; i++)
	xx[i] = s->x[i];
  
  // set the idx element to double x.
  xx[s->idx] = x;

  double dNegLL = mygarch::negLLF_struct(xx, s);
  
  return dNegLL;
}


//-------------------------------------------------------------------------------
// compute gradient of garch likelihood function

gsl_vector* mygarch::negLLF_grad(const unsigned n, const double* x,  garch_struct *s, const double *steps) {
  
  //double* grad = new double[n];
  gsl_vector* grad = gsl_vector_alloc(n);
  
  gsl_function F;
  double dx, dxerr;
  
  // <TODO> Inefficient to create a new structure.
  garch_struct ss;
  ss.n = s->n;
  ss.R = s->R;
  ss.M = s->M;
  ss.P = s->P;
  ss.Q = s->Q;
  ss.maxRMPQ = s->maxRMPQ;
  ss.iters = s->iters;
  ss.yy = s->yy;
  ss.ee = s->ee;
  ss.hh = s->hh;
  
  double xx[n];
  for (int i = 0; i < n; i++)
	xx[i] = x[i];
  ss.x = xx;
  
  for (int i = 0; i < n; i++) {
	ss.idx = i;
	F.function = &negLLF_par;
	F.params = &ss;
	gsl_deriv_central(&F, x[i], steps[i], &dx, &dxerr);
	
	//grad[i] = dx;
	gsl_vector_set(grad, i, dx);
  }

  mygarch::garch_struct_free(ss);
  
  return grad;
}


//-------------------------------------------------------------------------------
//

double negLLF_nlopt(unsigned n, const double *x, double *grad, void *params) {

  garch_struct* s = (garch_struct *) params;
  s->fxCalls++;

  bool isBadGrad = false;
  
  if (grad) {

	assert(0 && "Gradient not yet implemented");
	/*
	// Step sizes for differentiation
	double steps[7] = {1.0e-8, 1.0e-6, 1.0e-6, 1.0e-9, 1.0e-8, 1.0e-6, 1.0e-6};
	
	gsl_vector* g = mygarch::negLLF_grad(n, x, s, steps);
	
	for (int i = 0; i < n; i++) {
	  grad[i] = gsl_vector_get(g,i);
	  if (gsl_finite(grad[i]) != 1) {
		isBadGrad = true;
	  }
	}
	*/
  }
    
  //double dNegLL = mygarch::negLLF(n, x, s->y, s->R, s->M, s->P, s->Q);
  double dNegLL = mygarch::negLLF_struct(x, s);

  if (gsl_finite(dNegLL) != 1) { 
	//cout << "----> BAD_LL = " << dNegLL << ", setting to 1e20. pers = " << x[4] + x[5] << endl;
	//if (x[4] + x[5] >= 1.0)
	  return 1.0e20;
  }

  if (grad) {
	cout << "GRAD = " << arr2str(n, grad) <<  ", NEGLL = " << dNegLL << endl;
	if (isBadGrad) {
	  cout << "BAD GRAD = " << arr2str(n, grad) << endl;
	  // cout << "\tX    = " << arr2str(n, x) << endl;
	  cout << "\tpers = " << x[4] + x[5] << endl;
	  cout << "\tnegLL = " << dNegLL << endl;
	}
  }

  return dNegLL;
}



//-------------------------------------------------------------------------------
// Extract residuals from ARMA-GARCH fit

gsl_vector* garch_struct_res(garch_struct s) {

  int iPad = s.yy->size;
  
  gsl_vector* u = gsl_vector_alloc(iPad - s.maxRMPQ);

  for (int t = s.maxRMPQ; t < iPad ; t++) {

	// abs(h_t) just in case
	double u_t = gsl_vector_get(s.ee, t) / sqrt(abs(gsl_vector_get(s.hh, t)));
	
	gsl_vector_set(u, t-s.maxRMPQ, u_t);
  }

  return u;
}

//-------------------------------------------------------------------------------
//

void mygarch::garch_struct_free(garch_struct s) {
  gsl_vector_free(s.yy);
  gsl_vector_free(s.ee);
  gsl_vector_free(s.hh);
  gsl_vector_free(s.u);
  gsl_vector_free(s.x0);
  gsl_vector_free(s.x1);
}


//-------------------------------------------------------------------------------
//

void mygarch::garch_struct_print(garch_struct s) {
  cout << "-----------------------------------------------------------\n";
  cout << "t-ARMA(" << s.R << "," << s.M << ")-GARCH(" << s.P << "," << s.Q << ")" << endl;
  cout << "\talgo = " << s.algo << endl;
  cout << "\trunTime = " << s.runTime << endl;
  cout << "\tExited with status = " << s.exitStatus << endl;
  cout << "\tnegLL = " << s.negLL << endl;
  cout << "\tx0 = " << arr2str(s.n, s.x0->data) << endl;
  cout << "\tx1 = " << arr2str(s.n, s.x1->data) << endl;
  cout << "\titers = " << s.iters << endl;
  cout << "\tfxCalls = " << s.fxCalls << endl;
  cout << "\tiObs = " << s.yy->size << endl;
  cout << "-----------------------------------------------------------\n";
}



//-------------------------------------------------------------------------------
//
/*
void test_scaling() {
  
  string dir = "/Users/jimmiegoode/Dropbox/Projects/Glimm/data/scalingSeries/";

  for (int i = 1; i <= 7; i++) {

	gsl_matrix* mnRet = IO::importCSVmatrix(dir + boost::lexical_cast<string>(i)+".csv");
	gsl_vector* vnRet = gsl_vector_alloc(mnRet->size1);
	gsl_matrix_get_col(vnRet, mnRet, 0);

	double t;
	mygarch::fit_nlopt(vnRet, 1, 1, 1, 1);
	t = clock();
	float dt = ((float)t)/CLOCKS_PER_SEC;
	
	cout << "DT = " << dt << endl;
	
	gsl_vector_free(vnRet);
	gsl_matrix_free(mnRet);
  }
  cout << endl;
}
*/

//-------------------------------------------------------------------------------
//

void garch0(const int P, const int Q, const double unconVar, vec &garch, vec &arch, double &K) {

  for (int i = 0; i < P; i++)
	garch(i) = 0.85/max(((double) P), 1.0);

  double diff = 0.90 - sum(garch);
  
  for (int i = 0; i < Q; i++)
	arch(i) = diff/max(((double) Q), 1.0);

  if (unconVar <= 0)
	K = 1e-3; // A reasonable assumption for daily returns
  else
	K = unconVar*(1 - sum(garch) - sum(arch));
  
}


//-------------------------------------------------------------------------------
//   Compute initial estimates of the autoregressive (AR) and moving
//   average (MA) coefficients of a stationary/invertible univariate ARMA
//   model. Estimates of the model constant and the variance of the
//   innovations process are also provided. This function provides initial
//   coefficient estimates suitable for further refinement via maximum 
//   likelihood estimation.

void arma0(const gsl_vector *y, const int R, const int M, double &constant,
		   double &variance, vec &AR, vec &MA) {

  using namespace arma;
  
  int n = y->size;
  
  vec vnY(n);

  // <TODO> move this copy of y into main fit
  for (int i = 0; i < n; i++) {
	vnY(i) = gsl_vector_get(y, i);
  }

  // Check for the case of no ARMA model. In this case, compute the sample
  // mean and variance and exit:

  if (R + M == 0) {
	constant = mean(vnY);
	variance = var(vnY);
	return;
  }
  
  // Estimate the AR coefficients of a general ARMA(R,M) model:
  
  if (R > 0) {
	  
	vec correlation(R + M + 1);
	for (int i = 0; i <= R+M; i++) {
	  vec c  = cor(vnY.rows(i, n-1), vnY.rows(0, n-1-i));
	  correlation(i) = c(0);
	}

	// Compute the autocovariance sequence of the y(t) process. The variance
	// (zeroth-lag autocovariance) is found in the first element:
	  
	variance = var(vnY);
	vec covariance = correlation * variance;

	// i = abs(M:-1:M-R+1) + 1; % covariance(k) = covariance(-k)
    vector<double> idxs;	
	for(int i = M; i >= M-R+1; i--) {
	  idxs.push_back(abs(i) + 1);
	}
	vec vnI(idxs.size());
	for (int i = 0; i < vnI.n_rows; i++)
	  vnI(i) = idxs[i];
	
	vec c1(vnI.n_rows);
	vec c2(vnI.n_rows);
	vec c3(vnI.n_rows);

	if (M > 0) {
	  
	  // For ARMA processes, the matrix C of covariances derived from the 
      // estimated autocovariancRe sequence is Toeplitz, but nonsymmetric.
      // The AR coefficients are then found by solving the modified Yule-
      // Walker equations.

	  for (int i = 0; i < vnI.n_rows; i++) {
		c1(i) = covariance(M + 1 + i - 1);
		c2(i) = covariance(vnI(i) - 1);
		c3(i) = covariance(M + 2 + i -1);
	  }
	  
	  mat C = toeplitz(c1, c2);
  
	  if (R == 1) {
		AR = c3/C;
	  } else {
		AR = solve(C, c3);
	  }
	  
	} else {

	  if (R == 1) {		

		AR(0) = correlation(1);
		
	  } else {
		
		// For AR processes, the matrix C of covariances derived from the 
        // estimated auto-covariance sequence is Toeplitz and symmetric.
        // The AR coefficients are found by solving the Yule-Walker
        // equations.

		cout << "AR only\n" ;
		
		mat C = toeplitz(covariance.rows(0, R-1));
		AR = solve(C, covariance.rows(1, R));
		
	  }
	}
	
	// Ensure the AR process is stationary. If not, then set all ARMA
	// coefficients to 0. This ensures the subsequent optimization will begin
	// with a stationary/invertible ARMA model for the conditional mean.

    const int m = AR.n_rows + 1;
	double a[m];
	for (int i = 0; i < m-1; i++) {
	  a[i] = -AR(AR.n_rows-1-i);
	}
	a[m-1] = 1.0;

	double z[2*(m-1)];     
	gsl_poly_complex_workspace *w  = gsl_poly_complex_workspace_alloc(m);
	gsl_poly_complex_solve(a, m, w, z); 
	gsl_poly_complex_workspace_free(w);

	for (int i = 0; i < m-1; i++) {
	  double r = z[2*i];   // real part
	  double c = z[2*i+1]; // complex part
	  gsl_complex zi = gsl_complex_rect(r,c);
	  double ziabs = gsl_complex_abs(zi);

	  if (ziabs >= 1.0) {
		for (int i = 0; i < AR.n_rows; i++)
		  AR(i) = 0.0;
		for (int i = 0; i < MA.n_rows; i++)
		  MA(i) = 0.0;
		
		variance = var(vnY);
		constant = mean(vnY);
		return;
	  }
	  //printf ("z%d = %+.18f %+.18f, abs = %g, abs2 = %g\n", i, r, c, ziabs, sqrt(r*r + c*c));	   
	}	
  }

  // Filter the ARMA(R,M) input series y(t) with the estimated AR coefficients 
  // to obtain a pure MA process. If the input moving-average model order M is
  // zero (M = 0), then the filtered output is just a pure innovations process
  // (i.e., an MA(0) process); in this case the innovations variance estimate
  // is just the sample variance of the filtered output. If M > 0, then
  // compute the autocovariance sequence of the MA process and continue.

  vec temp(R+1);
  temp(0) = 1.0;
  for (int i = 1; i <= R; i++)
	temp(i) = -AR(i-1);

  // Convolution is same as matlab's filter() function (FIR),
  // except the last M terms dropped. This is a Finite Input Response (FIR).
  vec xx = conv(temp, vnY); 
  vec x = xx.rows(0, xx.n_rows-1-R);
  
  // x.print("x = ");
  // cout << x.n_rows << ", " << x.n_cols << endl;
  // temp.print("temp = ");

  constant = mean(x);

  if (M == 0) {
	variance = var(x);
	return;
  }

  // Covariance of an MA(M) process

  vec correlation(M+1);
  for (int i = 0; i <= M; i++) {
	vec c  = cor(x.rows(i, n-1), x.rows(0, n-1-i));
	correlation(i) = c(0);
  }

  vec c = correlation * var(x);

  // Estimate the variance of the white noise innovations process e(t)
  // and the MA coefficients of a general ARMA(R,M) model. The method of
  // computation is outlined in equation A6.2.4 (page 221) of [1].

  for (int i = 0; i < M; i++)
	MA(i) = 0.0;
  
  vec MA1(M);        // Saved MA coefficients from previous iteration
  int counter = 1;   // Iteration counter
  double tol = 0.05; // Convergence tolerance
  
  while ((norm(MA - MA1, 2) > tol) && (counter < 100)) {

	MA1 = MA;

	// Estimate the variance of the innovations process e(t):

	double den = 1.0;
	for (int i = 0; i < M; i++)
	  den += MA(i)*MA(i);
	
	variance = c(0)/den;

	if (abs(variance) < tol)
	  break;

	for (int j = M; j >= 1; j--) {

	  //Estimate the moving-average coefficients. The MA coefficients are
	  // the negative of those appearing in equation A6.2.4 (page 221) of
	  // [1]. This is due to the convention of entering coefficient values,
	  // via GARCHSET, exactly as the equation is written.
	  
	  double sum = c(j) * (1/variance);

	  if (M - j != 0) {
		
		vec t1(M-j);
		vec t2(M-j);

		for (int i = 0; i < M-j; i++)
		  t1(i) = -MA(i);

		int idx = 0;
		for (int i = j; i < M; i++) {
		  t2(idx) = MA(i);
		  idx++;
		}
	  
		sum += dot(t1, t2);
	  }

	  MA(j-1) = sum; 
	}
	counter++;
  }
}


//-------------------------------------------------------------------------------
//

void mygarch::getInitialPars(const int n, const gsl_vector* y,
							 const int R, const int M, const int P, const int Q,
							 double x[], double lb[], double ub[]) {


  double GARCH_TOLERANCE = 2e-7;
  
  vec AR(R);
  vec MA(M);
  double constant, variance, K;

  arma0(y, R, M, constant, variance, AR, MA);
  
  vec GARCH(P);
  vec ARCH(Q);
  
  garch0(P, Q, variance, GARCH, ARCH, K);

  // // initialize to something easy to check
  // for (int i = 0; i < MA.n_rows; i++) {
  // 	MA(i) = i+1;
  // }
  // for (int i = 0; i < AR.n_rows; i++) {
  // 	AR(i) = -(i+1);
  // }
  // for (int i = 0; i < GARCH.n_rows; i++) {
  // 	GARCH(i) = 2*(i+1);
  // }
  // for (int i = 0; i < ARCH.n_rows; i++) {
  // 	ARCH(i) = 3*(i+1);
  // }

  // C
  x[0] = constant;
  lb[0] = -0.5;
  ub[0] =  0.5;

  // AR
  for (int i = 1; i < 1+R; i++) {
	x[i] = AR(i-1);
	lb[i] = -((double) R);
	ub[i] =  ((double) R);
  }
  //cout << "x_AR = " << arr2str(n, x) << endl;

  // MA
  for (int i = 1+R; i < 1+R+M; i++) {
	x[i] = MA(i-(1+R));
	lb[i] = -((double) M);
	ub[i] =  ((double) M);
  }
  //cout << "x_MA = " << arr2str(n, x) << endl;

  // K
  x[1+R+M] = variance;
  lb[1+R+M] = GARCH_TOLERANCE;
  ub[1+R+M] = 1.0;
  //cout << "x_K = " << arr2str(n, x) << endl;

  int i0 = 1+R+M+1;
  
  // GARCH
  for (int i = i0; i < i0+P; i++) {
	x[i] = GARCH(i-i0);
	lb[i] = 0.0;
	ub[i] = 1.0;
  }
  //cout << "x_GARCH = " << arr2str(n, x) << endl;

  i0 = 1+R+M+1+P;
  
  // ARCH
  for (int i = i0; i < i0+Q; i++) {
	x[i] = ARCH(i-i0);
	lb[i] = 0.0;
	ub[i] = 1.0;
  }
  //cout << "x_ARCH = " << arr2str(n, x) << endl;

  // DoF
  x[n-1] = 5.0;
  lb[n-1] =  2.0 + GARCH_TOLERANCE;
  ub[n-1] =  100.0;
  
}

//-------------------------------------------------------------------------------
// Fit the model

garch_struct mygarch::fit_nlopt(const gsl_vector* y, const int R, const int M, const int P, const int Q) {
  
  if (y == NULL)
	assert(0 && "Input vector y is null.");

  int maxRMPQ = 0;

  // get maximum model order
  gsl_vector *orders = gsl_vector_alloc(4);
  gsl_vector_set(orders, 0, R);
  gsl_vector_set(orders, 1, M);
  gsl_vector_set(orders, 2, P);
  gsl_vector_set(orders, 3, Q);
  maxRMPQ = gsl_vector_max(orders);
  gsl_vector_free(orders);

  int iArmaGarchPars = (1 + R + M) + (1 + P + Q);
  int iInnovPars = 1; 
  int iPars  = iArmaGarchPars + iInnovPars;
  int n = iPars;
  
  // <TODO> Free these?
  vector<double>  x(iPars); // parameter vector
  vector<double> lb(iPars); // lower param bounds
  vector<double> ub(iPars); // upper param bounds

  double xa[n], lba[n], uba[n];

  // get initial guesses for params
  mygarch::getInitialPars(n, y, R, M, P, Q, xa, lba, uba);
     
  for (int i = 0; i < n; i++) {
	x[i]  =  xa[i];
	lb[i] = lba[i];
	ub[i] = uba[i];
  }

  // initialize padded series for returns, innovations, and variances
  int iPad = y->size + maxRMPQ;
  gsl_vector *yy = gsl_vector_alloc(iPad);
  gsl_vector *ee = gsl_vector_alloc(iPad);
  gsl_vector *hh = gsl_vector_alloc(iPad);

  // Pad y with maxRMPQ entries of E[y];
  double yMean = gsl_stats_mean(y->data, y->stride, y->size);
  
  for (int i = 0; i < maxRMPQ; i++)
	gsl_vector_set(yy, i, yMean);
  
  for (int i = maxRMPQ; i < iPad; i++) 
	gsl_vector_set(yy, i, gsl_vector_get(y, i-maxRMPQ));

  // pack parameters and variables into structure
  garch_struct s;
  s.n = iPars;
  s.R = R;
  s.M = M;
  s.P = P;
  s.Q = Q;
  s.maxRMPQ = maxRMPQ;
  s.iters = 0;
  s.fxCalls = 0;
  s.yy = yy;
  s.ee = ee;
  s.hh = hh;
	
  nlopt::algorithm algo = nlopt::LN_NELDERMEAD;
  //nlopt::algorithm algo = nlopt::LD_LBFGS;
  //nlopt::algorithm algo = nlopt::LD_MMA;
  nlopt::opt opt(algo, iPars);
    
  opt.set_min_objective(negLLF_nlopt, &s);
  opt.set_xtol_rel(1e-6);
  opt.set_ftol_abs(1e-6);
  //opt.set_maxeval(3000);
  //opt.set_stopval(-9000.0);
  opt.set_lower_bounds(lb);
  opt.set_upper_bounds(ub);

  //double ctol = 1.0e-8;
  //opt.add_inequality_constraint(constraint_garch, NULL, ctol);
  //opt.add_inequality_constraint(constraint_arma_nlc_ar, NULL, ctol);
  //opt.add_inequality_constraint(constraint_arma_nlc_ma, NULL, ctol);
  
  double minf;
  double time = MPI::Wtime();
  nlopt::result nloptResult;
		
  nloptResult = opt.optimize(x, minf);
  
  s.exitStatus = (int) nloptResult;
  s.algo = nlopt::algorithm_name(algo);
  s.negLL = minf;

  gsl_vector* x0 = gsl_vector_alloc(n);
  gsl_vector* x1 = gsl_vector_alloc(n);

  for (int i = 0; i < n; i++) {
	gsl_vector_set(x0, i, xa[i]);
	gsl_vector_set(x1, i, x[i]);
  }

  s.x0 = x0;
  s.x1 = x1;
  s.u = garch_struct_res(s);
  s.runTime = MPI::Wtime() - time;

  // <TODO> x0 and x1 might be orphaned pointers. Should they be freed?
  
  // try {
  // } catch(std::exception& ex) {
  // 	cerr << "Exception caught: " << ex.what() << endl;
  // } catch(...) {
  // 	cerr << "Unknown exception caught" << endl;
  // }

  /*  
  delete [] xa;
  delete [] lba;
  delete [] uba;
  */
  
  return s;
}




//-------------------------------------------------------------------------------
// forecast next return with ARMA-GARCH given standardized residual's in "resid"

vec mygarch::forecast(garch_struct s, vec resid) {

  // Set parameters from x[] to named values.
  int R = s.R;
  int M = s.M;
  int P = s.P;
  int Q = s.Q;

  // get named parameters from parameter vector of estimates (x1)
  vec AR, MA, GARCH, ARCH;
  
  if (R > 0)  AR    = zeros(R);
  if (M > 0)  MA    = zeros(M);
  if (P > 0)  GARCH = zeros(P);
  if (Q > 0)  ARCH  = zeros(Q);

  double C = s.x1->data[0];

  for (int i = 1; i < 1+R; i++) 
  	AR(i-1) = s.x1->data[i];

  int i0 = 1+R;

  for (int i = i0; i < i0+M; i++)
  	MA(i-i0) = s.x1->data[i];

  double K = s.x1->data[1+R+M];

  i0 = 1+R+M+1;
  
  for (int i = i0; i < i0+P; i++)
  	GARCH(i-i0) = s.x1->data[i];

  i0 = 1+R+M+1+P;
  
  for (int i = i0; i < i0+Q; i++)
  	ARCH(i-i0) = s.x1->data[i];

  //double DoF = s.x1->data[(s.n)-1];

  //
  // Forecast one period ahead
  //

  int iPad = s.yy->size; // length of padded series
  int iObs = iPad - s.maxRMPQ; // length of series
  int end = iObs - 1; // end element of the next three vectors
  
  vec innov(iObs);
  vec sigmas(iObs);
  vec series(iObs);
  
  for (int i = s.maxRMPQ; i < iPad; i++) {
  	int idx = i - s.maxRMPQ;
  	innov(idx) = gsl_vector_get(s.ee, i);
  	sigmas(idx) = gsl_vector_get(s.hh, i);
	series(idx) = gsl_vector_get(s.yy, i);
  }

  double nextSigma = std::sqrt(K
							   + sum(ARCH  % pow(innov.rows(end-Q+1,end), 2))
							   + sum(GARCH % sigmas.rows(end-P+1,end)));

  double nextMean = C + sum(AR % series.rows(end-R+1,end)) + sum(MA % innov.rows(end-M+1,end));

  vec nextRet = nextMean + nextSigma*resid;

  //printf("nextMean = %g, nextSigma = %g\n", nextMean, nextSigma);
 
  return nextRet;
}




//-------------------------------------------------------------------------------
// Get the MPI message size for model

int mygarch::getMessageSize(garch_struct s) {
  int iCoeffs = s.n;
  return iNonCoeffs + iCoeffs + 3*s.maxRMPQ + s.u->size;
}



//-------------------------------------------------------------------------------
// Convert structure to array, keeping only what's necessary for forecasting.
// This array is the main MPI message buffer, so it should be as small as possible.

vec mygarch::garch_struct_vec(garch_struct s) {

  int maxRMPQ = s.maxRMPQ; 
  int iL = mygarch::getMessageSize(s);

  vec x(iL);

  x(0) = (double) s.n;
  x(1) = (double) s.R;
  x(2) = (double) s.M;
  x(3) = (double) s.P;
  x(4) = (double) s.Q;
  x(5) = (double) s.maxRMPQ;
  x(6) = s.negLL;
  x(7) = (double) s.exitStatus;
  x(8) = (double) s.fxCalls;
  x(9) = (double) s.iters;
  x(10) = s.runTime;

  int iT = s.yy->size;
  int c = iNonCoeffs;
	
  for (int i = 0; i < s.n; i++) {
	x(c) = gsl_vector_get(s.x1, i);
	c++;
  }
  for (int i = iT-maxRMPQ; i < iT; i++) {
	x(c) = gsl_vector_get(s.ee, i);	
	c++;
  }
  for (int i = iT-maxRMPQ; i < iT; i++) {
	x(c) = gsl_vector_get(s.hh, i);	
	c++;
  }
  for (int i = iT-maxRMPQ; i < iT; i++) {
	x(c) = gsl_vector_get(s.yy, i);
	c++;
  }

  // include standardized residuals (all of them)
  for(int i = 0; i < s.u->size; i++) {
	x(c) = gsl_vector_get(s.u, i);
	c++;
  }

  return x;
}


//-------------------------------------------------------------------------------
// forecast next return with ARMA-GARCH given standardized residual's in "resid"

void mygarch::forecast_fromVec(vec x, vec resid,
							   vec &nextRet, double &nextMean, double &nextSigma) {

  //int n = x(0);
  int R = x(1);
  int M = x(2);
  int P = x(3);
  int Q = x(4);
  int maxRMPQ = x(5);

  vec AR, MA, GARCH, ARCH;
  
  if (R > 0)  AR    = zeros(R);
  if (M > 0)  MA    = zeros(M);
  if (P > 0)  GARCH = zeros(P);
  if (Q > 0)  ARCH  = zeros(Q);

  vec innov(maxRMPQ);
  vec sigmas(maxRMPQ);
  vec series(maxRMPQ);


  // get named ARMA-GARCH parameters
  int i0 = mygarch::iNonCoeffs; // starting idx of coefficients
  
  double C = x(i0);
  
  i0 = i0+1;
  
  for (int i = i0; i < i0+R; i++) 
  	AR(i-i0) = x(i);

  i0 = i0+R;

  for (int i = i0; i < i0+M; i++)
  	MA(i-i0) = x(i);

  double K = x(i0+M);

  i0 = i0+M+1;
  
  for (int i = i0; i < i0+P; i++)
  	GARCH(i-i0) = x(i);

  i0 = i0+P;
  
  for (int i = i0; i < i0+Q; i++)
  	ARCH(i-i0) = x(i);

  //double DoF = x(i0+Q);

  i0 = i0+Q+1;

  for (int i = i0; i < i0+maxRMPQ; i++)
	innov(i-i0) = x(i);

  i0 = i0 + maxRMPQ;
  
  for (int i = i0; i < i0+maxRMPQ; i++)
	sigmas(i-i0) = x(i);

  i0 = i0 + maxRMPQ;
  
  for (int i = i0; i < i0+maxRMPQ; i++)
	series(i-i0) = x(i);

  int end = maxRMPQ - 1; // end element idx
  
  nextSigma = std::sqrt(K
						+ sum(ARCH  % pow(innov.rows(end-Q+1,end), 2))
						+ sum(GARCH % sigmas.rows(end-P+1,end)));

  nextMean = C
	+ sum(AR % series.rows(end-R+1,end))
	+ sum(MA % innov.rows(end-M+1,end));

  nextRet = nextMean + nextSigma*resid;

  //printf("nextMean = %g, nextSigma = %g\n", nextMean, nextSigma);
  //return nextRet;
}


//-------------------------------------------------------------------------------
// test for above method

void struct2array_test(garch_struct s) {

  mygarch::garch_struct_print(s);

  vec arr = mygarch::garch_struct_vec(s);

  // arr.t().print("arr = ");

  cout << "msg_size = " << arr.n_rows << endl;
  
  int iT = s.yy->size;
  int i = iT - 2;
  
  cout << "\n" << gsl_vector_get(s.ee,i) << ", " << gsl_vector_get(s.hh,i) 
	   << ", " << gsl_vector_get(s.yy,i) << "\n";

  i = iT - 1;
  
  cout << gsl_vector_get(s.ee,i) << ", " << gsl_vector_get(s.hh,i) 
	   << ", " << gsl_vector_get(s.yy,i) << "\n\n";

  vec resid, f1;

  resid << 0.1 << endr
		<< 0.23 << endr;

  double nextMean, nextSigma;
  
  mygarch::forecast_fromVec(arr, resid, f1, nextMean, nextSigma);

  f1.print("f1 = ");

  vec f2 = mygarch::forecast(s, resid);

  f2.print("f2 = ");
}


//-------------------------------------------------------------------------------
// <main>

// int main() {

//   string fname = "/Users/jimmiegoode/Dropbox/Projects/Glimm/data/scalingSeries/1.csv";
 
//   gsl_matrix* mnRet = IO::importCSVmatrix(fname);
//   gsl_vector* vnRet = gsl_vector_alloc(mnRet->size1);
//   gsl_matrix_get_col(vnRet, mnRet, 0);

//   garch_struct g = mygarch::fit_nlopt(vnRet, 1, 1, 1, 1);
//   // mygarch::garch_struct_print(g);
//   struct2array_test(g);

//   return 0;
  
//   DistPars d = myskewt::skewedstudenttfit_bfgs(COLD, g.u);
//   myskewt::DistPars_print(d);
  
//   int nSim = 10;
//   vec gamma, mu;
//   double df = 5.0;
//   mat C;

//   gamma << 0.1 << endr
//   		<< -0.1 << endr;

//   mu << 0.01 << endr
//   	 << 0.02 << endr;
 
//   C << 1.5 << 0.1 << endr
//   	<< 0.1 << 2.0 << endr;

//   Stats::gsl_rng_init(1346788577);
//   mat resid = myskewt::mvskewtrnd(gamma, mu, df, C, nSim);
//   vec next = mygarch::forecast(g, resid.col(0));
  
//   return 0;

//   //Stats::mvnrnd_test();
//   //Stats::mvskewtrnd_test();
//   //myskewt::mvskewtrnd_test();

//   /*
//   double t = MPI::Wtime();
//   Stats::gsl_rng_init(0);
//   test_smallWorld();
//   cout << "\nElapsed Time = " << MPI::Wtime() - t << " seconds.\n\n";
//   */
  
//   return 0;
  
  
//   ////test_llf(vnRet, 285);
//   ////test_scaling();
//   //cout << "DT_test_1 = " << MPI::Wtime() - t << endl;

//   // Moments mom = xMomentSummary(s.u);
//   // printMoments(mom);
  
//   /*
//   int R, M;
//   double constant, variance;
  
//   R = 0; M = 3;
//   vec AR(R);
//   vec MA(M);
  
//   arma0(vnRet, R, M, constant, variance, AR, MA);

//   AR.print("AR = ");
//   MA.print("MA = ");
//   cout << "constant = " << constant << endl;
//   cout << "variance = " << variance << endl;
//   */

//   /*
//   // getInitialPars(const int n, const gsl_vector* y,
//   // 							 const int R, const int M, const int P, const int Q,
//   // 							 double x[], double lb[], double ub[]
  
//   int R = 1;
//   int M = 1;
//   int P = 1;
//   int Q = 1;
//   int n = (1 + R + M) + (1 + P + Q) + 1;

//   double x[n];
//   double lb[n];
//   double ub[n];
  
//   mygarch::getInitialPars(n, vnRet, R, M, P, Q, x, lb, ub);
//   cout << "x_0 = " << arr2str(n, x) << endl;
//   cout << "lb  = " << arr2str(n, lb) << endl;
//   cout << "ub  = " << arr2str(n, ub) << endl;
//   */
  
//   //double t;
//   //t = clock();
//   //cout << "DT_test_algo = " << ((float)t)/CLOCKS_PER_SEC << endl;

// }
