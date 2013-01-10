/*
  Author: Jimmie Goode
  Created: 2012-12-28
*/

#include "myskewt_lut.hpp"

myskewt_lut::myskewt_lut() {
  //cout << "constructed" << endl;
}

myskewt_lut::~myskewt_lut() {
}


//-------------------------------------------------------------------------------
// PDF of the first component of \widehat{X}

vec lut_pdf(double beta1, double beta2, double df, int d, int nSim, vec xGrid, mat A) {
  
  vec mu(d);
  mu.zeros();
 
  vec beta(d);
  beta.zeros();

  beta(0) = beta1;
  beta(1) = beta2;

  //mat X = myskewt::mvskewtrnd_1(beta, mu, df, Sigma, nSim);
  mat X = myskewt::mvskewtrnd_2(beta, mu, df, nSim, A);
  vec x1 = X.col(0);

  vec fx = Stats::ksdensity(x1, xGrid);

  return fx;
}


//-------------------------------------------------------------------------------
// Build a lookup table

void myskewt_lut::buildTable(vec beta1, vec beta2, double df, int d, int nSim, vec xGrid) {

  this->beta1 = beta1;
  this->beta2 = beta2;
  this->df = df;
  this->d = d;
  this->nSim = nSim;
  this->xGrid = xGrid;
  
  int ib1 = beta1.n_rows;
  int ib2 = beta2.n_rows;
  int ig  = xGrid.n_rows;

  mat tab(ib2*ig, ib1);

  cout << "Computing cholesky...";
  
  mat Sigma = ((df-2.0)/df) * eye(d,d);
  mat A = chol(Sigma);

  cout << "Done." << endl;
  
  for (int i = 0; i < ib1; i++) {

	int i0 = 0;
	double nBeta1 = beta1(i);

	printf("pct = %3.2f, beta1 = %3.4f\n", ((double) i)/((double) ib1)*100.0, nBeta1);

	for (int j = 0; j < ib2; j++) {

	  double nBeta2 = beta2(j);

	  vec f = lut_pdf(nBeta1, nBeta2, df, d, nSim, xGrid, A);
	  //vec f(ig);
	  //f = ones(ig) * (i0 + i);

	  tab(span(i0, i0+ig-1), i) = f;
	  i0 = i0+ig;
	}
  }
  this->table = tab;
}

//-------------------------------------------------------------------------------

string arr2str(vec a) {
  return "linspace("
	+ boost::lexical_cast<string>(a(0)) + ", "
	+ boost::lexical_cast<string>(a(a.n_rows-1)) + ", "
	+ boost::lexical_cast<string>(a.n_rows)
	+ ")";
}

//-------------------------------------------------------------------------------

void myskewt_lut::save(string path, file_type type, string ext) {

  //file_type type = csv_ascii;
  //string ext = ".csv";
  
  string info = "df = " + boost::lexical_cast<string>(df) + "\n"
	+ "d = " + boost::lexical_cast<string>(d) + "\n"
	+ "nSim = " + boost::lexical_cast<string>(nSim) + "\n"
	+ "xGrid = " + arr2str(xGrid) + "\n"
	+ "beta1 = " + arr2str(beta1) + "\n"
	+ "beta2 = " + arr2str(beta2);
  
  system(("mkdir " + path).c_str());

  FILE *fout;
  fout = fopen((path + "/info.txt").c_str(), "w+");
  fprintf(fout, "%s\n", info.c_str());
  fclose(fout);

  beta1.save(path + "/beta1" + ext, type);
  beta2.save(path + "/beta2" + ext, type);
  xGrid.save(path + "/xGrid" + ext, type);
  table.save(path + "/table" + ext, type);

  cout << "Saved LUT data to: " << path << endl;
}

//-------------------------------------------------------------------------------

void myskewt_lut::load(string folder, file_type type, string ext) {

  //file_type type = csv_ascii;
  //string ext = ".csv";

  cout << "Loading LUT from: '" << folder << "'...";
  
  // these aren't needed, set to null values
  df = 0;
  d  = 0;
  nSim = 0;

  // load arrays
  beta1.load(folder + "/beta1" + ext, type);
  beta2.load(folder + "/beta2" + ext, type);
  xGrid.load(folder + "/xGrid" + ext, type);
  table.load(folder + "/table" + ext, type);

  cout << "Done." << endl;
}

//-------------------------------------------------------------------------------
// print out the LUT members

void myskewt_lut::print() {
  beta1.t().print("beta1 = ");
  beta2.t().print("beta2 = ");
  xGrid.t().print("xGrid = ");
  table.print("Table = ");
}

//-------------------------------------------------------------------------------
// Get the bin index of x in the bins defined by midpoints in sorted vector X.

int findPoint(vec X, double x) {

  int n = X.n_rows;
  double inc = (X(1) - X(0))/2.0;

  if (x < X(0) - inc)
	return -1;
  
  if (x > X(n-1) + inc)
	return -2;
  
  for (int i = 0; i < n-1; i++) {
	if (x > X(i) && x <= X(i+1)) {
	  if ((x > X(i) - inc) && (x <= X(i) + inc))
		return i;
	  else
		return i+1;
	}
  } 
  return -3;
}

//-------------------------------------------------------------------------------
// find the density values for interpolation

vec myskewt_lut::lookupPDF(double nBeta1, double nBeta2) {

  int ib1 = findPoint(beta1, nBeta1);
  int ib2 = findPoint(beta2, nBeta2);
  int ig = xGrid.n_rows;
  
  int iCol = ib1;
  int i0 = ig*ib2; // starting row index for ib2 in column ib1

  vec f = table(span(i0, i0+ig-1), iCol);

  if (0) {
	cout << "-------------------------------------------------------------------------------\n";
	cout << "nBeta1 = " << nBeta1 << ", nBeta2 = " << nBeta2 << endl;
	cout << "ib1 (col idx) = " << ib1 << ", ib2 = " << ib2 << endl;
	f.t().print("f_lut = ");
	cout << "-------------------------------------------------------------------------------\n";
  }
  return f;
}


//-------------------------------------------------------------------------------
// VaR from lookup table

void myskewt_lut::lut_var(vec w, vec beta, vec mu, double df, mat A, double epsilon,
						  double &VaR, double &fVaR) {

  if (A(0, A.n_rows-1) != 0.0)
	assert(0 && "A should be lower triangular so that mnCov = A*A'");

  int d = w.n_rows;
  assert(d == mu.n_rows && d == beta.n_rows);

  mat Ainv     = inv(A);
  vec wHat     = A.t() * w;
  double normW = norm(wHat,2);
  vec betaHat  = Ainv * beta;
  vec e1       = wHat / normW;
  double b1    = dot(betaHat, e1);
  vec betaPara = b1 * e1;            // parallel part
  vec betaPerp = betaHat - betaPara; // perpendicular part

  lut_solveG(epsilon, betaPara(0), betaPerp(1), fVaR);
  VaR = normW*fVaR + dot(w,mu);
}

//-------------------------------------------------------------------------------
// Structure of paramters for evaluated the G(x) function.

struct G_params
{
  double epsilon;
  gsl_interp_accel *acc;
  gsl_spline *spline;
  vec xGrid;
  vec fx;
};

//-------------------------------------------------------------------------------
// The integrand corresponds to the tabulated PDF, f^(1)(X).

double integrand(double x, void *params) {
  
  struct G_params *p = (struct G_params *) params;

  double fxi;

  if ((x < p->xGrid(0)) || (x > p->xGrid(p->xGrid.n_rows-1)))
	fxi = 0.0;
  else
	fxi = gsl_spline_eval(p->spline, x, p->acc);

  //printf("f(%3.4f) = %3.4f\n", x, fxi);
  
  return fxi;
}

//-------------------------------------------------------------------------------
// G(x), returns (integral_{-inf}^{x} f(x) - epsilon)

double G(double x, void *params) {
  
  struct G_params *p = (struct G_params *) params;

  int numIntervals = 1000;
  
  gsl_integration_workspace *w = gsl_integration_workspace_alloc(numIntervals);
       
  double result, error;
  
  gsl_function F;
  F.function = &integrand;
  F.params = p;

  double lb = -10.0; // this should theoretically be -Inf
  double ub = x;
  
  if (ub <= lb)
	result = 0.0;
  else
	gsl_integration_qags(&F, lb, ub, 0, 1e-7, numIntervals, w, &result, &error); 

  //printf("integral(%3.4f, %3.4f) = % .18f\n", lb, ub, result);
     
  gsl_integration_workspace_free(w);

  return (result - p->epsilon);

}

//-------------------------------------------------------------------------------
// Solve G(fVaR) == 0

void myskewt_lut::lut_solveG(double epsilon, double nBeta1, double nBeta2, double &fVaR) {
  
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double r = 0;
  double x_lo = -9.0, x_hi = 9.0; // bounds for finding root
  gsl_function F;
  int verbose = 0;
  
  // Setp interpolation for evaluating f(x)
  vec fx =  myskewt_lut::lookupPDF(nBeta1, nBeta2);
   
  int iGrid = fx.n_rows;
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, iGrid);
  gsl_spline_init (spline, xGrid.memptr(), fx.memptr(), iGrid);
  
  G_params params;
  params.epsilon = epsilon;
  params.acc = acc;
  params.spline = spline;
  params.xGrid = xGrid;
  params.fx = fx;
     
  F.function = &G;
  F.params = &params;
     
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc (T);
  gsl_root_fsolver_set (s, &F, x_lo, x_hi);

  if (verbose) {
	printf ("using %s method\n", gsl_root_fsolver_name(s));
	printf ("%5s [%9s, %9s] %9s\n", "iter", "lower", "upper", "root");
  }

  do {
	iter++;
	status = gsl_root_fsolver_iterate (s);
	r = gsl_root_fsolver_root (s);
	x_lo = gsl_root_fsolver_x_lower (s);
	x_hi = gsl_root_fsolver_x_upper (s);
	status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
     
	if (status == GSL_SUCCESS) {
	  if (verbose)
		printf ("Converged:\n");
	}
    if (verbose) 
	  printf ("%5d [%.7f, %.7f] %.7f\n", iter, x_lo, x_hi, r);
	
  }
  while (status == GSL_CONTINUE && iter < max_iter);
     
  gsl_root_fsolver_free (s);

  
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);

  fVaR = r;
  //return status;
}



//-------------------------------------------------------------------------------

void myskewt_lut::tester(int c) {

  cout << "\n---> myskewt_lut tester: case " << c << " <---\n\n";
  
  switch(c) {
  case 1:
	{
	  // Test lookup table generation
	  
	  vec beta1 = linspace(-1,1,3);
	  vec beta2 = linspace(-1,1,3);
	  double df = 5.0;
      double d = 30;
	  int nSim = 10000;

	  vec xGrid = linspace(-10,10,5);
	  
	  myskewt_lut lut;
	  
	  lut.buildTable(beta1, beta2, df, d, nSim, xGrid);

	  lut.print();

	  // lut.save("../LUTs");

	  double nBeta1 = 0.6; //0.12813;
	  double nBeta2 = -0.6; //-0.44528;
	  
	  lut.lookupPDF(nBeta1, nBeta2);
	  
	}
	break;
  case 2:
	{
	  // Test findPoint
	  
	  vec xGrid = linspace(-10,10,5);
	  double x = 2.6;

	  int idx = findPoint(xGrid, x);
	  
	  xGrid.t().print("xGrid = ");
	  cout << "x = " << x << endl;
	  cout << "idx = " << idx << endl;
	  
	}
	break;
	
  case 3:
	{
	  // Test lut_solveG
	  
	  myskewt_lut lut;

	  lut.load("../LUTs/df=5_d=30_nSim=10000", csv_ascii, ".csv");
	  //lut.print();

	  double epsilon = 0.01;
	  //double nBeta1 = 0.12813;
	  //double nBeta2 = -0.44528;
	  //double fVaR;
	  //lut.lut_solveG(epsilon, nBeta1, nBeta2, fVaR);

	  vec w;
	  vec beta;
	  vec mu;
	  double df = 5;
	  mat A;
	  string p = "/Users/jimmiegoode/Documents/Glimm/LookupTables/";

	  w.load(p+"e_w.csv",csv_ascii);
	  beta.load(p+"e_beta.csv",csv_ascii);
	  mu.load(p+"e_mu.csv",csv_ascii);
	  A.load(p+"e_A.csv",csv_ascii);

	  double VaR, fVaR;
	  lut.lut_var(w, beta, mu, df, A, epsilon, VaR, fVaR);

	  cout << ">>" << VaR << ", " << fVaR << endl;

	}
	break;

  case 4:
	{
	  // test integration
	  
	  int i;
	  double xi, yi, x[10], y[10];
     
	  printf ("#m=0,S=2\n");
     
	  for (i = 0; i < 10; i++) {
		x[i] = i + 0.5 * sin (i);
		y[i] = i + cos (i * i);
		printf ("%g %g\n", x[i], y[i]);
	  }
     
	  printf ("#m=1,S=0\n");
     
	  {
		gsl_interp_accel *acc = gsl_interp_accel_alloc ();
		gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, 10);
     
		gsl_spline_init (spline, x, y, 10);
     
		for (xi = x[0]; xi < x[9]; xi += 0.01)
		  {
			yi = gsl_spline_eval (spline, xi, acc);
			printf ("%g %g\n", xi, yi);
		  }
		gsl_spline_free (spline);
		gsl_interp_accel_free (acc);
	  }
	}
	break;
	
  case 5:
	{
	  // test generation of table
	  myskewt_lut lut;

	  vec beta1 = linspace(-1,1,30);
	  vec beta2 = linspace(-1,1,30);
	  double df = 5.0;
      double d = 30;
	  int nSim = 10000;
	  vec xGrid = linspace(-10,10,100);
	  
	  lut.buildTable(beta1, beta2, df, d, nSim, xGrid);

	}
	break;

  case 6:
	{
	  vec x = linspace(-10,10,5);
	  x.save("testFile.bin", arma_binary);

	  vec y;
	  y.load("testFile.bin", arma_binary);

	  x.print("x = ");
	  y.print("y = ");
	  
	}
   
	
  } //switch
}

//-------------------------------------------------------------------------------
// <main>

// int main(int argc, char** argv) {
  
//   if (argc < 2)
// 	assert(0 && "Invalid Arguments");
  
//   int testCase = atoi(argv[1]);

//   myskewt_lut::tester(testCase);
  
//   return 0;
// }
