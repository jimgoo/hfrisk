#include <stdio.h>
#include <vector>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <mpi.h>
#include <armadillo>
#include "elemental.hpp"

#include "IO.hpp"
#include "Stats.hpp"
#include "mygarch.hpp"
#include "myskewt.hpp"
#include "Constants.h"
#include "mystable.hpp"
#include "myskewt_lut.hpp"

using namespace arma;
using namespace elem;

typedef double R;
typedef Complex<R> C;

//-------------------------------------------------------------------------------

mat mnGarch;

// string garchFile = "/Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/csi/patentData/csi_20030101_20120801_v3/logret.abin";
string garchFile = "/Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/csi/patentData/exports/20130105_ST/mnResults/20080103.h5";

int gRows, gCols, vnGarchSize;
double *vnGarch;

// GARCH estimation parameters
int oR = 1, oM = 1, oP = 1, oQ = 1, maxRMPQ; // orders
int iL = 1008; // lookback
int i0;

double gammaScale = 0.01;

int depStruc = 1; // {1 = MSST, 2 = ASSG}
int nSim = 10000; // number of MC samples

struct RunTimes
{
  double t_dep_est_mc;  // time to estimate MC dependence structure
  double t_dep_est_lut; // time to estimate LUT dependence structure
  double t_shrink;      // time to shrink
  double t_chol_mc;     // time for MC cholesky
  double t_chol_lut;    // time for LUT cholesky
  double t_VaR_mc;      // time to estimate VaR with MC given all required parameters
  double t_VaR_lut;     // time to estimate VaR with LUT given all required parameters
  double t_total_mc;    // time for everything with MC
  double t_total_lut;   // time for everything with LUT
};

RunTimes rt;


//-------------------------------------------------------------------------------

static void printRT(RunTimes rt) {
  
  cout << "t_dep_est_mc  = " << rt.t_dep_est_mc  << endl
	   << "t_dep_est_lut = " << rt.t_dep_est_lut << endl
	   << "t_chol_mc     = " << rt.t_chol_mc     << endl
	   << "t_chol_lut    = " << rt.t_chol_lut    << endl
	   << "t_VaR_mc      = " << rt.t_VaR_mc      << endl
	   << "t_VaR_lut     = " << rt.t_VaR_lut     << endl
	   << "t_total_mc    = " << rt.t_total_mc    << endl
	   << "t_total_lut   = " << rt.t_total_lut   << endl;
}

//-------------------------------------------------------------------------------

int main(int argc, char* argv[]) {
  
  Initialize(argc, argv);

  mpi::Comm comm = mpi::COMM_WORLD;
  const int commRank = mpi::CommRank(comm);

  try {
	
	// parse args
	for (int i = 1; i < argc-1; i++)
	  if (strncmp(argv[i],"-garchFile",20) == 0)
		garchFile = boost::lexical_cast<string>(argv[i+1]);

	// rank 0 setup
	if (commRank == 0) {

	  // Test if d=100,000 will fit in memory:
	  /*
	  mnGarch.load(garchFile, arma_binary);
	  cout << "size(mnGarch_0) = " << mnGarch.n_rows << " X " << mnGarch.n_cols << endl;
	  mnGarch = repmat(mnGarch, 1, 19);
	  cout << "size(mnGarch_1) = " << mnGarch.n_rows << " X " << mnGarch.n_cols << endl;
	  sleep(1000);
	  */

	  mnGarch.load(garchFile, hdf5_binary);
	  //mnGarch = mnGarch(span(0,gRows-1), span(0,gCols-1));
	  
	  gRows = mnGarch.n_rows;
	  gCols = mnGarch.n_cols;
	  vnGarchSize = gRows*gCols;

	  // extract Garch parameters and residuals

	  // set the number of elements in the MPI marginal message 
	  vec orders;
	  orders << oR << endr << oM << endr << oP << endr << oQ;
	  maxRMPQ = max(orders);
	  
	  int iCoeffs = (int) mnGarch(0,0);
	  int iNonCoeffs = mygarch::iNonCoeffs;
	  i0 = iNonCoeffs + iCoeffs + 3*maxRMPQ;
		
	} // commRank==0

	MPI_Bcast(&i0, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&iL, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&gRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&gCols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&vnGarchSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// allocate array on ALL processes
	vnGarch = new double[vnGarchSize];

	// fill array on the root process
	if (commRank == 0) { 						   
	  int idx = 0;
	  for (int c = 0; c < mnGarch.n_cols; c++) {
		for (int r = 0; r < mnGarch.n_rows; r++) {
		  vnGarch[idx] = mnGarch(r,c);
		  idx++;
		}
	  }
	  cout << "----> Broadcasting...\n";
	}

	MPI_Bcast(vnGarch, vnGarchSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
	// cout << "----> rank " << commRank
	// 	 << ", vnGarch[0] = " << vnGarch[0]
	// 	 << ", vnGarchSize = " << vnGarchSize
	// 	 << ", i0 = " << i0 << endl;

	Grid g(comm);
	DistMatrix<R> X(iL, gCols, g);
	DistMatrix<R> garchPars(i0, gCols, g);

	// fill X with data (residuals of GARCH)
	for (int j = 0; j < gCols; j++)
	  for (int i = i0; i <= gRows-4-1; i++)
		X.Set(i, j, vnGarch[i + j*gRows]);

	// fill garchPars with data
	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i <= i0-1; i++)
		garchPars.Set(i, j, vnGarch[i + j*gRows]);
		
	//X.Print("X = ");
	//garchPars.Print("garchPars = ");

	// set vectors gamma and mu
	int idxGamma = gRows - 4;
	int idxMu = gRows - 3;
	double df = 5.0;

	DistMatrix<R> stGamma(gCols, 1, g);
	DistMatrix<R> stMu(gCols, 1, g);

	// fill stGamma and stMu, scaling gamma as necessary
	for (int i = 0; i <= gCols; i++) {
	  stGamma.Set(i, 0, gammaScale * vnGarch[idxGamma + i*gRows]);
	  stMu.Set(i, 0, vnGarch[idxMu + i*gRows]);
	}
	
	
	//-------------------------------------------------------------------------------
	// Shrinkage
	//-------------------------------------------------------------------------------

	MPI_Barrier(MPI_COMM_WORLD);
	if (commRank == 0) {
	  rt.t_dep_est_mc = MPI::Wtime();
	}
	
	// get column sums
	DistMatrix<R> colSum(gCols, 1, g);
	Zeros(gCols, 1, colSum);
	
	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i < iL; i++)
		colSum.Set(j, 0, colSum.Get(j,0) + X.Get(i,j));

	// subtract column means
	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i < iL; i++)
		X.Set(i, j, X.Get(i,j) - colSum.Get(j,0)/iL);
	
	// Compute sample covariance
	//	(1/n) * X.t() * X
	DistMatrix<R> sample(gCols, gCols, g);
  
	double alpha = 1.0/((double) iL);
	
	Gemm(TRANSPOSE, NORMAL, alpha, X, X, (double)0, sample);

	//sample.Print("sample_dist = ");

	assert(sample.Width() == sample.Height()); // check square
	assert(sample.Width() == gCols);
	
	// Construct the prior:
	//	meanvar*mnEye + meancov*(ones(n,n) - diagmat(ones(n)));
	double meanvar = 0.0, meancov = 0.0, d = (double)gCols;

	for (int i = 0; i < gCols; i++)
	  for (int j = 0; j < gCols; j++)
		if (i == j)
		  meanvar += sample.Get(i,j);
		else
		  meancov += sample.Get(i,j);
	
	meanvar = meanvar/d;
	meancov = meancov/d/(d-1);

	//cout << "----> meanvar = " << meanvar << ", meancov = " << meancov << endl;

	// meanvar*eye
	DistMatrix<R> tmp1(gCols, gCols, g);
	Identity(gCols, gCols, tmp1);
	
	for (int i = 0; i < gCols; i++)
	  tmp1.Set(i, i, meanvar);

	// meancov*(ones(n,n) - diagmat(ones(n)))
	// This is zero on the diagonal and meancov everywhere else.
	DistMatrix<R> prior(gCols, gCols, g);
	
	for (int i = 0; i < gCols; i++)
	  for (int j = 0; j < gCols; j++)
		if (i == j)
		  prior.Set(i, j, 0.0);
		else
		  prior.Set(i, j, meancov);
	
	DistMatrix<R> eye(gCols, gCols, g);
	Identity(gCols, gCols, eye);

	// prior = tmp1 + prior
	Gemm(NORMAL, NORMAL, (double)1, tmp1, eye, (double)1, prior);

	//prior.Print("Prior = ");
	
	//tmp1 =  sample - prior
	tmp1 = prior;
    Gemm(NORMAL, NORMAL, (double)1, sample, eye, (double)-1, tmp1);

	const double c = std::pow(Norm(tmp1, FROBENIUS_NORM), 2.0);

	DistMatrix<R> Y(iL, gCols, g);
	for (int i = 0; i < iL; i++)
	  for (int j = 0; j < gCols; j++)
		Y.Set(i, j, std::pow(X.Get(i,j), 2.0));

	// tmp1 = Y' * Y
	Gemm(TRANSPOSE, NORMAL, (double)1, Y, Y, (double)0, tmp1);

	assert(tmp1.Width() == gCols && tmp1.Height() == gCols); // square?
	
	double Ysum = 0.0, samSum = 0.0;
	for (int i = 0; i < gCols; i++) {
	  for (int j = 0; j < gCols; j++) {
		Ysum += tmp1.Get(i, j);
		samSum += std::pow(sample.Get(i, j), 2.0);
	  }
	}

	const double p = (1.0/((double)iL))*Ysum - samSum;
	
	// Compute covariance of Y
	Zeros(gCols, 1, colSum);

	// compute column sums for the means
	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i < iL; i++)
		colSum.Set(j, 0, colSum.Get(j,0) + Y.Get(i,j));

	// subtract column means
	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i < iL; i++)
		Y.Set(i, j, Y.Get(i,j) - colSum.Get(j,0)/((double)iL));

	// tmp1 = 1/t * Y' * Y
    alpha = 1.0/((double)iL - 1.0);
	
	Gemm(TRANSPOSE, NORMAL, alpha, Y, Y, (double)0, tmp1);

	//tmp1.Print("cov(Y) = ");

	double rDiag = 0.0;
	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i < gCols; i++)
		rDiag += tmp1.Get(i,j);

	rDiag = (1.0/((double)gCols))*rDiag;

	const double k = (p - rDiag)/c;
	const double shrinkage = std::max(0.0, std::min(1.0, k/((double)iL)));

	// shrinkage*prior + (1-shrinkage)*sample
	Gemm(NORMAL, NORMAL, shrinkage, prior, eye, 1.0-shrinkage, sample);

	//sample.Print("Final = ");

	MPI_Barrier(MPI_COMM_WORLD);
	if (commRank == 0) {

	  rt.t_shrink = MPI::Wtime() - rt.t_shrink;
	  
	  //cout << "----> rt.t_shrink elem = " << rt.t_shrink << endl;
	  // t_shrink = MPI::Wtime();
	  // double ds;
	  // mat aS = Stats::cov2para(mnGarch, ds);
	  // t_shrink = MPI::Wtime() - t_shrink;
	  // cout << "----> t_shrink serial = " << t_shrink << endl;
	  //aS.print("armaFinal = ");
	}

	//-------------------------------------------------------------------------------
	// Cholesky
		
	//Gemm(NORMAL, NORMAL, ((df - 2.0)/df), sample, eye, (2.0*df)/((df - 2.0)*(df - 4.0)),
	double c1 = -2.0*df/((df - 2.0)*(df - 4.0));
	double c2 = ((df - 2.0)/df);
	Gemm(NORMAL, TRANSPOSE, c1, stGamma, stGamma, c2, sample);


	// start chol time
	MPI_Barrier(MPI_COMM_WORLD);
	if (commRank == 0) {
	  rt.t_chol_mc = MPI::Wtime();
	}
	
	// Cholesky
	Cholesky(LOWER, sample);

	// stop chol time
	MPI_Barrier(MPI_COMM_WORLD);
	if (commRank == 0) {
	  rt.t_chol_mc = MPI::Wtime() - rt.t_chol_mc;
	  rt.t_dep_est_mc = MPI::Wtime() - rt.t_dep_est_mc;
	}

	//-------------------------------------------------------------------------------
	// Simulate random deviates

	
	//-------------------------------------------------------------------------------
	//

	if (commRank == 0) {
	  
	  printRT(rt);
	  
	}
	
  } catch( ArgException& e ) {
    // do nothing
  } catch( exception& e ) {
	  ostringstream os;
	  os << "Process " << commRank << " caught exception: " << e.what() << endl;
	  cerr << os.str();
#ifndef RELEASE
	  DumpCallStack();
#endif
  }

  Finalize();
  
  return 0;
}



/*
  
Jimmies-MacBook-Pro:src jimmiegoode$ mpirun -np 1 ../bin/HFRisk -garchFile /Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/20080102.h5 
t_dep_est_mc  = 140.428
t_dep_est_lut = 0
t_chol_mc     = 9.65242
t_chol_lut    = 0
t_VaR_mc      = 0
t_VaR_lut     = 0
t_total_mc    = 0
t_total_lut   = 0

Jimmies-MacBook-Pro:src jimmiegoode$ mpirun -np 2 ../bin/HFRisk -garchFile /Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/20080102.h5 
t_dep_est_mc  = 244.037
t_dep_est_lut = 0
t_chol_mc     = 7.47791
t_chol_lut    = 0
t_VaR_mc      = 0
t_VaR_lut     = 0
t_total_mc    = 0
t_total_lut   = 0

*/
