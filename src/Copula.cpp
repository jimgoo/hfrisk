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

string garchFile = "/Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/csi/patentData/csi_20030101_20120801_v3/logret.abin";

int gRows, gCols, vnGarchSize;

//double *vnGarch; <TODO> adjustable sizes for vnGarch
double vnGarch[50];


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
	  
	  mnGarch.load(garchFile, arma_binary);
	  mnGarch = mnGarch(span(0,9), span(0,4));

	  //mnGarch.print("mnGarch = ");
	  //sum(mnGarch, 0).t().print("sum(mnGarch,0) = ");

	  gRows = mnGarch.n_rows;
	  gCols = mnGarch.n_cols;
	  vnGarchSize = gRows*gCols;
	  
	  //cout << "----> mnGarch size = " << gRows << " X " << gCols << endl;
	  //cout << "----> vnGarchSize = " << vnGarchSize << endl;

	  //vnGarch = new double[vnGarchSize]; //<TODO>
						   
	  int idx = 0;
	  for (int c = 0; c < mnGarch.n_cols; c++) {
		for (int r = 0; r < mnGarch.n_rows; r++) {
		  vnGarch[idx] = mnGarch(r,c);
		  idx++;
		}
	  }
	} // commRank==0

	MPI_Bcast(&gRows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&gCols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&vnGarchSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&vnGarch, vnGarchSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	cout << "----> rank " << commRank
		 << ", vnGarch[0] = " << vnGarch[0]
		 << ", vnGarchSize = " << vnGarchSize << endl;
	
	Grid g(comm);
	DistMatrix<R> X(gRows, gCols, g );

	//
	// fill X with data
	//
	
	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i < gRows; i++)
		X.Set(i, j, vnGarch[i + j*gRows]);

	//X.Print("X = ");

	//
	// get column sums
	//
	
	DistMatrix<R> colSum(gCols, 1, g);
	Zeros(gCols, 1, colSum);
	
	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i < gRows; i++)
		colSum.Set(j, 0, colSum.Get(j,0) + X.Get(i,j));

	//colSum.Print("colSum = ");

	//
	// subtract column means
	//

	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i < gRows; i++)
		X.Set(i, j, X.Get(i,j) - colSum.Get(j,0)/gRows);

	
	//X.Print("X_demean = ");
	
	//
	// Compute sample covariance
	//	(1/n) * X.t() * X
	//

	DistMatrix<R> sample(gCols, gCols, g);
  
	double alpha = 1.0/((double) gRows);
	
	Gemm(TRANSPOSE, NORMAL, alpha, X, X, (double)0, sample);

	//sample.Print("sample_dist = ");

	assert(sample.Width() == sample.Height()); // check square
	assert(sample.Width() == gCols);
	
	//
	// Construct the prior
	//	meanvar*mnEye + meancov*(ones(n,n) - diagmat(ones(n)));
	//

	double meanvar = 0.0, meancov = 0.0, d = (double)gCols;

	for (int i = 0; i < gCols; i++)
	  for (int j = 0; j < gCols; j++)
		if (i == j)
		  meanvar += sample.Get(i,j);
		else
		  meancov += sample.Get(i,j);
	
	meanvar = meanvar/d;
	meancov = meancov/d/(d-1);

	cout << "----> meanvar = " << meanvar << ", meancov = " << meancov << endl;

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

    //void Gemv(Orientation orientation, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& x, T beta, DistMatrix<T>& y)

	//void Gemm(Orientation orientationOfA, Orientation orientationOfB, T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, T beta, DistMatrix<T>& C)

	DistMatrix<R> eye(gCols, gCols, g);
	Identity(gCols, gCols, eye);

	// prior = tmp1 + prior
	Gemm(NORMAL, NORMAL, (double)1, tmp1, eye, (double)1, prior);

	//prior.Print("Prior = ");
	
	//tmp1 =  sample - prior
	tmp1 = prior;
    Gemm(NORMAL, NORMAL, (double)1, sample, eye, (double)-1, tmp1);

	//tmp1.Print("sample-prior = ");
	
	const double c = std::pow(Norm(tmp1, FROBENIUS_NORM), 2.0);

	DistMatrix<R> Y(gRows, gCols, g);
	for (int i = 0; i < gRows; i++)
	  for (int j = 0; j < gCols; j++)
		Y.Set(i, j, std::pow(X.Get(i,j), 2.0));

	// tmp1 = Y' * Y
	Gemm(TRANSPOSE, NORMAL, (double)1, Y, Y, (double)0, tmp1);

	assert(tmp1.Width() == gCols && tmp1.Height() == gCols); // square?
	
	double Ysum = 0.0, samSum;
	for (int i = 0; i < gCols; i++) {
	  for (int j = 0; j < gCols; j++) {
		Ysum += tmp1.Get(i, j);
		samSum += std::pow(sample.Get(i, j), 2.0);
	  }
	}

	double p = (1.0/((double)gRows))*Ysum - samSum;
	
	//cout << "p = " << p << endl;

	//
	// Compute covariance of Y
	//


	Y.Print("Y = ");
	
	Zeros(gCols, 1, colSum);

	// compute column sums for the means
	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i < gRows; i++)
		colSum.Set(j, 0, colSum.Get(j,0) + Y.Get(i,j));

	colSum.Print("colSum = ");

	// subtract column means
	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i < gRows; i++)
		Y.Set(i, j, Y.Get(i,j) - colSum.Get(j,0)/((double)gRows));

	// tmp1 = 1/t * Y' * Y
    alpha = 1.0/((double)gRows);
	
	Gemm(TRANSPOSE, NORMAL, alpha, Y, Y, (double)0, tmp1);

	tmp1.Print("cov(Y) = ");

	
	double rDiag = 0.0;
	for (int j = 0; j < gCols; j++)
	  for (int i = 0; i < gRows; i++)
		rDiag += tmp1.Get(i,j);

	rDiag = (1.0/((double)gCols))*rDiag;

	//cout << "rDiag = " << rDiag << endl;
	
	
	//-------------------------------------------------------------------------------
	
	if (commRank == 0) {
 
	  int t = mnGarch.n_rows;
	  int n = mnGarch.n_cols;

	  vec meanx = arma::mean(mnGarch, 0).t();  
	  mnGarch = mnGarch - repmat(meanx.t(), t, 1);

	  //mnGarch.print("mnGarch = ");
	  //mean(mnGarch2, 0).t().print("mean(mnGarch2) = ");
	  
	  mat aC = (1.0/t) * mnGarch.t() * mnGarch;

	  double ds;
	  mat aS = Stats::cov2para(mnGarch, ds);
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
static DistMatrix<R> DistMatCov(DistMatrix<R> X) {

  int t = X.Height();
  int n = X.Width();

  double dt = (double) t;
  
  DistMatrix<R> colSum(n, 1, g);
  Zeros(n, 1, colSum);

  // compute column sums for the means
  for (int j = 0; j < n; j++)
	for (int i = 0; i < gRows; i++)
	  colSum.Set(j, 0, colSum.Get(j,0) + X.Get(i,j));

  // subtract column means
  for (int j = 0; j < n; j++)
	for (int i = 0; i < gRows; i++)
	  X.Set(i, j, X.Get(i,j) - colSum.Get(j,0)/dt);

  // cov = 1/t * X' * X'
  DistMatrix<R> sample(n, n, g);
  double alpha = 1.0/dt;
	
  Gemm(TRANSPOSE, NORMAL, alpha, X, X, (double)0, sample);

  return sample
}
*/
