/*
  Author: Jimmie Goode
  Created: 2013-01-10
*/


#include "mygarch_par.hpp"


//-------------------------------------------------------------------------------
// <MASTER>

mygarch_par::mygarch_par(MPI_Comm comm_,
							  DistMatrix<double,STAR,VC> returns_,
							  DistMatrix<double,STAR,VC> &results_) {
  this->comm = comm_;
  this->returns = returns_;
  this->results = results_;
}

//-------------------------------------------------------------------------------
// <MASTER>

bool mygarch_par::get_next_work_item(double *work) {

  cout << "GARCH: get_next_work\n";

  if (iters > 5) {
	cout << "GARCH: max iters, breaking\n";
	return true;
  }
  
  return true;
}

//-------------------------------------------------------------------------------
// <SLAVE>

void mygarch_par::do_work(double *work, double *result, int size_work, int size_res) {

  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  
  cout << "GARCH: doing work, rank = " << myrank
	   << ", size_work = " << size_work << ", size_res = " << size_res << "...\n";
  
}

//-------------------------------------------------------------------------------
// <MASTER>

void mygarch_par::onWorkComplete(double *result) {
  cout << "GARCH: work_done\n";
}

//-------------------------------------------------------------------------------
// <MASTER>

void mygarch_par::onAllWorkComplete() {
  cout << "GARCH: all work done\n";
}


//-------------------------------------------------------------------------------
// <MASTER>
//
// DESC:
//	Fit all garch parameters for a given period
//
// INPUT:
// 	vector<string> files - list of HDF5 paths for each asset
//	int iBeg, iEnd - indices to start and end with
//  elem::Grid - Elemental process grid to construct output matrix

// OUTPUT:
//  DistMat mnResults: matrix of parameters and residuals

// void mygarch_par::fit(DistMatrix<double> returns, DistMatrix<double> &results) {

//   // initialilze the counter
//   currentAsset = 0;  
  
// }
