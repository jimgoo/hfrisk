/*
  Author: Jimmie Goode
  Created: 2013-01-10
*/


#include "mygarch_par.hpp"


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

// INPUT:
// 	vector<string> files: list of HDF5 paths
//	int iBeg, iEnd: indices to start and end with

// OUTPUT:
//  DistMat mnResults: matrix of parameters and residuals

//-------------------------------------------------------------------------------
// <MAIN>

int main(int argc, char **argv) {

  cout << "\n\nMain for mygarch_par.cpp\n\n";

  int myrank, ntasks;

  int size_work = 5;
  int size_res  = 3;

  mygarch_par p;

  Initialize( argc, argv );
  mpi::Comm comm = mpi::COMM_WORLD;
  myrank = mpi::CommRank(comm);
  ntasks = mpi::CommSize(comm);
  
  if (myrank == 0) {
	p.master(size_work, size_res);
  } else {	  
  	p.slave(size_work, size_res);
  }

  Finalize();

  return 0;

}
