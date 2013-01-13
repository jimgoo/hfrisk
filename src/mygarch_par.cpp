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

//-------------------------------------------------------------------------------
// <MAIN>

int main(int argc, char **argv) {

  cout << "\n\nMain for mygarch_par.cpp\n\n";

  // rank of this process
  int rank;

  // sizes of each group
  int size_world, size_garch, size_mat;

  // set group sizes <TODO>
  size_garch = 1;
  size_mat = 1;
  
  MPI_Group group_world, group_garch, group_mat;
  MPI_Comm comm_garch, comm_mat;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size_world);

  // check that group sizes sum to total world size
  assert(size_world == (size_garch + size_mat));

  // Set the inclusion rank arrays of form {minRank, maxRank, stride}.
  // The ordering is global = {garch, mat}.
  int range_garch[][3] = {0, size_garch-1, 1};
  int range_mat[][3] = {size_garch, size_garch+size_mat-1, 1};

  // set all groups
  MPI_Comm_group(MPI_COMM_WORLD, &group_world);
  MPI_Group_range_incl(group_world, 1, range_garch, &group_garch);
  MPI_Group_range_incl(group_world, 1, range_mat, &group_mat);
  
  // set communicators for each non-world group
  MPI_Comm_create(MPI_COMM_WORLD, group_garch, &comm_garch);
  MPI_Comm_create(MPI_COMM_WORLD, group_mat, &comm_mat);


  // <TODO>
  // Do all work here
  //	- Load HDF files into DistMat
  //	- Start backtest
  
  // free all communicators
  if (comm_garch != MPI_COMM_NULL) MPI_Comm_free(&comm_garch);
  if (comm_mat != MPI_COMM_NULL) MPI_Comm_free(&comm_mat);

  // free all groups
  MPI_Group_free(&group_mat);
  MPI_Group_free(&group_garch);
  MPI_Group_free(&group_world);

  // finalize
  MPI_Finalize();
  
  /*
  int myrank, ntasks;
  int size_work = 5;
  int size_res  = 3;

  mygarch_par p;

  if (myrank == 0) {
	p.master(size_work, size_res);
  } else {	  
  	p.slave(size_work, size_res);
  }

  if (comm_garch != MPI_COMM_NULL)
	MPI_Comm_free(&comm_garch);

  MPI_Group_free(&grp_garch);
  MPI_Group_free(&grp_world);
  
  MPI_Finalize();
  */
  return 0;

}
