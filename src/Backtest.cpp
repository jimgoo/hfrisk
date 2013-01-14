/*
  Author: Jimmie Goode
  Created: 2013-01-10
*/


#include "Backtest.hpp"

//-------------------------------------------------------------------------------

void checkErr(int i) {
  if (i != MPI_SUCCESS) {
	assert(0 && "MPI Error");
  }
}

//-------------------------------------------------------------------------------

void printRange(int r[][3], string name) {
  int i = 0;
  cout << name << " = [" << r[i][0] << ", " << r[i][1] << ", " << r[i][2] << "]" << endl;
}

//-------------------------------------------------------------------------------

int main(int argc, char **argv) {

  // rank of this process
  int rank, ie;

  // sizes of each group
  int size_world, size_garch, size_mat;

  // set group sizes <TODO>
  //size_garch = 6;
  //size_mat = 1;
  
  MPI_Group group_world, group_garch, group_mat;
  MPI_Comm comm_garch, comm_mat;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size_world);

  size_garch = size_world - 1;
  size_mat = 1;
  
  // check that group sizes sum to total world size
  assert(size_world == (size_garch + size_mat));

  // Set the inclusion rank arrays of form {minRank, maxRank, stride}.
  // The ordering is global = {garch, mat}.
  int range_garch[][3] = {0, size_garch-1, 1};
  int range_mat[][3] = {size_garch, size_garch+size_mat-1, 1};

  if (rank == 0) {
	printRange(range_garch, "range_garch");
	printRange(range_mat, "range_mat");
  }

  // set all groups
  MPI_Comm_group(MPI_COMM_WORLD, &group_world);
  MPI_Group_range_incl(group_world, 1, range_garch, &group_garch);
  MPI_Group_range_incl(group_world, 1, range_mat, &group_mat);
  
  // set communicators for each non-world group
  MPI_Comm_create(MPI_COMM_WORLD, group_garch, &comm_garch);
  MPI_Comm_create(MPI_COMM_WORLD, group_mat, &comm_mat);


  // print info for whichever group this process belongs to
  int grank, gsize;

  if (rank < size_garch) {

	MPI_Comm_rank(comm_garch, &grank);
	MPI_Comm_size(comm_garch, &gsize);
	cout << "GARCH Member: ranks = [" << rank << ", " << grank << "], groupSize = " << gsize << endl;
	  
  } else if (rank >= size_garch) {
	  
	MPI_Comm_rank(comm_mat, &grank);
	MPI_Comm_size(comm_mat, &gsize);
	cout << "MAT Member: ranks = [" << rank << ", " << grank << "], groupSize = " << gsize << endl;

	string fname = "/Users/jimmiegoode/Documents/Glimm/github/hfrisk/data/csi/patentData/csi_20030101_20120801_v3/logret.abin";

	file_type type = arma_binary;
	  
	Grid grid( comm_mat );

	//DistMatrix<double,STAR,VC> returns = IO::arma2distMat(fname, arma::arma_binary, grid);
	mat armaRet;
	armaRet.load(fname, type);
	armaRet = armaRet(span(0,5), span(0,2));

	armaRet.print("armaRet = ");
	  
	const int r = armaRet.n_rows;
	const int c = armaRet.n_cols;

	DistMatrix<double,STAR,VC> dmRet(r, c, grid);
	DistMatrix<double,STAR,VC> dmRes(r, c, grid);

	for (int i = 0; i < r; i++)
	  for (int j = 0; j < c; j++)
		dmRet.Set(i, j, armaRet(i,j));

	dmRet.Print("dmRet = ");

	mygarch_par gp(comm_garch, dmRet, &dmRes);
	
  }
  

  MPI_Barrier(MPI_COMM_WORLD);

  if (rank < size_garch) {
	MPI_Comm_rank(comm_garch, &grank);
	MPI_Comm_size(comm_garch, &gsize);

	// mygarch_par gp(comm_garch, dmRet, &dmRes);

	// if (grank == 0) {
	//   p.master();
	// } else {
	//   p.slave();
	// }
  }
 
  // free all communicators
  if (comm_garch != MPI_COMM_NULL) MPI_Comm_free(&comm_garch);
  if (comm_mat != MPI_COMM_NULL) MPI_Comm_free(&comm_mat);

  // free all groups
  MPI_Group_free(&group_mat);
  MPI_Group_free(&group_garch);
  MPI_Group_free(&group_world);

  // finalize
  MPI_Finalize();
  
  // int myrank, ntasks;
  // int size_work = 5;
  // int size_res  = 3;

  // mygarch_par p;

  // if (myrank == 0) {
  // 	p.master(size_work, size_res);
  // } else {	  
  // 	p.slave(size_work, size_res);
  // }

  // if (comm_garch != MPI_COMM_NULL)
  // 	MPI_Comm_free(&comm_garch);

  // MPI_Group_free(&grp_garch);
  // MPI_Group_free(&grp_world);
  
  // MPI_Finalize();
  
  return 0;
}


//-------------------------------------------------------------------------------

int main4(int argc, char **argv) {

  int rank[2], size[2], namelen;
  int xranks[] = { 0 };
  int iranks[] = { 6 };
  int rranks[][3] = {6,6,1};

  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Group group_world, group_garch;

  MPI_Comm comm_garch;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank[0]);
  MPI_Comm_size(MPI_COMM_WORLD, &size[0]);

  MPI_Get_processor_name(processor_name, &namelen);

  MPI_Comm_group(MPI_COMM_WORLD, &group_world);

  //MPI_Group_excl(group_world, 1, xranks, &group_garch);
  //MPI_Group_incl(group_world, 1, iranks, &group_garch);
  MPI_Group_range_incl(group_world, 1, rranks, &group_garch);

  MPI_Comm_create(MPI_COMM_WORLD, group_garch, &comm_garch);

  printf("Hello world! I'm rank %d of %d on %s\n", rank[0], size[0], processor_name);

  // If this rank is in the subgroup, get it's subgroup rank.
  if (rank[0] == 6) {
	MPI_Comm_rank(comm_garch, &rank[1]);
	MPI_Comm_size(comm_garch, &size[1]);
	cout << "rank[0] = " << rank[0] << ", rank[1] = " << rank[1] << endl;
	cout << "size[0] = " << size[0] << ", size[1] = " << size[1] << endl;
  }

 if (comm_garch != MPI_COMM_NULL)
	MPI_Comm_free(&comm_garch);

  MPI_Group_free(&group_garch);

  MPI_Group_free(&group_world);

  MPI_Finalize();

  return 0;
}
