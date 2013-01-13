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

int main(int argc, char **argv) {

  int rank[2], size[2], namelen, xranks[] = { 0 };
  int iranks[] = { 1 };

  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Group group_world, group_garch;

  MPI_Comm comm_garch;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank[0]);
  MPI_Comm_size(MPI_COMM_WORLD, &size[0]);

  MPI_Get_processor_name(processor_name, &namelen);

  MPI_Comm_group(MPI_COMM_WORLD, &group_world);

  //MPI_Group_excl(group_world, 1, xranks, &group_garch);
  MPI_Group_incl(group_world, 1, iranks, &group_garch);

  MPI_Comm_create(MPI_COMM_WORLD, group_garch, &comm_garch);

  printf("Hello world! I’m rank %d of %d on %s\n", rank[0], size[0], processor_name);

  if (rank[0]) {
	int rank1;
	MPI_Comm_rank(comm_garch, &rank1);
	cout << "rank1 = " << rank1 << endl;
  }

 if (comm_garch != MPI_COMM_NULL)
	MPI_Comm_free(&comm_garch);

  MPI_Group_free(&group_garch);

  MPI_Group_free(&group_world);

  MPI_Finalize();

  return 0;
}

//-------------------------------------------------------------------------------

int main3(int argc, char **argv) {

  // rank of this process
  int rank[2], ie;

  // sizes of each group
  int size_world, size_garch, size_mat;

  // set group sizes <TODO>
  size_garch = 1;
  size_mat = 1;
  
  MPI_Group group_world, group_garch, group_mat;
  MPI_Comm comm_garch, comm_mat;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank[0]);
  MPI_Comm_size(MPI_COMM_WORLD, &size_world);

  // check that group sizes sum to total world size
  assert(size_world == (size_garch + size_mat));

  // Set the inclusion rank arrays of form {minRank, maxRank, stride}.
  // The ordering is global = {garch, mat}.
  int range_garch[][3] = {0, size_garch-1, 1};
  ///int range_mat[][3] = {size_garch, size_garch+size_mat-1, 1};

  // set all groups
  MPI_Comm_group(MPI_COMM_WORLD, &group_world);

  if (range_garch[0][0] == range_garch[0][1]) {
	cout << "One processor for garch\n";
	
	int rnk[] = {1};
	//MPI_Group_incl(group_world, 2, rnk, &group_garch);
	MPI_Group_excl(group_world, 1, rnk, &group_garch);
	
  } else {

	//ie = MPI_Group_range_incl(group_world, 1, range_garch, &group_garch);
	assert(0 && "---> FAIL");
  }
  
  //// ie = MPI_Group_range_incl(group_world, 1, range_mat, &group_mat);
  
  // set communicators for each non-world group
  ie = MPI_Comm_create(MPI_COMM_WORLD, group_garch, &comm_garch);
  //// MPI_Comm_create(MPI_COMM_WORLD, group_mat, &comm_mat);

  // <TODO>
  // - Load HDF files into DistMat
  // - Start backtest

  if (rank[0]) {
	int rank1, rank2;
	//MPI_Comm_rank(comm_garch, &rank1);
	//MPI_Comm_rank(comm_garch, &rank2);
	cout << "rank1 = " << rank1 << ", rank2 = " << rank2 << endl;
  }
  
  // free all communicators
  if (comm_garch != MPI_COMM_NULL) MPI_Comm_free(&comm_garch);
  //// if (comm_mat != MPI_COMM_NULL) MPI_Comm_free(&comm_mat);

  // free all groups
  //// MPI_Group_free(&group_mat);
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



int main2(int argc, char **argv)
{

  int rank[2], size[2], namelen, xranks[] = { 0 };

  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Group mpi_group_world, group_slaves;

  MPI_Comm comm_slaves;

  int send_val, recv_val, send_val2, recv_val2;


  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank[0]);

  MPI_Comm_size(MPI_COMM_WORLD, &size[0]);

  MPI_Get_processor_name(processor_name, &namelen);


  MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_world);

  MPI_Group_excl(mpi_group_world, 1, xranks, &group_slaves);

  MPI_Comm_create(MPI_COMM_WORLD, group_slaves, &comm_slaves);

  printf("Hello world! I’m rank %d of %d on %s\n", rank[0], size[0], processor_name);

  if (rank[0]) {
	int rank1;
	MPI_Comm_rank(comm_slaves, &rank1);
	cout << "rank1 = " << rank1 << endl;
  }
  
  /*
  if (rank[0]) {

    MPI_Comm_rank(comm_slaves, &rank[1]);
    MPI_Comm_size(comm_slaves, &size[1]);

    printf("In the slave universe I’m rank %d of %d on %s\n", rank[1], size[1], processor_name);

    send_val = size[1];

    MPI_Reduce(&send_val, &recv_val, 1, MPI_INT, MPI_SUM, 0, comm_slaves);

    if (!rank[1])
	  printf("Slave leader received reduced value %d\n", recv_val);

  }

  send_val2 = size[0];

  MPI_Reduce(&send_val2, &recv_val2, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  if (!rank[0])
	printf("Master received reduced value %d\n", recv_val2);
  */
  
  if (comm_slaves != MPI_COMM_NULL)
	MPI_Comm_free(&comm_slaves);

  MPI_Group_free(&group_slaves);

  MPI_Group_free(&mpi_group_world);

  MPI_Finalize();

}
