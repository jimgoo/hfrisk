/*
  Author: Jimmie Goode
  Created: 2013-01-10
*/


#include "MasterSlave.hpp"

// MPI
#define WORKTAG  1
#define DIETAG   2


//-------------------------------------------------------------------------------
// <MASTER>

bool MasterSlave::get_next_work_item(double* work) {
  return false;
}

//-------------------------------------------------------------------------------
// <SLAVE>

void MasterSlave::do_work(double *work, double *result, int size_work, int size_res) {
  
}

//-------------------------------------------------------------------------------
// <MASTER>

void MasterSlave::onWorkComplete(double *result) {
  
}

//-------------------------------------------------------------------------------
// <MASTER>

void MasterSlave:: onAllWorkComplete() {

}

//-------------------------------------------------------------------------------
// <MASTER>

void MasterSlave::master(int size_work, int size_res) {

  iters = 0;
  
  double work[size_work];  
  double result[size_res];
  
  int myrank, ntasks, rank;
  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // Initialize the processor map
  int procMap[ntasks];

  for (int i = 0; i < ntasks; i++) {
	procMap[i] = -1;
  }
  
  // Seed the slaves; send one unit of work to each slave. 
  for (rank = 1; rank < ntasks; ++rank) {
	
    // Find the next item of work to do
	get_next_work_item(work);

    MPI_Send(&work,                 // message buffer 
             size_work,             // size of data item 
             MPI_DOUBLE,            // data item type is double 
             rank,                  // destination process rank 
             WORKTAG,               // user chosen message tag 
             MPI_COMM_WORLD);       // default communicator 	
  }

  // Loop over getting new work requests until there is no more work to be done
  bool hasData = get_next_work_item(work);
	
  while (hasData) {

    // Receive results from a slave
    MPI_Recv(&result, size_res, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
	onWorkComplete(result);
	
    MPI_Send(&work, size_work, MPI_DOUBLE, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);

    // Get the next unit of work to be done
    hasData = get_next_work_item(work);	
  }

  // There's no more work to be done, so receive all the outstanding results from the slaves. 
  for (rank = 1; rank < ntasks; ++rank) {
	
    MPI_Recv(&result, size_res, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	onWorkComplete(result);
  }
  
  // Tell all the slaves to exit by sending an empty message with the DIETAG.
  for (rank = 1; rank < ntasks; ++rank) {
	
    MPI_Send(0, 0, MPI_DOUBLE, rank, DIETAG, MPI_COMM_WORLD);
  }

  // Do the copula estimation now that we've got the marginals
  onAllWorkComplete();
  
}


//-------------------------------------------------------------------------------
// <SLAVE>

void MasterSlave::slave(int size_work, int size_res) {

  iters++;
  
  double work[size_work];
  double result[size_res];
  
  MPI_Status status;

  while (1) {

    // Receive a message from the master 
    MPI_Recv(&work, size_work, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    // Check the tag of the received message.
    if (status.MPI_TAG == DIETAG) {  
      return;
    }

	do_work(work, result, size_work, size_res);
	
    // Send the result back 
    MPI_Send(&result, size_res, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
  
}

//-------------------------------------------------------------------------------
// <MASTER>

void MasterSlave::distribute(int rank) {
  
  // MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // MPI_Bcast(&size_work, 1, MPI_INT, 0, MPI_COMM_WORLD);
  // MPI_Bcast(&size_res , 1, MPI_INT, 0, MPI_COMM_WORLD);  
  // MPI_Bcast(&vnRet, SIZE_RET, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}


//-------------------------------------------------------------------------------
//

/*

//-------------------------------------------------------------------------------
// <MAIN>

int main(int argc, char **argv) {

  cout << "\n\nMain for MasterSlave.cpp\n\n";

  int myrank, ntasks;

  int size_work = 5;
  int size_res  = 3;

  MasterSlave ms;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  ms.distribute(myrank);
  
  if (myrank == 0) {
  	ms.master(size_work, size_res);
  } else {	  
  	ms.slave(size_work, size_res);
  }

  MPI_Finalize();
  
  return 0;

}

*/
