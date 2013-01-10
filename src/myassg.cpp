
#include "mpi.h"
#include <armadillo>
#include "mystable.hpp"

#include <gsl/gsl_sf.h> // special function gamma for t-pdf
#include <gsl/gsl_poly.h> // root finding for arma0
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_sort.h>

// pars
#define SIZE_RET    20000  // 
#define SIZE_STKS   20     // number of stocks (n)
#define SIZE_PAIRS  210    // number of unique stock paris: n/2*(n+1)

// MPI
#define WORKTAG     1
#define DIETAG      2

int rows, cols;
int current_row = -1;

string dfile = "/Users/jimmiegoode/Documents/Glimm/Toolbox/ASSG/rand_ASSG_d=20.csv";
string saveAs = "/Users/jimmiegoode/Documents/Glimm/Toolbox/ASSG/mnSigmaHat.csv";

mat mnRet;

double vnRet[SIZE_RET];
double table[SIZE_PAIRS][2];

mat sigmaHat(SIZE_STKS, SIZE_STKS); //[SIZE_STKS][SIZE_STKS];
vec alphaHat(SIZE_STKS);

int verbose = 0;

int size_work = 2;
int size_res  = 4;

//-------------------------------------------------------------------------------
// <MASTER>

static int get_next_work_item(double* work) {

  current_row++;

  if (current_row >= SIZE_PAIRS) {
	return -1;
  }

  // fill the remaing work items with returns (oldest to newest)
  assert(size_work == 2);
  work[0] = table[current_row][0];
  work[1] = table[current_row][1];
  
  return 0;
}

//-------------------------------------------------------------------------------
// <MASTER>

static void onPairComplete(double* result) {
  
  // print out the results
  double pct =  ((double)current_row)/((double)SIZE_PAIRS)*100.00;

  int i = result[0];
  int j = result[1];
  double alpha = result[2];
  double sigma = result[3];
  
  sigmaHat(i,j) = sigma;

  if (i == j)
	alphaHat(i) = alpha;
  
  cout << setw(2) << pct << ": SigmaHat(" << i << "," << j << ") = " << sigma << endl;
  
  //if (i == j)
	//cout << "\talpha = " << alpha << endl;
}

static void onPairsComplete() {
  
  cout << "-----------------------------------------" << endl;
  cout << "copying lower diagonal...";
  for (int i = 1; i < cols; i++) {
	for (int j = 0; j < i; j++) {
	  sigmaHat(i,j) = sigmaHat(j,i);
	}
  }
  cout << "Done." << endl;
  
  double nAlphaHat = mean(alphaHat);
  sigmaHat.save(saveAs, csv_ascii);
  
  cout << "saved mnSigmaHat to: " << saveAs << endl;
  cout << "alphaHat = " << nAlphaHat << endl;
  cout << "-----------------------------------------" << endl;
}

//-------------------------------------------------------------------------------
// <MASTER>

static void master(void) {
  
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
  int ret = get_next_work_item(work);
	
  while (ret == 0) {

    // Receive results from a slave
    MPI_Recv(&result, size_res, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	onPairComplete(result);
	
    MPI_Send(&work, size_work, MPI_DOUBLE, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);

    // Get the next unit of work to be done
    ret = get_next_work_item(work);	
  }

  // There's no more work to be done, so receive all the outstanding results from the slaves. 
  for (rank = 1; rank < ntasks; ++rank) {
	
    MPI_Recv(&result, size_res, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	onPairComplete(result);
  }
  
  
  // Tell all the slaves to exit by sending an empty message with the DIETAG.
  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Send(0, 0, MPI_DOUBLE, rank, DIETAG, MPI_COMM_WORLD);
  }

  // Do the copula estimation now that we've got the marginals
  onPairsComplete();
  
}


//-------------------------------------------------------------------------------
// <SLAVE>

static void do_work(double* work, double* result) {
  int myrank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  int r = work[0];
  int c = work[1];

  assert(r < cols && c < cols);

  // copy indices 
  assert(size_res == 4);
  result[0] = work[0];
  result[1] = work[1];

  int iAlpha = 2;
  int iSigma = 3;
  
  double sum = 0.0;

  if (r == c) {

	vec X(rows);
	for (int i = 0; i < rows; i++) {
	  X(i) = vnRet[i + r*rows];
	}
	mystable::ts_struct fit = mystable::mle_nlopt(X, mystable::symAS);
	//mystable::ts_struct_print(fit);
	
	//cout << "rank " << myrank << " doing work: " << work[0] << "," << work[1] << "\n";

	double alpha = fit.pars[0];
	double sigma = fit.pars[1];
	
	result[iAlpha] = alpha;
	result[iSigma] = 2*sigma*sigma;

	/*
	  gsl_vector* p = gsl_vector_alloc(rows);

	  // sigma plus and minus 
	  for (int i = 0; i < rows; i++) {
	  gsl_vector_set(p, i, vnRet[i + r*rows]);
	  }
	
	  gsl_sort(p->data, p->stride, p->size);
	
	  double q1 = gsl_stats_quantile_from_sorted_data(p->data, p->stride, p->size, 0.72);
	  double q2 = gsl_stats_quantile_from_sorted_data(p->data, p->stride, p->size, 0.28);
	  double sigma = (q1 - q2)/1.654; // scale parameter
	  result[2] = 2*sigma*sigma; // 2*sigma^2 on the diagonal of Sigma matrix

	  gsl_vector_free(p);

	  if (verbose) {
	  cout << "rank " << myrank << " doing work: " << work[0] << "," << work[1]
	  << ", sigma = " << result[3] << "\n";
	  }
	*/
	
  } else {
	
	gsl_vector* pp = gsl_vector_alloc(rows);
	gsl_vector* pm = gsl_vector_alloc(rows);

	// sigma plus and minus 
	for (int i = 0; i < rows; i++) {
	  gsl_vector_set(pp, i, vnRet[i + r*rows] + vnRet[i + c*rows]);
	  gsl_vector_set(pm, i, vnRet[i + r*rows] - vnRet[i + c*rows]); 
	}
	
	gsl_sort(pp->data, pp->stride, pp->size);
	gsl_sort(pm->data, pm->stride, pm->size);
	
	double q1 = gsl_stats_quantile_from_sorted_data(pp->data, pp->stride, pp->size, 0.72);
	double q2 = gsl_stats_quantile_from_sorted_data(pp->data, pp->stride, pp->size, 0.28);
	double sigmaPlus = (q1 - q2)/1.654;

	q1 = gsl_stats_quantile_from_sorted_data(pm->data, pm->stride, pm->size, 0.72);
	q2 = gsl_stats_quantile_from_sorted_data(pm->data, pm->stride, pm->size, 0.28);
	double sigmaMinus = (q1 - q2)/1.654;

	double sigma = (sigmaPlus*sigmaPlus - sigmaMinus*sigmaMinus)/2.0;

	result[iSigma] = sigma;
	result[iAlpha] = -1.0;

	gsl_vector_free(pp);
	gsl_vector_free(pm);

	if (verbose) {
	  cout << "rank " << myrank << " doing work: " << work[0] << "," << work[1]
		   << ", q1 = " << q1 << ", q2 = " << q2
		   << ", sigmaPlus " << sigmaPlus << ", sigmaMinus " << sigmaMinus
		   << ", sigma = " << sigma << endl;
	}
	
  }
  
}


//-------------------------------------------------------------------------------
// <SLAVE>

static void slave(void) {

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

	do_work(work, result);
	
    // Send the result back 
    MPI_Send(&result, size_res, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
}

//-------------------------------------------------------------------------------
// <MASTER>

static bool distribute(int rank) {

  if (rank == 0) {
	rows = mnRet.n_rows;
	cols = mnRet.n_cols;
  }
  
  MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size_work, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size_res , 1, MPI_INT, 0, MPI_COMM_WORLD);  
  MPI_Bcast(&vnRet, SIZE_RET, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  //cout << ">> rank = " << rank << ", rows = " << rows << ", cols = " << cols << ", " << vnRet[3] << endl;

}

//-------------------------------------------------------------------------------

int main(int argc, char **argv) {
  
  int myrank, ntasks;
  double t_all;

  // initialize MPI
  MPI_Init(&argc, &argv);
  
  // find out my identity in the default communicator 
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);

  if (myrank == 0) {
	
	mnRet.load(dfile, csv_ascii);
	cout << "Loaded " << mnRet.n_rows << " by " << mnRet.n_cols << " matrix.\n";
	
	assert(mnRet.n_rows*mnRet.n_cols == SIZE_RET);
	
	int idx = 0;
	for (int c = 0; c < mnRet.n_cols; c++) {
	  for (int r = 0; r < mnRet.n_rows; r++) {
		vnRet[idx] = mnRet(r,c);
		idx++;
	  }
	}

	// table with all index combinations
	int n = mnRet.n_cols;
	assert(SIZE_PAIRS == n/2*(n+1));
	
	int c = 0;
	for (int i = 0; i < SIZE_STKS; i++) {
	  for (int j = i; j < SIZE_STKS; j++) {
		table[c][0] = i;
		table[c][1] = j;
		c++;
	  }
	}
  }
  
  int hasData = distribute(myrank);
  
  if (myrank == 0) {
  	master();
  } else {	  
  	slave();
  }

  MPI_Finalize();
  return 0;
}
