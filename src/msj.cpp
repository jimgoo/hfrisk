//
// Author: Jimmie Goode
// Created: 2012-09-01
//

//-------------------------------------------------------------------------------

#include "msj.h"

//-------------------------------------------------------------------------------
// MPI tag type constants

#define WORKTAG     1
#define DIETAG      2

//-------------------------------------------------------------------------------
// Elemental types

// Typedef our real and complex types to 'R' and 'C' for convenience
// typedef double ELEM_R;
// typedef Complex<R> ELEM_C;

//DistMatrix<double> mnResultsE;

//-------------------------------------------------------------------------------

string LINE = "=================================================================\n";

// current column in marginal estimation
int currentColumn = -1;

// sizes of things
int rows, cols, size_res, size_pars, size_work;

// model orders ARMA(R,M)-GARCH(P,Q)
int R = 1;
int M = 1;
int P = 1;
int Q = 1;
int maxRMPQ;

int iS;// = 5109;  // number of stocks
int iL = 1008;  // lookback window size
int t =  iL-1;  // idx of starting time (iL-1 <=> first possible day)
int iT;         // total number of periods
int iBeg, iEnd; // indices of the dates files to backtest over

mat mnRet;     // returns for current period
mat mnRetAll;  // returns for all periods

mat mnStart;   // hot-start parameter matrix
mat mnResults; // matrix of marginal estimation results
vec vnDates;   // dates corresponding to the rows of mnRetAll

FILE* fReport;       // report file
double MPI_Wtime();  // MPI format timestamp

double t_period = 0.; // time for all work on a given period
double t_garch  = 0.; // time for parallel GARCH estimations
// double t_adjust = 0.;
// double t_sample = 0.;

int margOnly = 0; // only fit the marginals
int doShrink = 1;
double gammaScale = 0.01;
int doCheckEigs = 1;

int depStruct = 1; // 1 = skewt, 2 = ASSG 
int innovType = 1; // 1 = skewt, 2 = stdAS, 3 = CTS, 4 = NTS

int marginalCount  = -1; // number of marginals to estimate (-1 for all)

string dataFile  = "";                          // data file of returns
string reportFile = "../data/exports/mnResults"; // output file for report

EST_METHOD estMethod = NMSIMPLEX;             // univariate skewed t MLE method
START_TYPE startType = COLD;                  // COLD or HOT
string startFile = "../data/mnStart.csv";     // file with starting copula parameters
int verbose = 0;                              // print out info
int iMaxT = -1;                               // maximum period index to run to

// Monte-Carlo parameters
int nSim = 10000;           // number of copula samples for MC VaR
double VaR_epsilon = 0.01;  // 100*(1-d)% confidence level for VaR and CVaR1

unsigned long int seed = 1; // RNG seed
std::ostringstream ss;

double cdfTol = 1.0e-6;
double beginDate = 20070104.;
double endDate = 20080101.;

int VaR_total = 0;
int CVaR_total = 0;

int isSimple = 1;

myskewt_lut lut;
string lut_path = "";

// file_type lut_type = hdf5_binary;
// string lut_ext = ".h5";
file_type lut_type = csv_ascii;
string lut_ext = ".csv";

vec wts;
int doLUT = 0;
int numCores;
int goBig = 0; 

//-------------------------------------------------------------------------------

// // sizes of each group
// int size_world, size_garch, size_mat;

// MPI_Group group_world, group_garch, group_mat;
// MPI_Comm comm_garch, comm_mat;


//-------------------------------------------------------------------------------
// get memory usage

static void process_mem_usage() {
  printf(LINE.c_str());
  system("free -m");
  printf(LINE.c_str());
}

//-------------------------------------------------------------------------------
// <main>

int main(int argc, char **argv) {
  
  //==========================================
  
  int rank;
  double t_all;

  //==========================================
  // initialize MPI (old way)
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numCores);

  //==========================================
 
  // Initialize( argc, argv );
  // mpi::Comm comm = mpi::COMM_WORLD;
  // myrank = mpi::CommRank(comm);
  // ntasks = mpi::CommSize(comm);

  // cout << "ntasks = " << ntasks << endl;
  // Grid g( comm );
  
  //==========================================
  // MPI with different groups
  /*
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
  */
  
  
  //==========================================

  // parse command line args
  if (rank == 0) {
	
	// print help and exit for "-h"
	if (argc == 2 && strncmp(argv[1], "-h", 2) == 0)  {
	  cout << "See msj.cpp for arguments." << endl;

	  MPI_Finalize();
	  
	  return 0;
	}

	// print help and exit for "-h"
	if (argc == 3 && strncmp(argv[1], "-test", 10) == 0)  {
	  
	  //==========================================
	  // Testing of individual classes only.
  
	  int testCase = atoi(argv[2]);
	  mystable::tester(testCase);
	  return 0;

	  MPI_Finalize();
	  
	  return 0;
	}

	if (argc == 2 && strncmp(argv[1], "-t", 2) == 0) {
	  
	  cout << "Building table...\n";

	  double t_lut = MPI::Wtime();

	  vec beta1 = linspace(-1,1,20);
	  vec beta2 = linspace(-1,1,20);
	  double df = 5.0;
      double d = 30;
	  int nSim = 10000;
	  
	  vec xGrid = linspace(-20,20,300);
	  
	  lut.buildTable(beta1, beta2, df, d, nSim, xGrid);

	  cout << "Elapsed time for table generation: " << MPI::Wtime()-t_lut << endl;

	  lut.save("../../../LUTs/test2", lut_type, lut_ext);
	  
	  MPI_Finalize();
	  return 0;
	}
	
	for (int i = 1; i < argc; i++) {
	  if (i + 1 != argc) {
		if (strncmp(argv[i],"-rf",3) == 0) {
		  reportFile = boost::lexical_cast<string>(argv[i+1]);	    
		} else if  (strncmp(argv[i],"-mc",3) == 0)  {
		  marginalCount = boost::lexical_cast<int>(argv[i+1]);
		} else if  (strncmp(argv[i],"-df",3) == 0)  {
		  dataFile = boost::lexical_cast<string>(argv[i+1]);
		} else if (strncmp(argv[i],"-st",3) == 0) {
		  startType = static_cast<START_TYPE> (boost::lexical_cast<int>(argv[i+1]));
		} else if (strncmp(argv[i],"-sf",3) == 0) {
		  startFile = boost::lexical_cast<string>(argv[i+1]);		  
		} else if (strncmp(argv[i],"-v",2) == 0) {
		  verbose = boost::lexical_cast<int>(argv[i+1]);
		} else if (strncmp(argv[i],"-estMethod",11) == 0) {
		  estMethod = static_cast<EST_METHOD> (boost::lexical_cast<int>(argv[i+1]));
		} else if (strncmp(argv[i],"-maxT",6) == 0) {
		  iMaxT = boost::lexical_cast<int>(argv[i+1]);
		} else if (strncmp(argv[i],"-nSim",7) == 0) {
		  nSim = boost::lexical_cast<int>(argv[i+1]);
		} else if (strncmp(argv[i],"-innovType",12) == 0) {
		  innovType = boost::lexical_cast<int>(argv[i+1]);
		} else if (strncmp(argv[i],"-margOnly",11) == 0) {
		  margOnly = boost::lexical_cast<int>(argv[i+1]);
		} else if (strncmp(argv[i],"-depStruct",13) == 0) {
		  depStruct = boost::lexical_cast<int>(argv[i+1]);
		} else if (strncmp(argv[i],"-doCheckEigs",25) == 0) {
		  doCheckEigs = boost::lexical_cast<int>(argv[i+1]);
		} else if (strncmp(argv[i],"-beginDate",15) == 0) {
		  beginDate = boost::lexical_cast<double>(argv[i+1]);
		} else if (strncmp(argv[i],"-endDate",15) == 0) {
		  endDate = boost::lexical_cast<double>(argv[i+1]);
		} else if (strncmp(argv[i],"-gammaScale",15) == 0) {
		  gammaScale = boost::lexical_cast<double>(argv[i+1]);
		} else if (strncmp(argv[i],"-lut_path", 12) == 0) {
		  lut_path = boost::lexical_cast<string>(argv[i+1]);
		} else if (strncmp(argv[i],"-doLUT",7) == 0) {
		  doLUT = boost::lexical_cast<int>(argv[i+1]);
		} else if (strncmp(argv[i],"-goBig",10) == 0) {
		  goBig = boost::lexical_cast<int>(argv[i+1]);
		}
	  }
	}
	
	// Print out runtime parameters.
	//<TODO> Check that parameters are actually valid.
	  
	ss << LINE
	   << "Parameters" << endl
	   << LINE
	   << "reportFile = "   << reportFile  << "\nmarginalCount = " << marginalCount
	   << "\ndataFile = "   << dataFile    << "\nstartType = "     << startType
	   << "\nstartFile = "  << startFile   << "\nverbose = "       << verbose
	   << "\nestMethod = "  << estMethod   << "\niMaxT = "         << iMaxT
	   << "\nnSim = "       << nSim        << "\ninnovType = "     << innovType
	   << "\nmargOnly = "   << margOnly    << "\ndepStruct = "     << depStruct
	   << "\ndoChkEigs = "  << doCheckEigs << "\nbeginDate = "     << (int)beginDate
	   << "\nendDate = "    << (int)endDate<< "\ngammaScale = "    << gammaScale
	   << "\nlut_path = "   << lut_path    << "\ndoLUT = "         << doLUT
	   << "\nnumCores = "   << numCores    << "\ngoBig = "         << goBig
	   << endl
	   << LINE;
	
	cout << ss.str();
	
    // Set random number generator seed
	Stats::gsl_rng_init(seed);
	
	// set the number of elements in the MPI marginal message 
	vec orders;
	orders << R << endr << M << endr << P << endr << Q;
	maxRMPQ = max(orders);

	// the size of the main MPI ARMA-GARCH message
	size_res = mygarch::iNonCoeffs + (1 + R + M) + (1 + P + Q) + 1 + 3*maxRMPQ + iL + 4;

	// import all returns
	cout << LINE;
	cout << "----> Loading return data..." << endl;
	
	//mnRetAll.load(dataFile + "_logret.csv", csv_ascii);
	//mnRetAll.load(dataFile + "/logret.abin", arma_binary);
	mnRetAll.load(dataFile + "/logret.h5", hdf5_binary);

	// replicate data for larger dimension testing
	if (goBig > mnRetAll.n_cols) {
	  int iReps = (int) goBig/mnRetAll.n_cols;
	  cout << "----> iReps = " << iReps << endl;
	  mnRetAll = repmat(mnRetAll, 1, iReps + 1);
	  mnRetAll = mnRetAll(span::all, span(0, goBig-1));
	}
	
	vnDates.load(dataFile + "/dates.csv");
	vnDates = vnDates(span(1,vnDates.n_rows-1));

	uvec idx0 = find(vnDates >= beginDate, 1, "first");
	uvec idx1 = find(vnDates <= endDate, 1, "last");

	cout << "----> Border dates: " << (int)vnDates(idx0(0)) << ", " << (int)vnDates(idx1(0)) << endl;

	// get column indices
	int iTkr;
	if (marginalCount > 0) {
	  iTkr = marginalCount-1;
	} else {
	  iTkr = mnRetAll.n_cols-1;
	}

	iBeg = idx0(0) - iL + 1;
	iEnd = idx1(0);

	vnDates = vnDates(span(iBeg,iEnd));
	mnRetAll = mnRetAll(span(iBeg,iEnd), span(0,iTkr));
	
	assert(vnDates.n_rows == mnRetAll.n_rows);
	
	cout << "size(mnRetAll) = " << mnRetAll.n_rows << ", " << mnRetAll.n_cols << endl;
	cout << "iBeg = " << iBeg << endl;
	cout << "iEnd = " << iEnd << endl;

	//iT = mnRetAll.n_rows;
	iT = vnDates.n_rows;
	iS = mnRetAll.n_cols;

	// set equal portfolio weights
	wts = 1.0/((double) iS) * ones((double) iS);

	cout << LINE;

	// open the VaR report file
	system(("mkdir " + reportFile).c_str());
	system(("mkdir " + reportFile + "/mnResults").c_str());
	fReport = fopen((reportFile + "/forc.csv").c_str(), "w+");
	fclose(fReport);

	// write the configuration to file
	FILE* fTime;
	fTime = fopen((reportFile + "/conf.txt").c_str(), "w+");	
	fprintf(fTime, "%s\n", (ss.str()).c_str());
	fclose(fTime);

	
	// load the lookup table
	if (doLUT) {
	  lut.load(lut_path, lut_type, lut_ext);
	}
	process_mem_usage();
	
	// start the timer for the whole shooting match
	t_all = MPI::Wtime();

  }

  // send time indices to all
  MPI_Bcast(&iT, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);


  while (1) {

  	int hasData = distribute(rank);

  	if (hasData == false) {

	  if (verbose) cout << "----> hasData = false\n";
	  
  	  if (rank == 0) {

		//if (verbose) cout << "----> closing report file\n";
		// close report file
		//fclose(fReport);
		
  		// export total runtime
  		t_all = MPI::Wtime() - t_all;
  	  }

	  if (verbose) cout << "----> finalizing MPI\n";
	  MPI_Finalize();
	  
  	  return 0;
  	}
	
  	if (rank == 0) {
  	  master();
  	} else {	  
  	  slave();
  	}
  }
}
 


//------------------------------------------------------------------------------- 
// update return matrix and distribute estimation parameters

static bool distribute(int rank) {

  // Return false if we've reached the end of the data
  if (t+1 >= iT)
	return false;
  
  if (rank == 0) {

	t_period = MPI::Wtime();
	t_garch = MPI::Wtime();
	
	// reset the current column
	currentColumn = -1;

    // set the lookback returns
	mnRet = mnRetAll(span(t-iL+1, t), span::all);
	t++;
	
	// set sizes of mnRet
	// <TODO><VIP> This will fail when the number of stocks changes from day to day
	rows = mnRet.n_rows; // iL
	cols = mnRet.n_cols; // iS

	// set the number of hot start parameters
	size_pars = 0;
	
	// Total length of work with hot-start parameters included
	size_work = rows + size_pars;

	// reset copula parameters
	mnResults = zeros(size_res, cols);

  }
  
  // Broadcast parameters to all ranks
  MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size_work, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size_res , 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&size_pars, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&startType, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&estMethod, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&depStruct, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast(&innovType, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&isSimple, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  return true;
	
}



//-------------------------------------------------------------------------------

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

	procMap[rank] = currentColumn; // update the mnRet index in the processor map
  }

  // Loop over getting new work requests until there is no more work to be done
  int ret = get_next_work_item(work);
	
  while (ret == 0) {

    // Receive results from a slave
    MPI_Recv(&result, size_res, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	//printf("received result> from rank %u, %f, %f\n", status.MPI_SOURCE, result[0], result[3]);
    //printProcMap(procMap, ntasks);
	
	onMarginalComplete(procMap[status.MPI_SOURCE], result);
	
    // Send the slave a new work unit
  	//printf("--> sending work 2> %f, %f\n", work[0], work[9]);

    MPI_Send(&work, size_work, MPI_DOUBLE, status.MPI_SOURCE, WORKTAG, MPI_COMM_WORLD);

	// update the mnRet index in the processor map
    procMap[status.MPI_SOURCE] = currentColumn;
	
    // Get the next unit of work to be done
    ret = get_next_work_item(work);	
  }

  // There's no more work to be done, so receive all the outstanding results from the slaves. 

  for (rank = 1; rank < ntasks; ++rank) {
	
    MPI_Recv(&result, size_res, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	
	onMarginalComplete(procMap[status.MPI_SOURCE], result);

	/*
	  if (rank == ntasks - 1) {
	  double vm, rss;
	  process_mem_usage(vm, rss);
	  printf("Final memory at largest rank node>> VM: %.2f, RSS: %.2f\n", vm, rss);
	  }
	*/
  }
  
  //printf("Sending DIETAGs to all\n");
  
  // Tell all the slaves to exit by sending an empty message with the DIETAG.
  for (rank = 1; rank < ntasks; ++rank) {
    MPI_Send(0, 0, MPI_DOUBLE, rank, DIETAG, MPI_COMM_WORLD);
  }

  // Do the copula estimation now that we've got the marginals
  onMarginalsComplete();
  
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

	if (isSimple) {
	  do_work_simple(work, result);
	} else {
	  do_work_copula(work, result);
	}
	
    // Send the result back 
    MPI_Send(&result, size_res, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

  }
}

//-------------------------------------------------------------------------------
// <MASTER>

static int get_next_work_item(double* work) {
  
  // increment the current column number
  currentColumn++;

  if (currentColumn >= cols) {
  	return -1;
  }

  // Import returns from HDF5 file
  /*
  string h5dir = dataFile + "/hdf5";
  vec ret;
  ret.load(h5dir + "/EQ_" + boost::lexical_cast<string>(currentColumn) + ".h5", hdf5_binary);
  ret = ret(span(iBeg,iEnd));
  
  for (int i = 0; i < ret.n_rows; i++)
	work[i] = ret(i);
  */

  switch(startType) {
	
  case COLD:
	for (int i = 0; i < size_pars; i++)
	  work[i] = NAN_PAR;
	break;
	
  case HOT:

	// Note column size of mnet is row size of mnStart,
	// thus the use of currentColumn for the row in mnStart.	
	for (int i = 0; i < size_pars; i++)
	  work[i] = mnStart(currentColumn, i);
	
	break;
  }

  //<HERE>
  // fill the remaing work items with returns (oldest to newest)
  for (int r = size_pars; r < size_work; r++) {
	work[r] = mnRet(r-size_pars, currentColumn);
  }
  
  //printf("getting work> currentColumn = %u, rows = %u, cols = %u, r = %u, work.length = %lu\n",
  //		 currentColumn, rows, cols, r, sizeof(work)/sizeof(double));
  
  return 0;
}

//-------------------------------------------------------------------------------
// <SLAVE>

static void do_work_simple(double* work, double* result) {

  // set the returns vector, omitting the first size_pars elements
  gsl_vector* retVec = gsl_vector_calloc(rows);

  bool hasNonzeros = false;
  
  for (int i = size_pars; i < size_work; i++)  {
	gsl_vector_set(retVec, i-size_pars, work[i]);
	if (work[i] != 0.0)
	  hasNonzeros = true;
  }

  if (hasNonzeros == false) {
	gsl_vector_free(retVec);
	return;
  }
  

  garch_struct g;
  
  try {

	// fit GARCH
	g = mygarch::fit_nlopt(retVec, R, M, P, Q);

	// convert structure to vector
	vec x = mygarch::garch_struct_vec(g);

	// set GARCH values to result
	for (int i = 0; i < x.n_rows; i++)
	  result[i] = x(i);

	// starting index of innovation distribution's parameters
	int i0 = x.n_rows; //mygarch::getMessageSize(g);

  
	//
	// 1 = skewed t, 2 = stdAS
	//
	switch (innovType) {
	case 1:
	  {
		// fit skewed t to the residuals
		DistPars d = myskewt::fit_nlopt(g.u);
		result[i0]   = d.gamma;
		result[i0+1] = d.mu;
		result[i0+2] = d.df;
		result[i0+3] = d.sigma;
	  }
	  break;
	case 2:
	  { 
		// convert GSL vector to ARMA vector
		vec gvec(g.u->size);
		for (int i = 0; i < g.u->size; i++)
		  gvec(i) = gsl_vector_get(g.u, i);
	  
		mystable::ts_struct s = mystable::mle_nlopt(gvec, mystable::symAS);
		result[i0]   = s.pars[0]; //alpha
		result[i0+1] = s.pars[1]; //sigma
		result[i0+2] = s.pars[2]; //mu
	  }
	  break;
	default:
	  assert(0 && "Unknown dependence type");
	}
	
  } catch (...) {

	cout << "----> garchfit warning, return vector exported" << endl;

	string errFile = "e_bad_retVec.csv"; //_" + boost::lexical_cast<string>(vnDates(t-1)) + ".csv";
	
	IO::exportGslVector(retVec, errFile);
	
	// free returns vector
	gsl_vector_free(retVec);

	return;
  }

  
  // free the GARCH structure
  mygarch::garch_struct_free(g);
  
  // free returns vector
  gsl_vector_free(retVec);
  
}


//-------------------------------------------------------------------------------
// <SLAVE>

static void do_work_copula(double* work, double* result) {

}


//-------------------------------------------------------------------------------
// <MASTER>
// Set the parameters for this marginal in CopulaPars

static void onMarginalComplete(int idx, double* result) {
  
  // print out the results
  double pct =  ((double)currentColumn)/((double)cols)*100.00;

  if (verbose)
	printf("msj.cpp>> %3.3f | %u | %g | %g  | %g  | %g  | %g  | %g  | %g | %g | %g | %g | %g\n",
		   pct, idx, result[0], result[1], result[2], result[3], result[4], result[5], result[6],
		   result[7], result[8], result[9], result[10]);	
	
  // add result to results matrix
  for (int i = 0; i < size_res; i++) {
	mnResults(i,idx) = result[i];
  }
}


//-------------------------------------------------------------------------------
// <MASTER>

struct RunTimes
{
  double t_dep_est_mc;  // time to estimate MC dependence structure
  double t_dep_est_lut; // time to estimate LUT dependence structure
  double t_shrink_mc;   // time to shrink for MC
  double t_shrink_lut;  // time to shrink for LUT
  double t_chol_mc;     // time for MC cholesky
  double t_chol_lut;    // time for LUT cholesky
  double t_VaR_mc;      // time to estimate VaR with MC given all required parameters
  double t_VaR_lut;     // time to estimate VaR with LUT given all required parameters
  double t_total_mc;    // time for everything with MC
  double t_total_lut;   // time for everything with LUT
};



//-------------------------------------------------------------------------------
// <MASTER>
// This uses the multivariate skew t distribution to directly model the
// skewed t GARCH innovations.

static RunTimes riskForecast_simple(double &VaR_mc, double &VaR_lut) {

  // RunTime structure
  RunTimes rt;
  
  // check sizes
  assert(size_res == mnResults.n_rows);
  assert(cols == mnResults.n_cols);
  //assert(cols == mnRet.n_cols);
  assert(maxRMPQ > 0);
 
  int iCoeffs = (int) mnResults(0,0);
  int iNonCoeffs = mygarch::iNonCoeffs;
  int i0 = iNonCoeffs + iCoeffs + 3*maxRMPQ;
  double shrinkage = 0.0;
 
  mat mnGarchRes = zeros(iL, cols); // mnRet.n_rows
  mat mnGarchPars = zeros(i0, cols);
  
  for (int i = 0; i < cols; i++) {
	mnGarchRes.col(i)  = mnResults(span(i0, size_res-4-1), i);
	mnGarchPars.col(i) = mnResults(span(0, i0-1), i);
  }	

  // Get MC sampled dependent GARCH residuals using specified dependence structure
  if (depStruct != 1) {
	assert(0 && "Unknown or unsupported dependence structure.");
  }

  // Allocate matrix before starting MC timer, as this is used for LUT as well.
  mat C(iS, iS);
  mat resid(nSim, iS); // MC sampled GARCH residuals
  vec nextMeans(iS);   // forecasted GARCH means
  vec nextSigmas(iS);  // forecasted GARCH standard deviations
  vec portRets(nSim);  // MC sampled portfolio returns
  
  switch (innovType) {
  case 1:
	{
	  // set vectors gamma and mu
	  int idxGamma = size_res - 4;
	  int idxMu = size_res - 3;
	
	  vec gamma = mnResults.row(idxGamma).t();
	  vec mu = mnResults.row(idxMu).t();
	  double df = 5.0;
  
	  // Scale gamma - necessary for keeping C positive definite.
	  gamma = gammaScale*gamma;
  
	  // begin timing
	  rt.t_total_mc = MPI::Wtime();

	  // estimate Sigma
	  rt.t_dep_est_mc = MPI::Wtime();
	  myskewt::estimateSigma(mnGarchRes, gamma, mu, df, doShrink, shrinkage, C);
	  rt.t_dep_est_mc = MPI::Wtime() - rt.t_dep_est_mc; // stop the dep. est. timer

	  // cholesky
	  rt.t_chol_mc = MPI::Wtime();
	  C = chol(C);
	  rt.t_chol_mc = MPI::Wtime() - rt.t_chol_mc;

	  // Forecast VaR with Monte Carlo (<TODO> time chol separately)
	  rt.t_VaR_mc = MPI::Wtime(); 
	  myskewt::portSample(nSim, gamma, mu, df, mnGarchPars, wts,
							  C, resid, nextMeans, nextSigmas, portRets);

	  VaR_mc = Stats::quantile(portRets, VaR_epsilon);
	  
	  rt.t_VaR_mc = MPI::Wtime() - rt.t_VaR_mc;

	  // end total time for MC
	  rt.t_total_mc = MPI::Wtime() - rt.t_total_mc;

	  //
	  // VaR with LUT method
	  //
	  // Technically we should add the time it took to compute nextSigmas and nextMeans (~0.6s),
	  // plus the allocation of memory for matrix C, which we update in place below.
	  //
	  if (doLUT == false) {
		break;
	  }

	  // begin total time for LUT
	  rt.t_total_lut = MPI::Wtime();
  
	  // Calculate scaled covariance of residuals (these were nested with chol() call to save memory).
	  //mat mnResScaled = mnGarchRes % repmat(nextSigmas.t(), mnGarchRes.n_rows, 1);
	  //mat mnCov_lut = Stats::cov2para(mnResScaled, shrinkage);

	  rt.t_dep_est_lut = MPI::Wtime();
	  C = Stats::cov2para(mnGarchRes % repmat(nextSigmas.t(), mnGarchRes.n_rows, 1), shrinkage);
	  rt.t_dep_est_lut = MPI::Wtime() - rt.t_dep_est_lut;
	
	  // time chol_lut
	  rt.t_chol_lut = MPI::Wtime();
	  C = chol(C);
	  rt.t_chol_lut = MPI::Wtime() - rt.t_chol_lut;

	  // time VaR_lut
	  rt.t_VaR_lut = MPI::Wtime();

	  C = C.t(); // transpose of chol() is reqiured (this isn't timed in chol())
  
	  // scale gamma by the GARCH forecasted stddevs
	  vec gammaScaled = gamma % nextSigmas;
  
	  double fVaR;
	  lut.lut_var(wts, gammaScaled, nextMeans, df, C, VaR_epsilon, VaR_lut, fVaR);

	  rt.t_VaR_lut = MPI::Wtime() - rt.t_VaR_lut;
	  rt.t_total_lut = MPI::Wtime() - rt.t_total_lut;
  
	}
	break;
  case 2:
	{
	  // set vectors gamma and mu
	  int idxAlpha = size_res - 4;
	  int idxSigma = size_res - 3;
	  int idxMu = size_res - 2;
	
	  vec alpha = mnResults.row(idxAlpha).t();
	  vec sigma = mnResults.row(idxSigma).t();
	  vec mu    = mnResults.row(idxMu).t();

	  // estimate scalar alpha
	  double alphaHat = mystable::assg_alphaEst(alpha);
	  
	  // estimate matrix Sigma
	  mystable::assg_dispersionEst(mnGarchRes, alphaHat, sigma, mu, C);

	  // cholesky
	  C = chol(C);

	  // VaR MC
	  mystable::portSample(nSim, alphaHat, mu, mnGarchPars, wts, C,
						   resid, nextMeans, nextSigmas, portRets);

	  VaR_mc = Stats::quantile(portRets, VaR_epsilon);
	  
	}
	break;
  default:
	{
	  assert(0 && "Unknown innovation type.");
	}
	break;
  }
  

  return rt;
}



//-------------------------------------------------------------------------------
// <MASTER>
// This uses a copula for the dependence structure.

static void riskForecast_copula(double &VaR, double &CVaR) {

}


//-------------------------------------------------------------------------------
// <MASTER>

void onMarginalsComplete() {

  t_garch = MPI::Wtime() - t_garch;

  double date = vnDates(t-1);

  // save results to binary file
  if (margOnly == 1) {
	double t0 = MPI::Wtime();
	string file = reportFile + "/mnResults/" + boost::lexical_cast<string>(date) + ".h5";
	mnResults.save(file, hdf5_binary);
	cout << "----> Exported mnResults to " << file << " (" << MPI::Wtime()-t0<< " sec)\n";
  }

  // time for the entire non-GARCH part
  double t_forc = MPI::Wtime();
  
  // forecast risk
  double VaR_mc, VaR_lut;

  RunTimes rt;
	
  if (margOnly == 1) {

	// skip the VaR compuation
	VaR_mc = 0.0;
	VaR_lut = 0.0;
	
  } else {

	if (verbose) {
	  cout << "----> Forecasting Risk for " << (long) date << ":\n";
	  process_mem_usage();
	}

	if (isSimple)
	  rt = riskForecast_simple(VaR_mc, VaR_lut);
	else
	  riskForecast_copula(VaR_mc, VaR_lut);
	
  }
  
  // Note t is the index of the NEXT returns b/c it was already incremented in the distribute() method.
  double nextRet = sum(mnRetAll(t, span::all).t() % wts); 

  t_period = MPI::Wtime() - t_period;
  t_forc = MPI::Wtime() - t_forc;

  // Append to report file
  fReport = fopen((reportFile + "/forc.csv").c_str(), "a");
  
  //fprintf(fReport,"%u,%9.0f,%g,%g,%u,%g,%g,%g,%g,%g\n",
  //		  t-1, date, VaR_mc, VaR_lut, nextRet, t_garch, t_adjust, t_sample, t_forc, t_period);

  fprintf(fReport, "%u, %9.0f, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
  		  t-1, date, VaR_mc, VaR_lut, nextRet, t_garch,
		  rt.t_dep_est_mc, rt.t_dep_est_lut,
		  rt.t_chol_mc, rt.t_chol_lut,
		  rt.t_VaR_mc, rt.t_VaR_lut,
		  rt.t_total_mc, rt.t_total_lut);

  printf("%u, %9.0f, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
  		  t-1, date, VaR_mc, VaR_lut, nextRet, t_garch,
		  rt.t_dep_est_mc, rt.t_dep_est_lut,
		  rt.t_chol_mc, rt.t_chol_lut,
		  rt.t_VaR_mc, rt.t_VaR_lut,
		  rt.t_total_mc, rt.t_total_lut);
  
  fclose(fReport);

  // cout << setw(10) << t-1
  // 	   << setw(10) << (int) date
  // 	   << setw(15) << VaR_mc
  // 	   << setw(15) << VaR_lut
  // 	   << setw(15) << nextRet
  // 	   << setw(15) << t_garch
  // 	   << setw(15) << t_adjust
  // 	   << setw(15) << t_sample
  // 	   << setw(15) << t_forc
  // 	   << setw(15) << t_period
  // 	   << setw(5)  << VaR_total << endl; 
  
}
