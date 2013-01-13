  /*
  int rank, size, namelen, xranks[] = { 0 };
  int send_val, recv_val, send_val2, recv_val2;
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Group mpi_group_world, group_slaves;
  MPI_Comm comm_slaves;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Get_processor_name(processor_name, &namelen);

  MPI_Comm_group(MPI_COMM_WORLD, &mpi_group_world);
  MPI_Group_excl(mpi_group_world, 1, xranks, &group_slaves);
  MPI_Comm_create(MPI_COMM_WORLD, group_slaves, &comm_slaves);

  if (comm_slaves != MPI_COMM_NULL)
	MPI_Comm_free(&comm_slaves);

  MPI_Group_free(&group_slaves);
  MPI_Group_free(&mpi_group_world);
  MPI_Finalize();
  */

//-------------------------------------------------------------------------------
// Test on a small universe

void test_smallWorld() {

  //----------------------------------------------

  int nSim = 10000;
  string fname = "/Users/jimmiegoode/Dropbox/Projects/Glimm/data/build_CSI_indexData_SPX.csv";
  //mnRetSmall.csv";
  const int R = 1;
  const int M = 1;
  const int P = 1;
  const int Q = 1;
  double VaR_epsilon = 0.01;

  EST_METHOD method = NMSIMPLEX;
  START_TYPE startType = COLD;

  //----------------------------------------------
  
  gsl_matrix* mnRet = IO::importCSVmatrix(fname);
  int iObs = mnRet->size1;
  int iS = mnRet->size2;

  std::vector<garch_struct> gs;
  std::vector<DistPars> ds;

  gsl_vector* vnRet = gsl_vector_alloc(iObs);

  // fit ARMA-GARCH and distribution parameters
  for (int i = 0; i < iS; i++) {
	
	printf("Progress: %3.1f, i = %u\n", ((double) (i+1)/iS)*100.0, i);
	//printProgressBar((int) boost::math::round((i+1)/iS*100.0));
	
	gsl_matrix_get_col(vnRet, mnRet, i);

	// fit GARCH
  	garch_struct g = mygarch::fit_nlopt(vnRet, R, M, P, Q);

	// fit skewed t to the residuals
	DistPars d;
	
	if (method == NMSIMPLEX)
	  d = myskewt::skewedstudenttfit_nmsimplex(startType, g.u);
	else if (method == BFGS)
	  d = myskewt::skewedstudenttfit_bfgs(startType, g.u);
	else
	  assert(0 && "Bad method");
	
	//mygarch::garch_struct_print(g);
	//myskewt::DistPars_print(d);
    
	gs.push_back(g);
	ds.push_back(d);
  }

  // fit copula
  vec gamma(iS);
  vec mu(iS);
  mat res(iObs,iS); // ARMA-GARCH residuals
  
  for (int i = 0; i < iS; i++) {
	gamma[i] = ds[i].gamma;
	mu[i] = ds[i].mu;
	
	for (int t = 0; t < iObs; t++)
	  res(t,i) = gsl_vector_get(gs[i].u, t);
  }

  // compute copula matrix. Z ~ N(0,C) in stochastic form.
  double df = 5.0;
  mat C1  = cov(res)*((df - 2.0)/df) - gamma*gamma.t() * (2.0*df)/((df - 2.0)*(df - 4.0));

  // Adjust C to be semi-possitive
  cx_vec eigval;
  cx_mat eigvec;

  // <TODO> the faster eig_sym can be used after forcing C to be symmetric
  eig_gen(eigval, eigvec, C1); 

  // P*D*P'
  cx_mat C2 = eigvec * diagmat(abs(eigval)) * eigvec.t();

  // convert complex-compatible matrix (cx_mat) to real matrix (mat)
  mat C = conv_to<mat>::from(C2);

  // simulate Monte Carlo sample from the multivariate residual distribution
  mat resid = myskewt::mvskewtrnd(gamma, mu, df, C, nSim);

  // forecast stock returns using ARMA-GARCH with copula residuals
  mat stockRets(nSim, iS);

  for (int i = 0; i < iS; i++) {
	stockRets.col(i) = mygarch::forecast(gs[i], resid.col(i));
  }

  // set equal portfolio weights
  vec wts = 1.0/((double) iS) * ones(iS);

  // portfolio returns for each path
  vec portRets = stockRets*wts;

  // GSL quantile for VaR
  gsl_vector* pr = gsl_vector_alloc(portRets.n_rows);

  for (int i = 0; i < portRets.n_rows; i++)
	gsl_vector_set(pr, i, portRets(i));

  gsl_sort(pr->data, pr->stride, pr->size);
  double VaR = gsl_stats_quantile_from_sorted_data(pr->data, pr->stride, pr->size, VaR_epsilon);
  gsl_vector_free(pr);
  
  cout << "\n\nVaR = " << VaR << "\n\n";
  
  // double median = gsl_stats_median_from_sorted_data (pr->data, pr->stride, pr->size);
  // double upperq = gsl_stats_quantile_from_sorted_data (pr->data, pr->stride, pr->size, 0.75);
  // double lowerq = gsl_stats_quantile_from_sorted_data (pr->data, pr->stride, pr->size, 0.25);
  // printf ("The median is %g\n", median);
  // printf ("The upper quartile is %g\n", upperq);
  // printf ("The lower quartile is %g\n", lowerq);
  
  // Print debug info
  if (debug) {
	
	// save
	string p = "/Users/jimmiegoode/Dropbox/Projects/Glimm/Verification/";
	gamma.save(p + "e_gamma", csv_ascii);
	mu.save(p + "e_mu", csv_ascii);
	res.save(p + "e_res", csv_ascii);
	C1.save(p + "e_C1", csv_ascii);
	C2.save(p + "e_C2", csv_ascii);
	C.save(p + "e_C", csv_ascii);

	resid.save(p + "e_resid", csv_ascii);
	stockRets.save(p + "e_stockRets", csv_ascii);
	wts.save(p + "e_wts", csv_ascii);
	portRets.save(p + "e_portRets", csv_ascii);
	
	// print
	// gamma.print("gamma = ");
	// mu.print("mu = ");
	// C.print("C = ");
	// //C2.print("Cadj = ");
	// wts.print("wts = ");
	// stockRets.print("stockRets = ");
	// portRets.print("portRets = ");
  }
 
  //printf("\n\nFreeing data structures...\n\n");
  
  // free structures
  for (int i = 0; i < iS; i++) {
	mygarch::garch_struct_free(gs[i]);
	//gsl_vector_free(ds[i].y);
  }
  
  // free return vectors
  gsl_matrix_free(mnRet);
  gsl_vector_free(vnRet);
}


/*
	index_zero = max(find(cdfi < 1.0e-5));
	index_one = min(find(cdfi > 1 - 1.0e-5));

	cdfi(index_one + 1:end) = [];
	x(index_one + 1:end) = [];

	cdfi(1:index_zero - 1) = [];
	x(1:index_zero - 1) = []; 

	if min(pvalue) < min(cdfi) | max(pvalue) > max(cdfi)
    in = find(pvalue > max(cdfi));
    pvalue(in) = max(cdfi);
    in = find(pvalue < min(cdfi));
    pvalue(in) = min(cdfi);
    warning('pvalue outside cdfi grid. Corrected.');
	end    

  */

  /*
  uvec idx0 = find(cdfi < tol, 1, "last");
  if (idx0.n_rows > 0) {
	assert(idx0.n_rows == 1);
	cdfi.shed_rows(0, idx0(0));
	x.shed_rows(0, idx0(0));
  }
  
  uvec idx1 = find(cdfi > 1.-tol, 1, "first");
  if (idx1.n_rows > 0) {
	assert(idx1.n_rows == 1);
	cdfi.shed_rows(idx1(0), iN-1);
	x.shed_rows(idx1(0), iN-1);
  }
  assert(cdfi.n_rows == x.n_rows);
  iN = cdfi.n_rows;
  */
  /* 
  cout << "Saving ...\n";
  cdfi.save("/Users/jimmiegoode/Documents/Glimm/Toolbox/alphastable/testing/e_cdfi.csv", csv_ascii);
  x.save("/Users/jimmiegoode/Documents/Glimm/Toolbox/alphastable/testing/e_x.csv", csv_ascii);
  */

  // vector<double> cdfi2;
  // vector<double> x2;

  // for (int i = 1; i < iN; i++) {
  // 	if (cdfi(i) - cdfi(i-1) < 0.0) {
  // 	  continue;
  // 	}
  // 	if (cdfi(i) < 0.0 || cdfi(i) > 1.0) {
  // 	  continue;
  // 	}
  // 	cdfi2.push_back(cdfi(i));
  // 	x2.push_back(x(i));
  // }
  
  // iN = cdfi2.size();
  // double vnCdfi[iN];
  // double vnX[iN];

  // std::copy(cdfi2.begin(), cdfi2.end(), vnCdfi);
  // std::copy(x2.begin(), x2.end(), vnX);

  //uvec idxKeep = find(cdfi > tol && cdfi < 1.-tol);

/*
//-------------------------------------------------------------------------------

double negLLF_stdCTS(unsigned n, const double *x, double *grad, void *data) {

  fxcalls++;
  
  gsl_vector *yy = (gsl_vector *) data;
  vec pdf = pdf_FFT(yy->size, 3, yy->data, x, mystable::stdCTS);
  double negLLF = -1.0*as_scalar(sum(log(pdf)));

  cout << arr2str(n,x) << ", negLLF = " << negLLF << " --> " << checkNegLLF(negLLF) << endl;
  
  return checkNegLLF(negLLF);
}

//-------------------------------------------------------------------------------

double negLLF_stdNTS(unsigned n, const double *x, double *grad, void *data) {

  fxcalls++;
  
  gsl_vector *yy = (gsl_vector *) data;
  vec pdf = pdf_FFT(yy->size, 3, yy->data, x, mystable::stdNTS);
  double negLLF = -1.0*as_scalar(sum(log(pdf)));
  cout << arr2str(n,x) << ", negLLF = " << negLLF << " --> " << checkNegLLF(negLLF) << endl;

  return checkNegLLF(negLLF);
}

//-------------------------------------------------------------------------------

double negLLF_stdAS(unsigned n, const double *x, double *grad, void *data) {

  fxcalls++;
  
  gsl_vector *yy = (gsl_vector *) data;
  vec pdf = pdf_FFT(yy->size, 2, yy->data, x, mystable::stdAS);
  double negLLF = -1.0*as_scalar(sum(log(pdf)));
  cout << arr2str(n,x) << ", negLLF = " << negLLF << " --> " << checkNegLLF(negLLF) << endl;

  return checkNegLLF(negLLF);
}
*/

/*

// cx_vec chf_CTS(cx_vec u, double alpha_r, double C_r, double lmp_r, double lmm_r, double mm_r) {
 
//   double mReal = mm_r - xGamma(1.0 - alpha_r) * C_r * (std::pow(lmp_r, alpha_r-1.0)
// 													   - std::pow(lmm_r, alpha_r-1.0));
//   cx_double m(mReal,0);  
//   cx_double onei(0,1);
//   cx_double alpha(alpha_r,0);
//   cx_double C(C_r,0);
//   cx_double lmp(lmp_r,0);
//   cx_double lmm(lmm_r,0);

//   //cx_double gam(xGamma(-alpha_r), 0);
//   //cx_double pow1(pow(lmp_r - onei % u, alpha), 0);

//   /*
//   cx_vec phi = arma::exp(dt % (onei % u % m + C % xGamma(-alpha_r) % (
// 					  pow(lmp_r - onei % u, alpha) - pow(lmp,alpha)
// 					  + pow(lmm + onei % u, alpha) - pow(lmm, alpha))));
//   */

//   cx_vec phi = arma::exp(dt * (onei * u * m + C * xGamma(-alpha_r) * (
// 					  pow(lmp_r - onei * u, alpha) - pow(lmp,alpha)
// 					  + pow(lmm + onei * u, alpha) - pow(lmm, alpha))));
  
  
//   return phi;
// }




  
  // T_GARCH case

  double startTime = MPI::Wtime();
  
  garchpars pars(1,1,1,1);
  pars.initPars();
  pars.setReturnData(retVec);
  GARCH::fit(pars);

  double runTime = MPI::Wtime() - startTime;
  
  result[0] = pars.C;
  result[1] = gsl_vector_get(pars.AR,0);
  result[2] = gsl_vector_get(pars.MA,0);
  result[3] = pars.K;
  result[4] = gsl_vector_get(pars.GARCH,0);
  result[5] = gsl_vector_get(pars.ARCH,0);
  result[6] = pars.DoF;
  result[7] = pars.iters;
  result[8] = pars.LL;
  result[9] = pars.exitFlag;
  result[10] = runTime;
  */  

  /*
  // set initial guess for parameters
  DistPars distPars0;
  switch(marginalType) {
  case ST:
	switch(startType) {
	case COLD:
	  distPars0.gamma = NAN_PAR;
	  distPars0.mu    = NAN_PAR;
	  distPars0.df    = NAN_PAR;
	  distPars0.sigma = NAN_PAR;
      break;
	case HOT:
	  distPars0.gamma = work[ST_GAMMA];
	  distPars0.mu    = work[ST_MU];
	  distPars0.df    = work[ST_DF];
	  distPars0.sigma = work[ST_SIGMA];
	  break;
	}

	//printf("\ndistPars: ");
	//printPars(distPars0);
	//printf("work[...] = %g, %g, %g, %g\n", work[ST_GAMMA], work[ST_MU], work[ST_DF], work[ST_SIGMA]);
	//exportDoubleArray(work, "../data/reports/work.csv", size_work);

	// set global fit variables
	setData(retVec, distPars0);
	
	// fit the marginal
	DistPars distPars;
	switch (estMethod) {
	case BFGS:
	   distPars = skewedstudenttfit_bfgs(startType); 
	   break;
	case NMSIMPLEX:
	  distPars = skewedstudenttfit_nmsimplex(startType);
	  break;
	}
	  
	result[0] = distPars.gamma;
	result[1] = distPars.mu;
	result[2] = distPars.df;
	result[3] = distPars.sigma;
	result[4] = distPars.iters;
	result[5] = distPars.LLF;
	result[6] = distPars.status;
    break;
	
  case AS:
	break;
  }
  */
/*
//<N>
static void multivariatenormalrnd(gsl_rng* rng , gsl_matrix* A,int nCol,int nSim);
static void skewedmvtrnd(gsl_rng * rng ,double * mu, double * gamma, int df,gsl_matrix * Z);
*/


/*
// structure to hald multivariate copula parameters
struct CopulaPars {
int dim;
double * gamma;  // vector of marginal gamma's
double * mu;     // vector of marginal mu's 
double * df;     // vector for marginal degrees of freedom
double * sigma;  // vector for marginal sigmas
double ** C;     // matrix for copula covariance
double * iters;
double * LLF;
double * status;
};
*/


/*
// Declare the supported options.
po::options_description desc("Allowed options");
desc.add_options()
("help", "produce help message")
("m", po::value<int>(),
"sets number of marginals to estimate, default is -1 for all")
("f", po::value<string>(), "name of file to export results to");

po::variables_map vm;
po::store(po::parse_command_line(argc, argv, desc), vm);
po::notify(vm);    

if (vm.count("help")) {
cout << desc << "\n";
return 1;
}

if (vm.count("m")) {
cout << "Marginal count was set to " << vm["m"].as<int>() << endl;
} else {
cout << "Marginal count was not set.\n";
}
if (vm.count("f")) {
cout << "Report file was set to " << vm["f"].as<string>() << endl;
reportFile = vm["f"].as<string>();
} else {
reportFile = "report.csv";
}
*/

/*** <N><TODO> merge with other methods ***/
/*
//setup of gsl rng
  
srand(time(NULL));

gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
  
gsl_rng_env_setup();
  
gsl_rng_set(rng, time(NULL));   
  
static void multivariatenormalrnd(rng , mnRet,nCol,nSim);
  
static void skewedmvtrnd(rng ,&mu, &gamma,df,Z);
*/


  /*
	gsl_matrix* mnIn = importCSVmatrix(dataFile);
	mnRet = gsl_matrix_alloc(251, 2);
	for (int i = 0; i < mnRet->size1; i++) {
	gsl_matrix_set(mnRet, i, 0, gsl_matrix_get(mnIn, i, 1));
	gsl_matrix_set(mnRet, i, 1, gsl_matrix_get(mnIn, i, 2));
	}
  */
  
  /*
  // int idxs[3] = {1993, 1222, 1219};
  gsl_matrix* mnIn = importCSVmatrix("../data/ret_252_days_2127_stocks.csv");
  mnRet = gsl_matrix_alloc(251, 3);

  for (int i = 0; i < mnRet->size1; i++) {
  gsl_matrix_set(mnRet, i, 0, gsl_matrix_get(mnIn, i, idxs[0]));
  gsl_matrix_set(mnRet, i, 1, gsl_matrix_get(mnIn, i, idxs[1]));
  gsl_matrix_set(mnRet, i, 2, gsl_matrix_get(mnIn, i, idxs[2]));
  }
  */
