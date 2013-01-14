#ifndef __MYGARCH_PAR_HPP__
#define __MYGARCH_PAR_HPP__

#include <iostream> 
#include <vector>
#include <cassert>

#include <mpi.h>
#include "elemental.hpp"
#include <armadillo>

#include "MasterSlave.hpp"


using namespace std;
using namespace elem;
//using namespace arma;

typedef double R;

class mygarch_par: public MasterSlave {

  int currentAsset;
  MPI_Comm comm;
  DistMatrix<double,STAR,VC> returns;
  DistMatrix<double,STAR,VC> results;

public:
  mygarch_par(MPI_Comm comm, DistMatrix<double,STAR,VC> returns, DistMatrix<double,STAR,VC> &results);
  
  bool get_next_work_item(double *work);
  void do_work(double *work, double *result, int size_work, int size_res);
  void onWorkComplete(double *result);
  void onAllWorkComplete();

  //void fit(DistMatrix<double> returns, DistMatrix<double> &results);
};

#endif
