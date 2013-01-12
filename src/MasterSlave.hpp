#ifndef __MASTER_SLAVE_HPP__
#define __MASTER_SLAVE_HPP__

#include <iostream>

#include <mpi.h>


using namespace std;

class MasterSlave {

public:

  int iters;
  
  virtual bool get_next_work_item(double *work);
  virtual void do_work(double *work, double *result, int size_work, int size_res);
  virtual void onWorkComplete(double *result);
  virtual void onAllWorkComplete();
  
  void slave(int size_work, int size_res);
  void master(int size_work, int size_res);
  void distribute(int rank);

};

#endif
