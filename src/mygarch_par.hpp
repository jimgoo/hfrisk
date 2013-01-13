#ifndef __MYGARCH_PAR_HPP__
#define __MYGARCH_PAR_HPP__

#include <iostream>
#include <vector>
#include <cassert>

#include <mpi.h>
#include "elemental.hpp"

#include "MasterSlave.hpp"


using namespace std;
using namespace elem;

class mygarch_par: public MasterSlave {

  int a;

public:
  bool get_next_work_item(double *work);
  void do_work(double *work, double *result, int size_work, int size_res);
  void onWorkComplete(double *result);
  void onAllWorkComplete();
  
};

#endif
