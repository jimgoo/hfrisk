#ifndef __MSJ_H__
#define __MSJ_H__

#include <stdio.h>
#include <vector>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <mpi.h>
#include <armadillo>
//#include <boost/filesystem.hpp>
#include "elemental.hpp"

#include "IO.hpp"
#include "Stats.hpp"
#include "mygarch.hpp"
#include "myskewt.hpp"
#include "Constants.h"
#include "mystable.hpp"
#include "myskewt_lut.hpp"

using namespace arma;
using namespace elem;

static void master(void);
static void slave(void);
static int get_next_work_item(double* work);
static void do_work_simple(double* work, double* result);
static void do_work_copula(double* work, double* result);
static bool distribute(int rank);
static void onMarginalsComplete();
static void onMarginalComplete(int idx, double * result);

#endif
