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

#include "IO.h"
#include "Stats.hpp"
#include "mygarch.hpp"
#include "myskewt.hpp"
#include "Constants.h"
#include "mystable.hpp"
#include "myskewt_lut.hpp"

using namespace arma;

static void master(void);
static void slave(void);
static int get_next_work_item(double* work);
static void do_work_simple(double* work, double* result);
static void do_work_copula(double* work, double* result);
static bool distribute(int rank);
static void onMarginalsComplete();
static void onMarginalComplete(int idx, double * result);
//static void printProcMap(int * map, int size);
//tatic void freeAllMemory();

#endif


//INITIAL MEM >> VM: 190472.00, RSS: 5044.00
//               VM: 198612.00, RSS: 10052.00

//INITIAL MEM >> VM: 190300.00, RSS: 5044.00
//               VM: 198612.00, RSS: 10040.00

/*
// print the processor to data index mapping
static void printProcMap(int* map, int size) {

  FILE* out;
  out = fopen("procMap.txt", "w+");
  
  fprintf(out, "\n============== ProcMap ==============\n");
  for (int i = 0; i < size; i++) {
	fprintf(out,"proc %u -> %u\n", i, map[i]);
  }
  
  fclose(out); 
}
*/
