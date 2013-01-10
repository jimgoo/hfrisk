/*
  Author: Jimmie Goode
  Created: 2012-12-28
*/

#ifndef __MYSKEWT_LUT_HPP__
#define __MYSKEWT_LUT_HPP__

// std
#include <iostream>

// external
#include <armadillo>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h> // interp for pdfs
#include <gsl/gsl_interp.h> // interp for pdfs

// local
#include "myskewt.hpp"


using namespace std;


class myskewt_lut {
  
  vec beta1;
  vec beta2;
  double df;
  int d;
  int nSim;
  vec xGrid;
  mat table;
  
public:

  myskewt_lut();
  ~myskewt_lut();
  
  static void tester(int c);

  void buildTable(vec beta1, vec beta2, double df, int d, int nSim, vec xGrid);
  
  void save(string path, file_type type, string ext);
  void load(string folder, file_type type, string ext);
  
  vec lookupPDF(double nBeta1, double nBeta2);
  void print();
  void lut_solveG(double epsilon, double nBeta1, double nBeta2, double &fVaR);
  void lut_var(vec w, vec beta, vec mu, double df, mat A, double epsilon, double &VaR, double &fVaR);
  
};
  
#endif
