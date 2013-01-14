/*
  Author: Jimmie Goode
  Created: 2012-09-01
*/

// Utilities for importing and exporting data.

#include "IO.hpp"

using namespace std;

//-------------------------------------------------------------------------------
// Import a multicolumn CSV file to GSL matrix

gsl_matrix* IO::importCSVmatrix(string file) {
  
  string line;
  ifstream myfile(file.c_str());
  vector<double> data;
  vector<string> toks;

  typedef boost::tokenizer< boost::escaped_list_separator<char> > Tokenizer;

  if (!myfile.is_open()) {
	cout << "Unable to open file: " + file << endl;
  }

  // get the data dimensions
  int r = 0;
  int c = 0;

  while(myfile.good()) {
	getline (myfile, line);
	if (line == "") {
	  break;
	}
	r++;
	if (r == 1) {
	  Tokenizer tok(line);
	  toks.assign(tok.begin(), tok.end());
	  c = toks.size();
	}
  }
  
  // Allocate matrix and parse
  gsl_matrix *mnData = gsl_matrix_alloc(r,c);
  
  myfile.clear();
  myfile.seekg(0, ios::beg);

  int r2 = 0;

  while (myfile.good()) {
	getline (myfile,line);
	if (line == "") {
	  break;
	}
	r2++;
	Tokenizer tok(line);
	toks.assign(tok.begin(), tok.end());

	for (int i = 0; i < c; i++) {
	  boost::trim(toks[i]);
	  gsl_matrix_set(mnData, r2-1, i, boost::lexical_cast<double>(toks[i]));
	}
  }

  cout << "IO.cpp> Imported " << boost::lexical_cast<int>(r) << " x "
  	   << boost::lexical_cast<int>(c) << " matrix from '" << file << "'" << endl;

  myfile.close();

  return mnData;
}

//-------------------------------------------------------------------------------
// Export a gsl_matrix

void IO::exportGslMatrix(const gsl_matrix *mat, const string fname) {

  FILE* fout;
  fout = fopen(fname.c_str(), "w+");

  for (int r = 0; r < mat->size1; r++) {
	for (int c = 0; c < mat->size2; c++) {
	  fprintf(fout, "%5.10f,", gsl_matrix_get(mat, r, c));
	}
	fprintf(fout,"\n");
  }
  fclose(fout);

  if (1) {
  	cout << "IO.cpp> Exported " << boost::lexical_cast<int>(mat->size1) << " x "
  		 << boost::lexical_cast<int>(mat->size2) << " matrix to '" << fname << "'" << endl;
  } else {
	string msg = "IO.cpp> Export failed for " + fname;;
	cout << msg << endl;
	//assert(0 && msg)
  }
}


//-------------------------------------------------------------------------------
// utility function to see vectors

void IO::exportGslVector(gsl_vector* vec, string fname) {
    
  FILE* fout = fopen(fname.c_str(),"w+");
  gsl_vector_fprintf(fout, vec, "%f");
  fclose(fout);

  // if (ret == 0) {
  // 	cout << "IO.cpp> Exported length " << boost::lexical_cast<int>(vec->size)
  // 		 << " vector to '" << fname << "'" << endl;
  // } else {
  // 	string msg = "IO.cpp> Export failed for " + fname;
  // 	cout << msg << endl;
  // 	//assert(0 && msg);
  // }
}




//-------------------------------------------------------------------------------

void IO::printMatrix(gsl_matrix* mat) {
  for (int r = 0; r < mat->size1; r++) {
	for (int c = 0; c < mat->size2; c++) {
	  printf("%5.10f ", gsl_matrix_get(mat, r, c));
	}
	printf("\n");
  }
}


//-------------------------------------------------------------------------------
/*
// import skewed copula parameters
boost::tuple<gsl_matrix *,gsl_matrix *> IO::importCopulaPars(string fbase) {
  gsl_matrix *C = importCSVmatrix(fbase + ".cov");
  gsl_matrix *P = importCSVmatrix(fbase + ".par");

  return boost::tuple<gsl_matrix *, gsl_matrix *>(C, P);
}
*/


//-------------------------------------------------------------------------------
/*
void IO::exportDoubleArray(double * arr, string filename, int length) {

  FILE* outfile = fopen(filename.c_str(),"w+");

  for (int i = 0; i < length; i++) {
	fprintf(outfile,"%+10.20f \n", arr[i]);
  }  
  fclose(outfile);
}
*/


//-------------------------------------------------------------------------------
// Data is distributed by columns across each process
/*
DistMatrix<double, STAR, VC> IO::arma2distMat(string fname, arma::file_type type, elem::Grid grid) {

  mat armaRet;
  armaRet.load(fname, type);

  const int r = armaRet.n_rows;
  const int c = armaRet.n_cols;

  DistMatrix<double, STAR, VC> dmRet(r, c, grid);

  return dmRet;
}
*/
