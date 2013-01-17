/*
  Convert CSV file to HDF5 file using Armadillo.
*/

#include <armadillo>
#include <iostream>
#include <cassert>

using namespace std;
using namespace arma;

int main(int argc, char** argv) {

  if (argc < 3)
	assert(0 && "Too few arguments");

  string fin  = argv[1];
  string fout = argv[2];

  mat x;
  x.load(fin, csv_ascii);
  x.save(fout, hdf5_binary);

  return 0;
}

	
