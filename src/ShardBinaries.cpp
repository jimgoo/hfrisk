
// g++ -o ShardBinaries ShardBinaries.cpp -I/Users/jimmiegoode/Documents/Glimm/lib/armadillo-3.4.2/include -framework Accelerate

// g++ -o ShardBinaries ShardBinaries.cpp -I/nfs/user03/copula/20120323/lib/armadillo-3.4.1/include -I/nfs/user03/copula/20120323/lib/OpenBLAS-v0.2.3-0/xianyi-OpenBLAS-48f075c/install/include -L/nfs/user03/copula/20120323/lib/OpenBLAS-v0.2.3-0/xianyi-OpenBLAS-48f075c/install/lib -lopenblas /usr/lib/liblapack.so -lgfortran

// g++ -o ShardBinaries ShardBinaries.cpp -I/Users/jimmiegoode/Documents/Glimm/lib/hdf5-1.8.9/install/include -L/Users/jimmiegoode/Documents/Glimm/lib/hdf5-1.8.9/install/lib -I/Users/jimmiegoode/Documents/Glimm/lib/armadillo-3.4.2/include -lhdf5 -framework Accelerate 


#include <armadillo>
#include <iostream>
#include <cassert>

#include "IO.h"

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{

  if (argc < 3)
	assert(0 && "Too few arguments");

  string fin  = argv[1];
  string fout = argv[2];

  mat X;

  /*
  if (0) {
	  cout << "\n\nGenerating\n\n";
	  X = randu<mat>(1000,3);
	  X.save(fout, hdf5_binary);
  } else {
	  cout << "\n\nLoading\n\n";
	  X.load(fout, hdf5_binary);
  }
  
  //X.print("X = ");
  cout << X.n_rows << endl;
  cout << X.n_cols << endl;

  
  sum(X).print("sum(X) = ");
  //std(X).print("std(X) = ");
  */

  X.load(fin, arma_binary);

  system(("mkdir " + fout).c_str());
  
  for (int i = 0; i < X.n_cols; i++) {	  
	vec x = X.col(i);
	x.save(fout + "/EQ_" + boost::lexical_cast<string>(i) + ".h5", hdf5_binary);
	//+ ".abin", arma_binary); //+ ".h5", hdf5_binary);
  }
  
  return 0;
}

	
