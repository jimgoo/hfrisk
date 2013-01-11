
// g++ -o Csv2Binary Csv2Binary.cpp -I/Users/jimmiegoode/Documents/Glimm/lib/armadillo-3.4.2/include -framework Accelerate

// g++ -o Csv2Binary Csv2Binary.cpp -I/nfs/user03/copula/20120323/lib/armadillo-3.4.1/include -I/nfs/user03/copula/20120323/lib/OpenBLAS-v0.2.3-0/xianyi-OpenBLAS-48f075c/install/include -L/nfs/user03/copula/20120323/lib/OpenBLAS-v0.2.3-0/xianyi-OpenBLAS-48f075c/install/lib -lopenblas /usr/lib/liblapack.so -lgfortran

// g++ -o Csv2Binary Csv2Binary.cpp -I/Users/jimmiegoode/Documents/Glimm/lib/hdf5-1.8.9/install/include -L/Users/jimmiegoode/Documents/Glimm/lib/hdf5-1.8.9/install/lib -I/Users/jimmiegoode/Documents/Glimm/lib/armadillo-3.4.2/include -lhdf5 -framework Accelerate 


#include <armadillo>
#include <iostream>
#include <cassert>

using namespace std;
using namespace arma;

int main(int argc, char** argv) {

  if (argc < 3)
	assert(0 && "Too few arguments");

  //int testCase = atoi(argv[1]);
  string fin  = argv[1];
  string fout = argv[2];

  mat x;
  x.load(fin, csv_ascii);
  x.save(fout, arma_binary);

  return 0;
}

	
