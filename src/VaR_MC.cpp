
#include "elemental.hpp"
#include <armadillo>
#include <cstdio>

using namespace std;
using namespace elem;
using namespace arma;

// Typedef our real type to 'R' for convenience
typedef double R;

double t_chol;

int main( int argc, char* argv[] ) {
  
  Initialize( argc, argv );

  mpi::Comm comm = mpi::COMM_WORLD;
  const int commRank = mpi::CommRank( comm );

  try {
	
	  //const int rows = Input("-rows", "number of rows (observations)", 100);
	  //const int cols = Input("-cols", "number of cols (assets)", 100);
	  const int n = Input("-n", "number of assets", 100);
	  const int nSim = Input("-nSim", "number of simulations", 100000);
	  const bool print = Input("-print", "print matrices?", false);
	  
	  ProcessInput();
	  //PrintInputReport();

	  Grid g( comm );

	  DistMatrix<R> C( n, n, g );

	  vec diag = linspace(1.0, 1.0/n, n);

	  for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
		  double val = 1.0/(i + j + 1.0);
		  if (i == j)
			val += diag(i);
		  C.Set(i, j, val);
		}
	  }

	  if (print) {
		C.Print("C = ");
	  }

	  /*
	  MPI_Barrier(MPI_COMM_WORLD);
	  if (commRank == 0) {
		t_chol = MPI::Wtime();
	  }

	  Cholesky(LOWER, C);

	  MPI_Barrier(MPI_COMM_WORLD);
	  if (commRank == 0) {
		t_chol = MPI::Wtime() - t_chol;
		cout << "Cholesky runtime = " << t_chol << endl;
	  }
	  */

	  DistMatrix<C> U(g);
	  Uniform(U);

	  
	  
        
    } catch( ArgException& e ) {
        // There is nothing to do
    } catch( exception& e ) {
        ostringstream os;
        os << "Process " << commRank << " caught exception with message: "
           << e.what() << endl;
        cerr << os.str();
#ifndef RELEASE
        DumpCallStack();
#endif
    }

    Finalize();
    return 0;
}
