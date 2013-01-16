/*
test code with WORKING passing of grid in nested functions



*/



#include "elemental.hpp"
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <cmath>

using namespace std;
using namespace elem;

typedef double R;
typedef float F;

void Foo( Grid** grid, mpi::Comm comm );

void Bar( mpi::Comm comm );


int nprocs;
int commRank;

int main( int argc, char* argv[] )
{
  Initialize( argc, argv );

  mpi::Comm comm = mpi::COMM_WORLD;
  commRank = mpi::CommRank( comm );
  nprocs = mpi::CommSize(comm );
  const int n = Input("--size","NxN size of Matricies",10);
  const int df = Input("--df","Degrees of freedom",5);
  const int nSims = Input("--nSims","Number of Monte Carlo Simulations",10);
  const bool print = Input("--print","print matrices?",false);
  ProcessInput();
  PrintInputReport();

  Bar(comm);

     
  Finalize();
  return 0;
 }




void Foo( Grid** grid, mpi::Comm comm )
{
    *grid = new Grid( comm );
    DistMatrix<R> M (10,10,**grid);
}

void Bar( mpi::Comm comm )
{
    Grid* grid;
    Foo( &grid, comm );
    DistMatrix<R> N (10,10,*grid);
    std::cout << "My rank: " << grid->Rank() << std::endl;
    std::cout<<" N.Width(): = " << N.Width() <<std::endl;
    delete grid;
}
