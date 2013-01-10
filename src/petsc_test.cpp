
//#include <petscblaslapack.h>
#include <petscmat.h>
//#include <petsc.h>

#include <mpi.h>

#include <iostream>


/*
  
tacc-web.austin.utexas.edu/veijkhout/public_html/taccnotes/petsclinking/petsclinking.html
  
FROM the petsc install dir:

/Users/jimmiegoode/Documents/Glimm/lib/petsc-3.3-p3

ISSUE this command to get all the linking and include flags:

make getlinklibs
make getincludedirs


mpic++  -L/Users/jimmiegoode/Documents/Glimm/lib/petsc-3.3-p3/arch-darwin-c-debug/lib -L/Users/jimmiegoode/Documents/Glimm/lib/petsc-3.3-p3/arch-darwin-c-debug/lib -lpetsc -lX11 -lpthread -llapack -lblas -L/opt/local/lib -L/usr/lib/clang/3.1/lib/darwin -lclang_rt.osx -lmpi_f90 -lmpi_f77 -lm -lgfortran -L/opt/local/lib/gcc45/gcc/x86_64-apple-darwin11/4.5.4 -L/opt/local/lib/gcc45 -lm -lgcc_s.10.5 -lgcc_ext.10.5 -ldl -lmpi -lSystem -L/usr/bin/../lib/clang/3.1/lib/darwin -lclang_rt.osx -ldl -I/Users/jimmiegoode/Documents/Glimm/lib/petsc-3.3-p3/include -I/Users/jimmiegoode/Documents/Glimm/lib/petsc-3.3-p3/arch-darwin-c-debug/include -I/opt/local/include/openmpi petsc_test.cpp -o testp

mpirun -np 2 testp

*/

using namespace std;


static char help[] = "";

int main(int argc, char **argv) {

  int myrank, ntasks;

   // initialize MPI
  MPI_Init(&argc, &argv);
  
  // find out my identity in the default communicator 
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  
  cout << "hello from rank " << myrank << endl;


  // Petsc
  PetscInitialize(&argc, &argv, (char *)0, help);

  

  
  MPI_Finalize();

  return 0;
}

