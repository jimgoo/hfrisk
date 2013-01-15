/*
  Skewed T MVT RND using elemental
  X:= mu + Gamma*W + Z*sqrt(W)

  Alpha Stable RNG using elemental


  <TODO>
  * repmatrix function - needs optimization to use just one DistMatrix instead of two
  * Set seeds for production code to be randomized
  * MASSIVE memory deallocation

*/


#include "elemental.hpp"
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <cmath>

using namespace std;
using namespace elem;

//#include <random>
typedef double R;
typedef float F;

//function declarations
void invgamma_matrix(DistMatrix<R> &W, DistMatrix<R> &sqrtW, const int df);
void repmat(DistMatrix<R> &vector, DistMatrix<R> &repmatrix);
void mvnrnd(DistMatrix<R> &Sigma,DistMatrix<R> &Z,const int nSims, mpi::Comm Comm);
void dotproduct(DistMatrix<R> &A,DistMatrix<R> &B);
void skewedtrnd(DistMatrix<R>& X,DistMatrix<R>& mu,DistMatrix<R>& gamma, DistMatrix<R>& Sigma,const int df, const int nSims,mpi::Comm Comm);
void standard_normal_matrix(DistMatrix<R> &Normal);



int nprocs;
int commRank;

int main( int argc, char* argv[] ){
  Initialize( argc, argv );

  mpi::Comm comm = mpi::COMM_WORLD;
  commRank = mpi::CommRank( comm );
  nprocs = mpi::CommSize(comm );

  MPI_Bcast(&nprocs, 1, MPI_INT, 0, comm);
  //MPI::Comm::Bcast(nprocs , 1,mpi::int , 0);
  
  const int n = Input("--size","NxN size of Matricies",10);
  const int df = Input("--df","Degrees of freedom",5);
  const int nSims = Input("--nSims","Number of Monte Carlo Simulations",10);
  const bool print = Input("--print","print matrices?",false);
  ProcessInput();
  PrintInputReport();

  const Grid grid( comm );
  DistMatrix<R> mu(1, n, grid);
  DistMatrix<R> Sigma(n, n, grid);
  DistMatrix<R> gamma(1, n, grid);
  Uniform( 1, n, mu );
  Uniform( 1, n, gamma);
  standard_normal_matrix(Sigma);
  DistMatrix<R> X(nSims,n,grid);
  DistMatrix<R> SigmaT(n,n,grid);
  Copy(Sigma,SigmaT);
  DistMatrix<R> S(n,n,grid);
  Gemm(NORMAL,TRANSPOSE, 1.0, Sigma, SigmaT, 0.0, S);
  S.Print("Positive Definte Symmetric Sigma");
  skewedtrnd(X,mu,gamma, S, df, nSims,comm);
  Sigma.Empty();
  SigmaT.Empty();
  S.Empty();
  MPI_Barrier(comm);
  X.Print("Skewed MVT RND");

  //#ifndef RELEASE
  //DumpCallStack();
  //#endif
   
  Finalize();
  return 0;
}

// X: mu + (gamma/beta * W) * (Z*sqrt(W))
void skewedtrnd(DistMatrix<R>& X,DistMatrix<R>& mu,DistMatrix<R>& gamma, DistMatrix<R>& Sigma,const int df, const int nSims,mpi::Comm Comm){
		/*
		  
		  nVariables = size(gamma, 2);
		  Z = mvnrnd(zeros(1, nVariables),C,nSim);
		  W = df ./ gamrnd(df./2, 2, nSim, nVariables);
		  result = repmat(mu,nSim,1) + repmat(gamma,nSim,1) .* W  + sqrt(W) .* Z;

        
		*/   

  const Grid grid( Comm );
  const R alpha = 1;		

		//nVariables = size(gamma, 2);
  const int nVariables = gamma.Width();
 
  //Generate W and Sqrt W from InvG ~ (df/2, df/2)
  DistMatrix<R> W(nSims, nVariables,grid);
  DistMatrix<R> sqrtW(nSims, nVariables,grid);
  invgamma_matrix(W,sqrtW,df);
  
  DistMatrix<R> Z(nSims,Sigma.Width(), grid);
  mvnrnd(Sigma,Z,nSims,Comm);
  DistMatrix<R> GAMMA(nSims,nVariables,grid);
  DistMatrix<R> MU(nSims,nVariables,grid);
  repmat(gamma,GAMMA);	      
  repmat(mu,MU);
  
   //repmat(gamma,nSim,1) .* W
   dotproduct(GAMMA,W);
   //sqrt(W) .* Z;
   dotproduct(sqrtW,Z);
   //repmat(gamma,nSim,1) .*W
   Axpy(alpha,W,Z);
   //repmat(mu,nSim,1) + repmat(gamma,nSim,1) .* W  + sqrt(W) .* Z
   Axpy(alpha,MU,Z);

   //free matricies allocated   
   /*
   GAMMA.Empty();
   MU.Empty();
   W.Empty();
   sqrtW.Empty();
   */  
   Copy(Z,X);

   //Z.Empty();
}

//function to generate W and Sqrt W for skewed t mvt rnd
void invgamma_matrix(DistMatrix<R> &W, DistMatrix<R> &sqrtW, const int df){

  boost::random::gamma_distribution<> dist(df/2.0,2.0);  
  boost::random::mt19937 rng(static_cast<boost::uint32_t>(commRank));

  double temp;		       
  const int colShift = W.ColShift(); // first row we own
  const int rowShift = W.RowShift(); // first col we own
  const int colStride = W.ColStride();
  const int rowStride = W.RowStride();
  const int localHeight = W.LocalHeight();
  const int localWidth = W.LocalWidth();
      for( int iLocal=0; iLocal<localHeight; ++iLocal ){
	 for( int jLocal=0; jLocal<localWidth; ++jLocal ){
	   //const int i = colShift + iLocal*colStride;
	   //const int j = rowShift + jLocal*rowStride;
	  temp = (double)df/dist(rng);
	  W.SetLocal( iLocal, jLocal, temp  );
	  sqrtW.SetLocal( iLocal, jLocal, sqrt(temp));
      }
  }
}

  //generate N~(0,1) matrix/vector
void standard_normal_matrix(DistMatrix<R> &Normal){
  boost::random::normal_distribution<> dist(0.0,1.0);  
  boost::random::mt19937 rng(static_cast<boost::uint32_t>(commRank));
  double temp;			       
  const int colShift = Normal.ColShift(); // first row we own
  const int rowShift = Normal.RowShift(); // first col we own
  const int colStride = Normal.ColStride();
  const int rowStride = Normal.RowStride();
  const int localHeight = Normal.LocalHeight();
  const int localWidth = Normal.LocalWidth();
      for( int iLocal=0; iLocal<localHeight; ++iLocal ){
	 for( int jLocal=0; jLocal<localWidth; ++jLocal ){
	   //const int i = colShift + iLocal*colStride;
	   //const int j = rowShift + jLocal*rowStride;
	  temp = dist(rng);
	  Normal.SetLocal( iLocal, jLocal, temp);
      }
  }
}

void dotproduct(DistMatrix<R> &A,DistMatrix<R> &B){
  if(A.Height() != B.Height()){
    //ERROR!
  }
  if(A.Width() != B.Width()){
    //ERROR!
  }
  double temp,temp1;
 
  const int colShift = A.ColShift(); // first row we own
  const int rowShift = A.RowShift(); // first col we own
  const int colStride = A.ColStride();
  const int rowStride = A.RowStride();
  const int localHeight = A.LocalHeight();
  const int localWidth = A.LocalWidth();
    for( int iLocal=0; iLocal<localHeight; ++iLocal ){
      for( int jLocal=0; jLocal<localWidth; ++jLocal ){
		  const int i = colShift + iLocal*colStride;
		  const int j = rowShift + jLocal*rowStride;
      temp = A.GetLocal(i, j);
      temp1 = B.GetLocal(i,j);
      B.SetLocal(i, j, (temp*temp1));

      }
    }
}

void repmat(DistMatrix<R> &vector, DistMatrix<R> &repmatrix){
  //old 1 x d
  //new d x d
  //Must set dimensions of New before being passed to function
  double temp;
  const int colShift = repmatrix.ColShift(); // first row we own
  const int rowShift = repmatrix.RowShift(); // first col we own
  const int colStride = repmatrix.ColStride();
  const int rowStride = repmatrix.RowStride();
  const int localHeight = repmatrix.LocalHeight();
  const int localWidth = repmatrix.LocalWidth();
    for( int iLocal=0; iLocal<localHeight; ++iLocal ){
      for( int jLocal=0; jLocal<localWidth; ++jLocal ){
		  const int i = colShift + iLocal*colStride;
		  const int j = rowShift + jLocal*rowStride;
      temp = vector.GetLocal(0, j);
      repmatrix.SetLocal(i, j, temp);
    }
  }
}

void mvnrnd(DistMatrix<R> &Sigma,DistMatrix<R> &Z, const int nSims, mpi::Comm Comm){
  // WARNING if Cholesky is NOT postive definite the code segfaults!
  
  // const Orientation orientation = NORMAL;
  const R alpha = 1;
  const R beta = 0;
  Grid grid(Comm);
  DistMatrix<R> N(nSims,Sigma.Width(), grid);
  standard_normal_matrix(N);
  Cholesky(UPPER,Sigma);
  Gemm(NORMAL,NORMAL, alpha, Sigma, N, beta, Z);

  //N.Empty();
}
