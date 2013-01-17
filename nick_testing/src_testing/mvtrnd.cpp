/*
  Skewed T MVT RND using elemental
  X:= mu + Gamma*W + Z*sqrt(W)

  Alpha Stable RNG using elemental
  X: mu + Z*W


  <TODO>
  * repmatrix function - needs optimization to use just one DistMatrix instead of two
  * Set seeds for production code to be randomized
  

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
void mvnrnd(DistMatrix<R> &Sigma,DistMatrix<R> &Z,const int nSims,Grid** grid, mpi::Comm Comm);
void dotproduct(DistMatrix<R> &A,DistMatrix<R> &B);
void skewedtrnd(DistMatrix<R> &X,DistMatrix<R> &mu,DistMatrix<R> &gamma, DistMatrix<R> &Sigma,const int df, const int nSims,mpi::Comm Comm);
void standard_normal_matrix(DistMatrix<R> &Normal);

void calc_x_if(const R alpha,const R beta,const R c, const R  delta,DistMatrix<R> &x, DistMatrix<R> &phi,DistMatrix<R> &w);
void calc_x_else(const R alpha, const R beta,const R c, const R delta, DistMatrix<R> &x, DistMatrix<R> &phi,DistMatrix<R> &w);
void  mvSubGaussStablernd(DistMatrix<R> &X,DistMatrix<R> &Sigma, DistMatrix<R> &mu, const R alpha,const int nSims, mpi::Comm Comm);
void stablernd(const int M, const int N,const R alpha, const R beta,const R c,const R delta, DistMatrix<R> &x,Grid** grid, mpi::Comm Comm);



int nprocs;
int commRank;
int verbose= 1;

int main( int argc, char* argv[] ){
  Initialize( argc, argv );

  mpi::Comm comm = mpi::COMM_WORLD;
  commRank = mpi::CommRank( comm );
  nprocs = mpi::CommSize(comm );

  MPI_Bcast(&nprocs, 1, MPI_INT, 0, comm);
  MPI_Bcast(&verbose,1 ,MPI_INT, 0, comm);
  //MPI::Comm::Bcast(nprocs , 1,mpi::int , 0);
  
  const int n = Input("--size","NxN size of Matricies",10);
  const int df = Input("--df","Degrees of freedom",5);
  const int nSims = Input("--nSims","Number of Monte Carlo Simulations",10);
  const R alpha =  Input("--alpha","skewness parameter",1.4);

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
  


  //SKEWED T HERE
  skewedtrnd(X,mu,gamma, S, df, nSims,comm);
 
  MPI_Barrier(comm);
  X.Print("Skewed MVT RND");
 
  
  DistMatrix<R> mu2(1, n, grid);
  DistMatrix<R> Sigma2(n, n, grid);
  Uniform( 1, n, mu );
  standard_normal_matrix(Sigma2);
  DistMatrix<R> X2(nSims,n,grid);
  DistMatrix<R> SigmaT2(n,n,grid);
  Copy(Sigma2,SigmaT2);
  DistMatrix<R> S2(n,n,grid);
  Gemm(NORMAL,TRANSPOSE, 1.0, Sigma2, SigmaT2, 0.0, S2);
  
  
  DistMatrix<R> Xas(nSims,n,grid);

  //ALPHA STABLE SUB GAUSS HERE
  mvSubGaussStablernd(Xas,S2, mu2, alpha,nSims,comm);

  Xas.Print("Alpha Stable RND");
  
  cout<<"Success!"<<endl;

#ifndef RELEASE
  DumpCallStack();
#endif
   
  Finalize();
  return 0;
}



// X: mu + (gamma/beta * W) * (Z*sqrt(W))
void skewedtrnd(DistMatrix<R>& X,DistMatrix<R>& mu,DistMatrix<R>& gamma, DistMatrix<R>& Sigma,const int df, const int nSims,mpi::Comm Comm){

  Grid* grid;
  Grid* GGrid = new Grid(Comm);

  const R alpha = 1;		

  //nVariables = size(gamma, 2);
  const int nVariables = gamma.Width();
  
  //Generate W and Sqrt W from InvG ~ (df/2, df/2)
  DistMatrix<R> W(nSims, nVariables,*GGrid);
  DistMatrix<R> sqrtW(nSims, nVariables,*GGrid);

  invgamma_matrix(W,sqrtW,df);


  DistMatrix<R> Z(nSims,Sigma.Width(), *GGrid);

  mvnrnd(Sigma,Z,nSims,&grid,Comm);
  MPI_Barrier(Comm);

  DistMatrix<R> GAMMA(nSims,nVariables,*GGrid);
  DistMatrix<R> MU(nSims,nVariables,*GGrid);

  repmat(gamma,GAMMA);	      
  repmat(mu,MU);


  //repmat(gamma,nSim,1) .* W
  dotproduct(GAMMA,W);
  MPI_Barrier(Comm);

  //sqrt(W) .* Z;
  dotproduct(sqrtW,Z);
  MPI_Barrier(Comm);

  //repmat(gamma,nSim,1) .*W
  Axpy(alpha,W,Z);
  MPI_Barrier(Comm);

  //repmat(mu,nSim,1) + repmat(gamma,nSim,1) .* W  + sqrt(W) .* Z
  Axpy(alpha,MU,Z);
  MPI_Barrier(Comm);

  //free matricies allocated   
  GAMMA.Empty();
  MU.Empty();
  W.Empty();
  sqrtW.Empty();
  X = Z;
  Z.Empty();

  delete grid;
  delete GGrid;
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
      temp = A.GetLocal(iLocal, jLocal);
      temp1 = B.GetLocal(iLocal,jLocal);
      B.SetLocal(iLocal, jLocal, (temp*temp1));
    }
  }
}

void repmat(DistMatrix<R> &vector, DistMatrix<R> &repmatrix){
  //old 1 x d
  //new nSim x d
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
      temp = vector.GetLocal(0, iLocal);
      repmatrix.SetLocal(iLocal, jLocal, temp);
    }
  }
}

void mvnrnd(DistMatrix<R> &Sigma,DistMatrix<R> &Z, const int nSims,Grid** grid, mpi::Comm Comm){
  // WARNING if Cholesky is NOT postive definite the code segfaults!
  const R alpha = 1;
  const R beta = 0;
  //Grid grid(Comm);
  *grid = new Grid(Comm);
  DistMatrix<R> N(nSims,Sigma.Width(), **grid);
  standard_normal_matrix(N);
  Cholesky(UPPER,Sigma);
  Gemm(NORMAL,NORMAL, alpha, Sigma, N, beta, Z);
  N.Empty();
}


//w = -log(rand(m,n));
void log_uniform(DistMatrix<R> &U){

  //boost::random::gamma_distribution<> dist(df/2.0,2.0);  
  boost::random::uniform_01<> dist;
  boost::random::mt19937 rng(static_cast<boost::uint32_t>(commRank));

  double temp;	

  const int colShift = U.ColShift(); // first row we own
  const int rowShift = U.RowShift(); // first col we own
  const int colStride = U.ColStride();
  const int rowStride = U.RowStride();
  const int localHeight = U.LocalHeight();
  const int localWidth = U.LocalWidth();
  for( int iLocal=0; iLocal<localHeight; ++iLocal ){
    for( int jLocal=0; jLocal<localWidth; ++jLocal ){
      const int i = colShift + iLocal*colStride;
      const int j = rowShift + jLocal*rowStride;
      temp = dist(rng);

      U.SetLocal(iLocal, jLocal, log(temp));
    
    }
  }
}


//phi = (rand(m,n)-.5)*pi;
void phi_uniform(DistMatrix<R> &U){
      
  boost::random::uniform_01<> dist;
  boost::random::mt19937 rng(static_cast<boost::uint32_t>(commRank));

  R temp;	
  R pi = 4*atan(1);

  const int colShift = U.ColShift(); // first row we own
  const int rowShift = U.RowShift(); // first col we own
  const int colStride = U.ColStride();
  const int rowStride = U.RowStride();
  const int localHeight = U.LocalHeight();
  const int localWidth = U.LocalWidth();
  for( int iLocal=0; iLocal<localHeight; ++iLocal ){
    for( int jLocal=0; jLocal<localWidth; ++jLocal ){
      const int i = colShift + iLocal*colStride;
      const int j = rowShift + jLocal*rowStride;
      temp = dist(rng);

      U.SetLocal(iLocal, jLocal, ((temp-0.5)*pi) );
    
    }
  }
}

// stablernd(m, n, alpha, beta, c, delta, type)
//W = stablernd(nSim, 1, alpha/2, 1, cos(pi*alpha/4)^(2/alpha), 0);
void stablernd(const int M, const int N,const R alpha, const R beta,const R c,const R delta, DistMatrix<R> &x,Grid** grid, mpi::Comm Comm){
  
  *grid = new Grid(Comm);

  if (alpha < 0.1 || alpha > 2){
    //ERROR
  }

  if (abs(beta) > 1){
    //ERROR
  }
  
  const R pi = 4*atan(1);


  DistMatrix<R> w(M,N,**grid); 
  DistMatrix<R> phi(M,N,**grid);
  //DistMatrix<R> x(M,N,**grid);

  DistMatrix<R> xx(M,N,**grid);

  //Generate exponential w and uniform phi
  log_uniform(w);
  phi_uniform(phi);
  
  //cosphi = cos(phi);
  if (abs(alpha-1) > 1.0E-8){
  
    calc_x_if(alpha,beta,c,delta,xx,phi,w);
    
  }

  else{   
    calc_x_else(alpha, beta,c,delta,xx, phi,w);
  }
  
  phi.Empty();
  w.Empty();
  
  x = xx;

  xx.Empty();
}
 
void calc_x_if(const R alpha,const R beta,const R c, const R delta,DistMatrix<R> &x, DistMatrix<R> &phi,DistMatrix<R> &w){
 
  R temp; // phi
  R temp2; // w
  R aphi;
  R a1phi;
  const R pi = 4*atan(1);
  const R zeta = beta * tan(pi * alpha/2);

  const int colShift = x.ColShift(); // first row we own
  const int rowShift = x.RowShift(); // first col we own
  const int colStride = x.ColStride();
  const int rowStride = x.RowStride();
  const int localHeight = x.LocalHeight();
  const int localWidth = x.LocalWidth();
  for( int iLocal=0; iLocal<localHeight; ++iLocal ){
    for( int jLocal=0; jLocal<localWidth; ++jLocal ){
      const int i = colShift + iLocal*colStride;
      const int j = rowShift + jLocal*rowStride;
      temp = phi.GetLocal(iLocal,jLocal);//phi
      temp2 = w.GetLocal(iLocal,jLocal); //w
       aphi = alpha * temp;
       a1phi = (1-alpha)*temp;

       x.SetLocal(iLocal, jLocal, (delta + c *( ( (sin(aphi)+zeta * cos(aphi))/ cos(temp)) * -1 * pow( abs( ((cos(a1phi) + zeta * sin(a1phi))/ (temp2 * cos(temp))) ) ,((1-alpha)/alpha) ) + beta * tan(pi * alpha/2) ))  );


    }
  }
}

void calc_x_else(const R alpha, const R beta,const R c, const R delta, DistMatrix<R> &x, DistMatrix<R> &phi,DistMatrix<R> &w){
/*
   bphi = (pi/2) + beta * phi;
        x = (2/pi) * (bphi .* tan(phi) - beta * log((pi/2) * w ... 
            .* cosphi ./ bphi));
        if alpha ~= 1
            x = x + beta * tan(pi * alpha/2);


 */
R bphi;
R temp;
R temp2;
const R pi = 4*atan(1);

  const int colShift = x.ColShift(); // first row we own
  const int rowShift = x.RowShift(); // first col we own
  const int colStride = x.ColStride();
  const int rowStride = x.RowStride();
  const int localHeight = x.LocalHeight();
  const int localWidth = x.LocalWidth();
  for( int iLocal=0; iLocal<localHeight; ++iLocal ){
    for( int jLocal=0; jLocal<localWidth; ++jLocal ){
      const int i = colShift + iLocal*colStride;
      const int j = rowShift + jLocal*rowStride;
      temp = phi.GetLocal(0,jLocal);//phi
      temp2 = w.GetLocal(0,jLocal); //w
      bphi = (pi/2) + beta * temp;
      //cout<< (delta + c * (((2/pi) * (bphi * tan(temp) - beta * log((pi/2) * temp2 * cos(temp) / bphi))) + (beta * tan (pi *alpha/2))))<<endl;
      x.SetLocal(iLocal, jLocal,  (delta + c * (((2/pi) * (bphi * tan(temp) - beta * log((pi/2) * temp2 * cos(temp) / bphi))) + (beta * tan (pi *alpha/2)))));
    }
  }
} 

//result = mvSubGaussStablernd(nSim, alpha, Sigma, mu)
void  mvSubGaussStablernd(DistMatrix<R> &X,DistMatrix<R> &Sigma, DistMatrix<R> &mu, const R alpha,const int nSims, mpi::Comm Comm){

  Grid* grid;//grid for in stablernd
  Grid* GGrid = new Grid(Comm);//grid for here
    
  // nDim = size(Sigma, 2); 
  const int nDim = Sigma.Width();
  const R pi = 4*atan(1);
    
  // Z = mvnrnd(zeros(1, nDim), Sigma, nSim);
  DistMatrix<R> Z(nSims, nDim, *GGrid);
  mvnrnd(Sigma,Z,nSims,&grid,Comm);

  //W = stablernd(nSim, 1, alpha/2, 1, cos(pi*alpha/4)^(2/alpha), 0);
  DistMatrix<R> x(nSims, 1, *GGrid);

  stablernd(nSims, 1, (alpha/2), 1, pow(cos(pi*alpha/4),(2/alpha) ),0, x, &grid, Comm); //sqrt(W) returned
  MPI_Barrier(Comm);
  

  DistMatrix<R> W(nSims, nDim, *GGrid);
  repmat(x,W);

  x.Empty();
  DistMatrix<R> MU(nSims,nDim,*GGrid);
  repmat(mu,MU);
  
  // Z .* repmat(sqrt(W), 1, nDim)
  dotproduct(Z,W); 
  MPI_Barrier(Comm);
  Z.Empty();

  Axpy(1,MU,W);
  MPI_Barrier(Comm);
  MU.Empty();
  
  X = W;
  W.Empty();

  delete grid;
  delete GGrid;
}
