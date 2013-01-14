#--------------------------------------------------------------
mkdir lib
cd lib


#LIBDIR="/gpfs/home3/j/jgoode/lib"
LIBDIR="/nfs/user01/jimmie21/github/hfrisk/lib"

#--------------------------------------------------------------
# Download all libs

wget http://sourceforge.net/projects/arma/files/armadillo-3.6.1.tar.gz
wget http://ab-initio.mit.edu/nlopt/nlopt-2.3.tar.gz
wget http://sourceforge.net/projects/boost/files/boost/1.51.0/boost_1_51_0.tar.gz/download
wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.15.tar.gz
wget http://www.cmake.org/files/v2.8/cmake-2.8.9.tar.gz
wget ftp://ftp.gnu.org/gnu/gcc/gcc-4.7.2/gcc-4.7.2.tar.gz
wget http://www.hdfgroup.org/ftp/HDF5/current/src/hdf5-1.8.9.tar.gz
wget http://www.fftw.org/fftw-3.3.2.tar.gz

tar -xzf armadillo-3.6.1.tar.gz 
tar -xzf boost_1_51_0.tar.gz 
tar -xzf gsl-1.15.tar.gz 
tar -xzf nlopt-2.3.tar.gz
tar -xzf cmake-2.8.9.tar.gz
tar -xzf gcc-4.7.2.tar.gz
tar -xzf hdf5-1.8.9.tar.gz
tar -zxf fftw-3.3.2.tar.gz

#--------------------------------------------------------------
# GCC

# cd gcc-4.7.2/

# # get the old configure options
# gcc -v

# /gpfs/home3/j/jgoode/lib/gcc-4.7.2/configure --enable-threads=posix --prefix=/usr --with-local-prefix=/usr/local --infodir=/usr/share/info --mandir=/usr/share/man --enable-languages=c,c++,f77,objc,java,ada --disable-checking --libdir=/usr/lib --enable-libgcj --with-gxx-include-dir=/usr/include/g++ --with-slibdir=/lib --with-system-zlib --enable-shared --enable-__cxa_atexit --host=powerpc-suse-linux --build=powerpc-suse-linux --target=powerpc-suse-linux --enable-targets=powerpc64-suse-linux --enable-biarch


#--------------------------------------------------------------
# GSL

cd gsl-1.15/
mkdir install/gsl
./configure --prefix=${LIBDIR}/gsl-1.15/install/gsl
make
make install

#--------------------------------------------------------------
# NLOpt

cd ../
cd nlopt-2.3
mkdir install
./configure --prefix=${LIBDIR}/nlopt-2.3/install
make
make install

#--------------------------------------------------------------
# CMake (required for Armadillo installation)

cd ../
cd cmake-2.8.9/
mkdir install
./bootstrap --prefix=${LIBDIR}/cmake-2.8.9/install
make
make install
export PATH=${LIBDIR}/cmake-2.8.9/install/bin:$PATH

#--------------------------------------------------------------
# Boost (<TODO> Armadillo's config.hpp has a flag to use Boost)

cd ../
cd boost_1_51_0/
mkdir install
./bootstrap.sh
./b2 --prefix=${LIBDIR}/boost_1_51_0/install


#--------------------------------------------------------------
# FFTW

cd ../
cd fftw-3.3.2
mkdir installdir
#./configure --prefix=/Users/jimmiegoode/Documents/Glimm/lib/fftw-3.3.2/installdir
#./configure --prefix=/nfs/user03/copula/20120323/lib/fftw-3.3.2/installdir
./configure --prefix=${LIBDIR}/fftw-3.3.2/installdir
make
make install

#--------------------------------------------------------------
# HDF5 is used by Armadillo

cd ../
cd hdf5-1.8.9/
./configure --prefix=${LIBDIR}/hdf5-1.8.9/install --enable-fortran --enable-cxx


#--------------------------------------------------------------
# Armadillo (one needy bitch, got to do it by hand)

cd ../
cd armadillo-3.6.1/
mkdir install
#cmake -DCMAKE_INSTALL_PREFIX:PATH=${LIBDIR}/armadillo-3.6.1/install .
cmake . 
make
make install DESTDIR=${LIBDIR}/armadillo-3.6.1/install

# Now open:
#	armadillo-3.6.1/install/armadillo_bits/config.hpp
# Comment out these:
#	define ARMA_USE_WRAPPER
#	define ARMA_ATLAS_INCLUDE_DIR
# Uncomment these:
#	define ARMA_USE_HDF5
#	define ARMA_PRINT_HDF5_ERRORS
	

#--------------------------------------------------------------
# Galaxy needs the Armadillo lib dir to be included in the LD_LIBRARY_PATH variable.
# If it's not, then mpirun will fail. Galaxy also needs the PATH set as follows.

#export LD_LIBRARY_PATH=/usr/local/pkg/openmpi/lib:/nfs/user03/copula/20120323/lib/armadillo-3.4.1/install/usr/lib64:$LD_LIBRARY_PATH  

#export PATH=/usr/local/pkg/openmpi/bin:$PATH

#--------------------------------------------------------------
# Set number of threads for Intel ESSL BLAS/LAPACK on BlueGene
#OMP_NUM_THREADS

#--------------------------------------------------------------


#--------------------------------------------------------------
# Elemental - new version of CMake is required

# Galaxy
/nfs/user03/copula/20120323/lib/cmake-2.8.9/install/bin/cmake -D CMAKE_INSTALL_PREFIX=????  ..

# MAC
cmake -D CMAKE_INSTALL_PREFIX=/Users/jimmiegoode/Documents/Glimm/lib/elemental/elemental-0.77/build/ \
           -D CMAKE_CXX_COMPILER=/usr/bin/g++ \
           -D CMAKE_C_COMPILER=/usr/bin/gcc   \
           -D CMAKE_Fortran_COMPILER=/opt/local/bin/gfortran \
           -D MATH_LIBS="-D__ACCELERATE_ -framework Accelerate" ..
make
make install
