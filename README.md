# hfrisk

This code does parallel VaR and CVaR forecasting with Copula-ARMA-GARCH models, with applications to high-frequency market risk management (HFRisk). It was written for a [patent](http://appft1.uspto.gov/netacgi/nph-Parser?Sect1=PTO1&Sect2=HITOFF&d=PG01&p=1&u=/netahtml/PTO/srchnum.html&r=1&f=G&l=50&s1=20140214722.PGNR) that I worked on with James Glimm during my PhD.

Parallelization of univariate ARMA-GARCH estimation done via MPI. Testing was on an IBM Blue Gene/P.

The `configure.sh` script will try to install all dependencies in a folder called `lib`. 

The code has been tested with the following:

- FFTW 3.3.2: for computing density functions of CTS and NTS tempered stable and alpha stable distributions. 
- Armadillo 3.6.1: for easy linear algebra
- nlopt 2.3: for nonlinear optimization of ARMA-GARCH log-likelihood functions.
- GSL 1.5: for stats and optimization algorithms.
- Boost 1.52.0
- cmake 2.8.9
- gcc 4.7.2
- HDF5 1.8.9

Additionally, you'll need to link to your existing BLAS, LAPACK, and MPI libraries. You can build with:

```bash
cd src
make
```

