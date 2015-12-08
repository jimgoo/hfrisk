# hfrisk

This code does parallel VaR and CVaR forecasting with Copula-ARMA-GARCH models, with application to high-frequency risk management (HFRisk).

Parallelization of univariate ARMA-GARCH estimation done via MPI. Testing was on an IBM Blue Gene/P.

The `configure.sh` script will try to install all dependencies. Once this is done,

```bash
cd src
make
```

## Link to patent:

http://appft1.uspto.gov/netacgi/nph-Parser?Sect1=PTO1&Sect2=HITOFF&d=PG01&p=1&u=/netahtml/PTO/srchnum.html&r=1&f=G&l=50&s1=20140214722.PGNR.
