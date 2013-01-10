/*
  Author: Jimmie Goode
  Created: 2012-09-01
*/

#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <assert.h>
#include <string>
#include <sstream>

/* These conflict with boost libs (see <boost/math/constants/constants.hpp>) */
/* #define E		2.71828182845904523536029 */
/* #define Pi		3.14159265358979323846264 */
/* #define Degree	0.01745329251994329576924 */

// The null value of a mariginal distribution's parameters
#define NAN_PAR 666.0

#define stringify( name ) // name

using namespace std;

// Marginal distributions
enum MarginalType {ST = 0, AS = 1, T_GARCH = 2};

// the size of each type of marginals parameter estimates
enum RES_SIZE {NULL_RES_SIZE = 0, ST_RES_SIZE = 7, AS_RES_SIZE = 7,
			   T_GARCH_RES_SIZE = 11};

// the number of parameters for each marginal distribution
enum PAR_SIZE {NULL_PAR_SIZE = 0, ST_PAR_SIZE = 4, AS_PAR_SIZE = 4, T_GARCH_PAR_SIZE = 7};

// Skewed t parameter indices in MPI result
enum IDX_ST {ST_GAMMA = 0, ST_MU = 1, ST_DF = 2, ST_SIGMA = 3};

// Alpha-stable parameter indices in MPI result
enum IDX_AS {AS_ALPHA = 0, AS_BETA = 1, AS_SIGMA = 2, AS_MU = 3};

enum START_TYPE {COLD = 0, HOT = 1};

// for skewed t estimation
enum EST_METHOD {BFGS = 0, NMSIMPLEX = 1};


//////////////////////////////////////////////////////////////////////////


// Fetch the result size for a given marginal type. Includes estimation info like iterations and LLF.
inline RES_SIZE getResultSize(MarginalType marginalType) {
  switch (marginalType) {
  case ST:
	return ST_RES_SIZE;
  case AS:
	return AS_RES_SIZE;
  case T_GARCH:
	return T_GARCH_RES_SIZE;
  default:
	assert(0 && "Unknown marginal type.\n");
	return NULL_RES_SIZE;
  }
}

// fetch the number of actually parameters in the marginal model
inline PAR_SIZE getParamCount(MarginalType marginalType) {
  switch (marginalType) {
  case ST:
	return ST_PAR_SIZE;
  case AS:
	return AS_PAR_SIZE;
  case T_GARCH:
	return T_GARCH_PAR_SIZE;
  default:
	assert(0 && "Unknown marginal type.\n");
	return NULL_PAR_SIZE;
  }
}


#endif
