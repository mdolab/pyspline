#ifndef PY_SPLINE_H
#define PY_SPLINE_H

/*
  C/C++ header file for the pySpline fortran routines.
*/

#include "TACSObject.h"

#ifdef __cplusplus
extern "C" {
#endif
  extern double bvalu_( double * TX, double * BCOEF, 
			int * NX, int * KX, 
			int * IDX, double * XVAL, int * INBV,
			double * WORK );  

  extern double b2val_( double * XVAL, double * YVAL, 
			int * IDX, int * IDY, 
			double * TX, double * TY, 
			int * NX, int * NY, 
			int * KX, int * KY, 
			double * BCOEF, double * WORK );

  extern double b3val_( double * XVAL, double * YVAL, double * ZVAL,
			int * IDX, int * IDY, int * IDZ,
			double * TX, double * TY, double * TZ,
			int * NX, int * NY, int * NZ, 
			int * KX, int * KY, int * KZ,
			double * BCOEF, double * WORK );

  extern TacsComplex cbvalu_( double * TX, TacsComplex * BCOEF, 
			      int * NX, int * KX, 
			      int * IDX, double * XVAL, int * INBV,
			      TacsComplex * WORK );  

  extern TacsComplex cb2val_( double * XVAL, double * YVAL, 
			      int * IDX, int * IDY, 
			      double * TX, double * TY, 
			      int * NX, int * NY, 
			      int * KX, int * KY, 
			      TacsComplex * BCOEF, TacsComplex * WORK );
  
  extern TacsComplex cb3val_( double * XVAL, double * YVAL, double * ZVAL,
			      int * IDX, int * IDY, int * IDZ,
			      double * TX, double * TY, double * TZ,
			      int * NX, int * NY, int * NZ, 
			      int * KX, int * KY, int * KZ,
			      TacsComplex * BCOEF, TacsComplex * WORK );
#ifdef __cplusplus
}
#endif

#ifdef TACS_USE_COMPLEX
#define BVALU cbvalu_
#define B2VAL cb2val_
#define B3VAL cb3val_
#else
#define BVALU bvalu_
#define B2VAL b2val_
#define B3VAL b3val_
#endif

#endif
