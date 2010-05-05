#ifndef PY_SPLINE_H
#define PY_SPLINE_H

/*
  C/C++ header file for the pySpline fortran routines.
*/

#ifdef __cplusplus
extern "C" {
#endif

  /* 
     Evaluate the curve 
  */

  extern void eval_curve_( double * u, 
			   double * tu,
			   int * ku,
			   double * coef, 
			   int * nctlu,
			   int * ndim, double * val );
  /* 
     Evaluate a surface 
  */
  
  extern void eval_surface_( double * u, double * v, 
			     double * tu, double * tv, 
			     int * ku, int * kv,
			     double * coef, 
			     int * nctlu, int * nctlv, 
			     int * ndim, double * val );

  /*
    Evaluate a volume 
  */

  extern void eval_volume_( double * u, double * v, double * w,
			    double * tu, double * tv, double * tw, 
			    int * ku, int * kv, int * kw,
			    double * coef, 
			    int * nctlu, int * nctlv, int * nctlw,
			    int * ndim, double * val );

  /*
    Functions from cmlib that are useful
  */
  
  extern void bknot_( double * X, int * N, int * K, double * T );

  extern void intrv_( double * XT, int * LXT, double * X, int * ILO,
		      int * ILEFT, int * MFLAG );

  extern void bspvn_( double * T, int * JHIGH, int * K, int * INDEX, 
		      double * X, int * ILEFT, double * VNIKX, 
		      double * WORK, int * IWORK );

  extern void b2ink_( double * X, int * NX, 
		      double * Y, int * NY,
		      double * FCN, int * LDF,
		      int * KX, int * KY, double * TX, double * TY, 
		      double * BCOEFF, double * WORK, int * IFLAG );

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
#ifdef __cplusplus
}
#endif

#endif
