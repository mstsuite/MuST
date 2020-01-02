//-*-C++-*-
// ****************************************************************************
// * C++ wrapper for BLAS                                                     *
// *                                                                          *
// * Thomas Schulthess, ORNL, October 1999                                    *
// * Richard Thigpen, ORNL, June 2003                                         *
// ****************************************************************************


#ifndef PSIMAG_BLAS
#define PSIMAG_BLAS

#include <complex>

/** \file BLAS.h
 *  \author Thomas C. Schulthess and Richard N. Thigpen
 */


/** \brief Namespace for psimag wrappers of BLAS functions
 */
namespace BLAS {

// ============================================================================
// = Level 3 BLAS             GEMM
// ============================================================================
extern "C" void sgemm_(const char*,const char*,int*,int*,int*,const float*,
		       const float*,int*,const float*,int*,
		       const float*,float*,int*);
extern "C" void dgemm_(const char*,const char*,int*,int*,int*,const double*,
		       const double*,int*,const double*,int*,
		       const double*,double*,int*);
extern "C" void cgemm_(const char*,const char*,int*,int*,int*,const std::complex<float>*,
		       const std::complex<float>*,int*,const std::complex<float>*,int*,
		       const std::complex<float>*,std::complex<float>*,int*);
extern "C" void zgemm_(const char*,const char*,int*,int*,int*,const std::complex<double>*,
		       const std::complex<double>*,int*,const std::complex<double>*,int*,
		       const std::complex<double>*,std::complex<double>*,int*);
//*****************************************************************************
//*                           SYMM 
//*****************************************************************************

extern "C" void ssymm_(const char*,const char*,int*,int*,const float*,const float*,int*,
		       const float*,int*,const float*,const float*,int*);
extern "C" void dsymm_(const char*,const char*,int*,int*,const double*,const double*,int*,
		       const double*,int*,const double*,const double*,int*);
extern "C" void csymm_(const char*,const char*,int*,int*,const std::complex<float>*,
		       const std::complex<float>*,int*,
		       const std::complex<float>*,int*,const std::complex<float>*,
		       const std::complex<float>*,int*);
extern "C" void zsymm_(const char*,const char*,int*,int*,const std::complex<double>*,
		       const std::complex<double>*,int*,
		       const std::complex<double>*,int*,const std::complex<double>*,
		       const std::complex<double>*,int*);

//*****************************************************************************
//*                           HEMM
//*****************************************************************************

extern "C" void chemm_(const char*,const char*,int*,int*,const std::complex<float>*,
		       const std::complex<float>*,int*,const std::complex<float>*,
		       int*,const std::complex<float>*,std::complex<float>*,int*);
extern "C" void zhemm_(const char*,const char*,int*,int*,const std::complex<double>*,
		       const std::complex<double>*,int*,const std::complex<double>*,
		       int*,const std::complex<double>*,std::complex<double>*,int*);

// ****************************************************************************
// *                          SYRK
// ****************************************************************************

extern "C" void ssyrk_(const char*,const char*,int*,int*,const float*,const float*,int*,
		       const float*,float*,int*);
extern "C" void dsyrk_(const char*,const char*,int*,int*,const double*,const double*,int*,
		       const double*,double*,int*);
extern "C" void csyrk_(const char*,const char*,int*,int*,const std::complex<float>*,
		       const std::complex<float>*,int*,
		       const std::complex<float>*,std::complex<float>*,int*);
extern "C" void zsyrk_(const char*,const char*,int*,int*,const std::complex<double>*,
		       const std::complex<double>*,int*,
		       const std::complex<double>*,std::complex<double>*,int*);

// ****************************************************************************
// *                          HERK
// ****************************************************************************
extern "C" void cherk_(const char*,const char*,int*,int*,const std::complex<float>*,
		       const std::complex<float>*,int*,
		       const std::complex<float>*,std::complex<float>*,int*);

extern "C" void zherk_(const char*,const char*,int*,int*,const std::complex<double>*,
		       const std::complex<double>*,int*,
		       const std::complex<double>*,std::complex<double>*,int*);
// ****************************************************************************
// *                          SYR2K
// ****************************************************************************
extern "C" void ssyr2k_(const char*,const char*,int*,int*,const float*,const float*,int*,
		       const float*,int*,const float*,float*,int*);
extern "C" void dsyr2k_(const char*,const char*,int*,int*,const double*,const double*,int*,
		       const double*,int*,const double*,double*,int*);
extern "C" void csyr2k_(const char*,const char*,int*,int*,const std::complex<float>*,const std::complex<float>*,
		       int*,const std::complex<float>*,int*,const std::complex<float>*,
		       std::complex<float>*,int*);
extern "C" void zsyr2k_(const char*,const char*,int*,int*,const std::complex<double>*,const std::complex<double>*,
		       int*,const std::complex<double>*,int*,const std::complex<double>*,
		       std::complex<double>*,int*);
// ****************************************************************************
// *                          HER2k
// ****************************************************************************
extern "C" void cher2k_(const char*,const char*,int*,int*,const std::complex<float>*,const std::complex<float>*,
		       int*,const std::complex<float>*,int*,const std::complex<float>*,
		       std::complex<float>*,int*);
extern "C" void zher2k_(const char*,const char*,int*,int*,const std::complex<double>*,const std::complex<double>*,
		       int*,const std::complex<double>*,int*,const std::complex<double>*,
		       std::complex<double>*,int*);
// ****************************************************************************
// *                          TRMM
// ****************************************************************************
extern "C" void strmm_(const char*,const char*,char*,char*,int*,int*,const float*,const float*,
		       int*,float*,int*);
extern "C" void dtrmm_(const char*,const char*,char*,char*,int*,int*,const double*,const double*,
		       int*,double*,int*);
extern "C" void ctrmm_(const char*,const char*,char*,char*,int*,int*,const std::complex<float>*,
		       const std::complex<float>*,int*,std::complex<float>*,int*);
extern "C" void ztrmm_(const char*,const char*,char*,char*,int*,int*,const std::complex<double>*,
		       const std::complex<double>*,int*,std::complex<double>*,int*);
// ****************************************************************************
// *                          TRSM
// ****************************************************************************
extern "C" void strsm_(const char*,const char*,char*,char*,int*,int*,const float*,const float*,
		       int*,float*,int*);
extern "C" void dtrsm_(const char*,const char*,char*,char*,int*,int*,const double*,const double*,
		       int*,double*,int*);
extern "C" void ctrsm_(const char*,const char*,char*,char*,int*,int*,const std::complex<float>*,
		       const std::complex<float>*,int*,std::complex<float>*,int*);
extern "C" void ztrsm_(const char*,const char*,char*,char*,int*,int*,const std::complex<double>*,
		       const std::complex<double>*,int*,std::complex<double>*,int*);
// ****************************************************************************
// *    Level 2 BLAS          GEMV
// ****************************************************************************
extern "C" void sgemv_(const char*,int*,int*,const float*,const float*,int*,
		       const float*,int*,const float*,float*,int*);
extern "C" void dgemv_(const char*,int*,int*,const double*,const double*,int*,
		       const double*,int*,const double*,double*,int*);
extern "C" void cgemv_(const char*,int*,int*,const std::complex<float>*,const std::complex<float>*,
		       int*,const std::complex<float>*,int*,const std::complex<float>*,
		       std::complex<float>*,int*);  
extern "C" void zgemv_(const char*,int*,int*,const std::complex<double>*,const std::complex<double>*,
		       int*,const std::complex<double>*,int*,const std::complex<double>*,
		       std::complex<double>*,int*);
// ****************************************************************************
// *                          GBMV
// ****************************************************************************
extern "C" void sgbmv_(const char*,int*,int*,int*,int*,const float*,const float*,int*,
		       const float*,int*,const float*,float*,int*);
extern "C" void dgbmv_(const char*,int*,int*,int*,int*,const double*,const double*,int*,
		       const double*,int*,const double*,double*,int*);
extern "C" void cgbmv_(const char*,int*,int*,int*,int*,const std::complex<float>*,
		       const std::complex<float>*,int*,const std::complex<float>*,int*,
		       const std::complex<float>*,std::complex<float>*,int*);
extern "C" void zgbmv_(const char*,int*,int*,int*,int*,const std::complex<double>*,
		       const std::complex<double>*,int*,const std::complex<double>*,int*,
		       const std::complex<double>*,std::complex<double>*,int*);
			 
// ****************************************************************************
// *                          HEMV
// ****************************************************************************
extern "C" void chemv_(const char*,int*,const std::complex<float>*,const std::complex<float>*,
		       int*,const std::complex<float>*,int*,const std::complex<float>*,
		       std::complex<float>*,int*);
extern "C" void zhemv_(const char*,int*,const std::complex<double>*,const std::complex<double>*,
		       int*,const std::complex<double>*,int*,const std::complex<double>*,
		       std::complex<double>*,int*);			 
// ****************************************************************************
// *                         HBMV
// ****************************************************************************
extern "C" void chbmv_(const char*,int*,int*,const std::complex<float>*,const std::complex<float>*,
		       int*,const std::complex<float>*,int*,const std::complex<float>*,
		       std::complex<float>*,int*);
extern "C" void zhbmv_(const char*,int*,int*,const std::complex<double>*,const std::complex<double>*,
		       int*,const std::complex<double>*,int*,const std::complex<double>*,
		       std::complex<double>*,int*);			 
// ****************************************************************************
// *                         HPMV
// ****************************************************************************
extern "C" void chpmv_(const char*,int*,const std::complex<float>*,const std::complex<float>*,
		       const std::complex<float>*,int*,const std::complex<float>*,
		       std::complex<float>*,int*);
extern "C" void zhpmv_(const char*,int*,const std::complex<double>*,const std::complex<double>*,
		       const std::complex<double>*,int*,const std::complex<double>*,
		       std::complex<double>*,int*);
// ******************************************************************************
// *                         SYMV
// ******************************************************************************
extern "C" void ssymv_(const char*,int*,const float*,const float*,int*,const float*,int*,const float*,
		       float*,int*);
extern "C" void dsymv_(const char*,int*,const double*,const double*,int*,const double*,int*,
		       const double*,double*,int*);
// ******************************************************************************
// *                         SBMV
// ******************************************************************************
extern "C" void ssbmv_(const char*,int*,int*,const float*,const float*,int*,const float*,
		       int*,const float*,float*,int*);
extern "C" void dsbmv_(const char*,int*,int*,const double*,const double*,int*,const double*, int*,
		       const double*,double*,int*);
// ******************************************************************************
// *                         SPMV
// ******************************************************************************
extern "C" void sspmv_(const char*,int*,const float*,const float*,const float*,int*,const float*,
		       float*,int*);
  extern "C" void dspmv_(const char*,int*,const double*,const double*,const double*,int*,const double*,
		       double*,int*);
// ******************************************************************************
// *                         TRMV
// ******************************************************************************
extern "C" void strmv_(const char*,const char*,const char*,int*,const float*,int*,float*,int*);
extern "C" void dtrmv_(const char*,const char*,const char*,int*,const double*,int*,double*,int*);
extern "C" void ctrmv_(const char*,const char*,const char*,int*,const std::complex<float>*,int*,
		       std::complex<float>*,int*);
extern "C" void ztrmv_(const char*,const char*,const char*,int*,const std::complex<double>*,int*,
		       std::complex<double>*,int*);

// ******************************************************************************
// *                         TBMV
// ******************************************************************************
extern "C" void stbmv_(const char*,const char*,const char*,int*,int*,const float*,int*,float*,int*);
extern "C" void dtbmv_(const char*,const char*,const char*,int*,int*,const double*,int*,double*,int*);
extern "C" void ctbmv_(const char*,const char*,const char*,int*,int*,const std::complex<float>*,int*,
		       std::complex<float>*,int*);
extern "C" void ztbmv_(const char*,const char*,const char*,int*,int*,const std::complex<double>*,int*,
		       std::complex<double>*,int*);
// ******************************************************************************
// *                         TPMV
// ******************************************************************************
extern "C" void stpmv_(const char*,const char*,const char*,int*,const float*,float*,int*);
extern "C" void dtpmv_(const char*,const char*,const char*,int*,const double*,double*,int*);
extern "C" void ctpmv_(const char*,const char*,const char*,int*,const std::complex<float>*,
		       std::complex<float>*,int*);
extern "C" void ztpmv_(const char*,const char*,const char*,int*,const std::complex<double>*,
		       std::complex<double>*,int*);
// ******************************************************************************
// *                         TRSV
// ******************************************************************************
extern "C" void strsv_(const char*,const char*,const char*,int*,const float*,int*,float*,int*);
extern "C" void dtrsv_(const char*,const char*,const char*,int*,const double*,int*,double*,int*);
extern "C" void ctrsv_(const char*,const char*,const char*,int*,const std::complex<float>*,int*,
		       std::complex<float>*,int*);
extern "C" void ztrsv_(const char*,const char*,const char*,int*,const std::complex<double>*,int*,
		       std::complex<double>*,int*);
// ******************************************************************************
// *                         TBSV
// ******************************************************************************
extern "C" void stbsv_(const char*,const char*,const char*,int*,int*,const float*,int*,float*,int*);
extern "C" void dtbsv_(const char*,const char*,const char*,int*,int*,const double*,int*,double*,int*);
extern "C" void ctbsv_(const char*,const char*,const char*,int*,int*,const std::complex<float>*,int*,
		       std::complex<float>*,int*);
extern "C" void ztbsv_(const char*,const char*,const char*,int*,int*,const std::complex<double>*,int*,
		       std::complex<double>*,int*);
// ******************************************************************************
// *                         TPSV
// ******************************************************************************
extern "C" void stpsv_(const char*,const char*,const char*,int*,const float*,float*,int*);
extern "C" void dtpsv_(const char*,const char*,const char*,int*,const double*,double*,int*);
extern "C" void ctpsv_(const char*,const char*,const char*,int*,const std::complex<float>*,
		       std::complex<float>*,int*);
extern "C" void ztpsv_(const char*,const char*,const char*,int*,const std::complex<double>*,
		       std::complex<double>*,int*);
// ******************************************************************************
// *                         GER
// ******************************************************************************
extern "C" void sger_(int*,int*,const float*,const float*,int*,const float*,int*,
		       float*,int*);
extern "C" void dger_(int*,int*,const double*,const double*,int*,const double*,int*,
		       double*,int*);
// ******************************************************************************
// *                         GERU
// ******************************************************************************
extern "C" void cgeru_(int*,int*,const std::complex<float>*,const std::complex<float>*,
		       int*,const std::complex<float>*,int*,std::complex<float>*,int*);
extern "C" void zgeru_(int*,int*,const std::complex<double>*,const std::complex<double>*,
		       int*,const std::complex<double>*,int*,std::complex<double>*,int*);
// ******************************************************************************
// *                         GERC
// ******************************************************************************
extern "C" void cgerc_(int*,int*,const std::complex<float>*,const std::complex<float>*,
		       int*,const std::complex<float>*,int*,std::complex<float>*,int*);
extern "C" void zgerc_(int*,int*,const std::complex<double>*,const std::complex<double>*,
		       int*,const std::complex<double>*,int*,std::complex<double>*,int*);
// ******************************************************************************
// *                         HER
// ******************************************************************************
extern "C" void cher_(const char*,int*,const std::complex<float>*,const std::complex<float>*,
		      int*,std::complex<float>*,int*);
extern "C" void zher_(const char*,int*,const std::complex<double>*,const std::complex<double>*,
		      int*,std::complex<double>*,int*);
// ******************************************************************************
// *                         HPR
// ******************************************************************************
extern "C" void chpr_(const char*,int*,const std::complex<float>*,const std::complex<float>*,
		      int*,std::complex<float>*);
extern "C" void zhpr_(const char*,int*,const std::complex<double>*,const std::complex<double>*,
		      int*,std::complex<float>*);
// ******************************************************************************
// *                         HER2
// ******************************************************************************
extern "C" void cher2_(const char*,int*,const std::complex<float>*,const std::complex<float>*,
		       int*,const std::complex<float>*,int*,std::complex<float>*,int*);
extern "C" void zher2_(const char*,int*,const std::complex<double>*,const std::complex<double>*,
		       int*,const std::complex<double>*,int*,std::complex<double>*,int*);
// ******************************************************************************
// *                         HPR2
// ******************************************************************************
extern "C" void chpr2_(const char*,int*,const std::complex<float>*,const std::complex<float>*,
		       int*,const std::complex<float>*,int*,std::complex<float>*);
extern "C" void zhpr2_(const char*,int*,const std::complex<double>*,const std::complex<double>*,
		       int*,const std::complex<double>*,int*,std::complex<double>*);
// ******************************************************************************
// *                         SYR
// ******************************************************************************
  extern "C" void ssyr_(const char*,int*,const float*,const float*,int*,float*,int*);
  extern "C" void dsyr_(const char*,int*,const double*,const double*,int*,double*,int*);
// ******************************************************************************
// *                         SPR
// ******************************************************************************
  extern "C" void sspr_(const char*,int*,const float*,const float*,int*,float*);
  extern "C" void dspr_(const char*,int*,const double*,const double*,int*,double*);
// ******************************************************************************
// *                         SYR2
// ******************************************************************************
  extern "C" void ssyr2_(const char*,int*,const float*,const float*,int*,const float*,
			 int*,float*,int*);
  extern "C" void dsyr2_(const char*,int*,const double*,const double*,int*,const double*,
			 int*,double*,int*);
// ******************************************************************************
// *                         SPR2
// ******************************************************************************
  extern "C" void sspr2_(const char*,int*,const float*,const float*,int*,const float*,
			 int*,float*);
  extern "C" void dspr2_(const char*,int*,const double*,const double*,int*,const double*,
			 int*,double*);
// ******************************************************************************
// *Level 1 BLAS
// ******************************************************************************
  
extern "C" void srotg_(float*,float*,float*,float*);
extern "C" void drotg_(double*,double*,double*,double*);
  
extern "C" void srotmg_(float*,float*,float*,float*,float*);
extern "C" void drotmg_(double*,double*,double*,double*,double*);
 
extern "C" void srot_(int*,float*,int*,float*,int*,const float*,const float*);
extern "C" void drot_(int*,double*,int*,double*,int*,const double*,const double*);
  
extern "C" void srotm_(int*,float*,int*,float*,int*,const float*);
extern "C" void drotm_(int*,double*,int*,double*,int*,const double*);

  extern "C" void sswap_(int*,float*,int*,float*,int*);
  extern "C" void dswap_(int*,double*,int*,double*,int*);
  extern "C" void cswap_(int*,std::complex<float>*,int*,std::complex<float>*,int*);
  extern "C" void zswap_(int*,std::complex<double>*,int*,std::complex<double>*,int*);


extern "C" void saxpy_(int*,const float*,const float*,int*,float*,int*);
extern "C" void daxpy_(int*,const double*,const double*,int*,double*,int*);
extern "C" void caxpy_(int*,const std::complex<float>*,const std::complex<float>*,int*,
		       std::complex<float>*,int*);
extern "C" void zaxpy_(int*,const std::complex<double>*,const std::complex<double>*,int*,
		       std::complex<double>*,int*);

extern "C" void scopy_(int*,const float*,int*,float*,int*);
extern "C" void dcopy_(int*,const double*,int*,double*,int*);
extern "C" void ccopy_(int*,const std::complex<float>*,int*,std::complex<float>*,int*);
extern "C" void zcopy_(int*,const std::complex<double>*,int*,std::complex<double>*,int*);

extern "C" void sscal_(int*,const float*,float*,int*);
extern "C" void dscal_(int*,const double*,double*,int*);
extern "C" void cscal_(int*,const std::complex<float>*,std::complex<float>*,int*);
extern "C" void zscal_(int*,const std::complex<double>*,std::complex<double>*,int*);

// ============================================================================
inline void GEMM(char c1,char c2,int sX,int sY,int sZ,const float &a,
		 const float* x,int sx,
		 const float* y,int sy,
		 const float&b,float* z,int sz) {
  sgemm_(&c1,&c2,&sX,&sY,&sZ,&a,x,&sx,y,&sy,&b,z,&sz);
}

inline void GEMM(char c1,char c2,int sX,int sY,int sZ,const double &a,
		 const double* x,int sx,
		 const double* y,int sy,
		 const double&b,double* z,int sz) {
  dgemm_(&c1,&c2,&sX,&sY,&sZ,&a,x,&sx,y,&sy,&b,z,&sz);
}

inline void GEMM(char c1,char c2,int sX,int sY,int sZ,const std::complex<float> &a,
		 const std::complex<float>* x,int sx,
		 const std::complex<float>* y,int sy,
		 const std::complex<float>&b,std::complex<float>* z,int sz) {
  cgemm_(&c1,&c2,&sX,&sY,&sZ,&a,x,&sx,y,&sy,&b,z,&sz);
}

inline void GEMM(char c1,char c2,int sX,int sY,int sZ,const std::complex<double> &a,
		 const std::complex<double>* x,int sx,
		 const std::complex<double>* y,int sy,
		 const std::complex<double>&b,std::complex<double>* z,int sz) {
  zgemm_(&c1,&c2,&sX,&sY,&sZ,&a,x,&sx,y,&sy,&b,z,&sz);
}

// ***************************************************************************
inline void SYMM(char c1,char c2,int sX,int sY,const float &a,const float* x,
		 int sx,const float* y,int sy, const float &b,float* z,int sz){
  ssymm_(&c1,&c2,&sX,&sY,&a,x,&sx,y,&sy,&b,z,&sz); 
}

inline void SYMM(char c1,char c2,int sX,int sY,const double &a,const double* x,
		 int sx,const double* y,int sy, const double &b,double* z,int sz){
  dsymm_(&c1,&c2,&sX,&sY,&a,x,&sx,y,&sy,&b,z,&sz);
}
inline void SYMM(char c1,char c2,int sX,int sY,const std::complex<float> &a,
		 const std::complex<float>* x,int sx,
		 const std::complex<float>* y,int sy, 
		 const std::complex<float> &b,std::complex<float>* z,int sz){
  csymm_(&c1,&c2,&sX,&sY,&a,x,&sx,y,&sy,&b,z,&sz);
}
inline void SYMM(char c1,char c2,int sX,int sY,const std::complex<double> &a,
		 const std::complex<double>* x,int sx,
		 const std::complex<double>* y,int sy, 
		 const std::complex<double> &b,std::complex<double>* z,int sz){
  zsymm_(&c1,&c2,&sX,&sY,&a,x,&sx,y,&sy,&b,z,&sz);
}
// ---------------------------------------------------------------------------
inline void HEMM(char c1,char c2,int sX,int sY, const std::complex<float> &a,
		 const std::complex<float>* x,int sx,
		 const std::complex<float>* y,int sy,
		 const std::complex<float> &b, std::complex<float>* z, int sz){
  chemm_(&c1,&c2,&sX,&sY,&a,x,&sx,y,&sx,&b,z,&sz);
}
inline void HEMM(char c1,char c2,int sX,int sY,const std::complex<double> &a,
		 const std::complex<double>* x,int sx,
		 const std::complex<double>* y,int sy,
		 const std::complex<double> &b,std::complex<double>* z,int sz){
  zhemm_(&c1,&c2,&sX,&sY,&a,x,&sx,y,&sx,&b,z,&sz);
}
// **************************************************************************
inline void SYRK(char UPLO,char TRANS,int N,int K,const float &ALPHA,
		 const float* A,int LDA,const float &BETA,float* C,int LDC){
  ssyrk_(&UPLO,&TRANS,&N,&K,&ALPHA,A,&LDA,&BETA,C,&LDC);
}
inline void SYRK(char UPLO,char TRANS,int N,int K,const double &ALPHA,
		 const double* A,int LDA,const double &BETA,double* C,int LDC){
  dsyrk_(&UPLO,&TRANS,&N,&K,&ALPHA,A,&LDA,&BETA,C,&LDC);
}
inline void SYRK(char UPLO,char TRANS,int N,int K,const std::complex<float> &ALPHA,
		 const std::complex<float>* A,int LDA,const std::complex<float> &BETA,
		 std::complex<float>* C,int LDC){
  csyrk_(&UPLO,&TRANS,&N,&K,&ALPHA,A,&LDA,&BETA,C,&LDC);
}
inline void SYRK(char UPLO,char TRANS,int N,int K,const std::complex<double> &ALPHA,
		 const std::complex<double>* A,int LDA,const std::complex<double> &BETA,
		 std::complex<double>* C,int LDC){
  zsyrk_(&UPLO,&TRANS,&N,&K,&ALPHA,A,&LDA,&BETA,C,&LDC);
}

// ***************************************************************************
inline void HERK(char UPLO,char TRANS,int N,int K,const std::complex<float> &ALPHA,
		 const std::complex<float>* A,int LDA,const std::complex<float> &BETA,
		 std::complex<float>* C,int LDC){
  cherk_(&UPLO,&TRANS,&N,&K,&ALPHA,A,&LDA,&BETA,C,&LDC);
}
inline void HERK(char UPLO,char TRANS,int N,int K,const std::complex<double> &ALPHA,
		 const std::complex<double>* A,int LDA,const std::complex<double> &BETA,
		 std::complex<double>* C,int LDC){
  zherk_(&UPLO,&TRANS,&N,&K,&ALPHA,A,&LDA,&BETA,C,&LDC);
}
// ***************************************************************************
inline void SYR2K(char uplo,char trans,int n, int k,const float &alpha,
		  const float* A,int lda,const float* B,int ldb,const float &beta,
		  float* C,int ldc){
  ssyr2k_(&uplo,&trans,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}
inline void SYR2K(char uplo,char trans,int n,int k,const double &alpha,
		  const double* A,int lda,const double* B,int ldb,const double &beta,
		  double* C,int ldc){
  dsyr2k_(&uplo,&trans,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}
inline void SYR2k(char uplo,char trans,int n,int k,const std::complex<float> &alpha,
		  const std::complex<float>* A,int lda,
		  const std::complex<float>* B,int ldb,
		  const std::complex<float> &beta,std::complex<float>* C,int ldc){
  csyr2k_(&uplo,&trans,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}
inline void SYR2k(char uplo,char trans,int n,int k,const std::complex<double> &alpha,
		  const std::complex<double>* A,int lda,
		  const std::complex<double>* B,int ldb,
		  const std::complex<double> &beta,std::complex<double>* C,int ldc){
  zsyr2k_(&uplo,&trans,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}
// ***************************************************************************
inline void HER2k(char uplo,char trans,int n,int k,const std::complex<float> &alpha,
		  const std::complex<float>* A,int lda,
		  const std::complex<float>* B,int ldb,
		  const std::complex<float> &beta,std::complex<float>* C,int ldc){
  cher2k_(&uplo,&trans,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}
inline void HER2k(char uplo,char trans,int n,int k,const std::complex<double> &alpha,
		  const std::complex<double>* A,int lda,
		  const std::complex<double>* B,int ldb,
		  const std::complex<double> &beta,std::complex<double>* C,int ldc){
  zher2k_(&uplo,&trans,&n,&k,&alpha,A,&lda,B,&ldb,&beta,C,&ldc);
}
// ********************************************************************************
inline void TRMM(char side,char uplo,char transa,char diag,int m,int n,const float &alpha,
		 const float* A,int lda,float* B,int ldb){
  strmm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,A,&lda,B,&ldb);
}
inline void TRMM(char side,char uplo,char transa,char diag,int m,int n,const double &alpha,
		 const double* A,int lda,double* B, int ldb){
  dtrmm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,A,&lda,B,&ldb);
}
inline void TRMM(char side,char uplo,char transa,char diag,int m,int n,
	         const std::complex<float> &alpha,const std::complex<float>* A,int lda,
		 std::complex<float>* B,int ldb){
  ctrmm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,A,&lda,B,&ldb);
}
inline void TRMM(char side,char uplo,char transa,char diag,int m,int n,
		 const std::complex<double> &alpha,const std::complex<double>* A, int lda,
		 std::complex<double>* B, int ldb){
 ztrmm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,A,&lda,B,&ldb);
}
// ********************************************************************************
inline void TRSM(char side,char uplo,char transa,char diag,int m,int n,const float &alpha,
		 const float* A,int lda,float* B,int ldb){
  strsm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,A,&lda,B,&ldb);
}
inline void TRSM(char side,char uplo,char transa,char diag,int m,int n,const double &alpha,
		 const double* A,int lda,double* B, int ldb){
  dtrsm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,A,&lda,B,&ldb);
}
inline void TRSM(char side,char uplo,char transa,char diag,int m,int n,
	         const std::complex<float> &alpha,const std::complex<float>* A,int lda,
		 std::complex<float>* B,int ldb){
  ctrsm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,A,&lda,B,&ldb);
}
inline void TRSM(char side,char uplo,char transa,char diag,int m,int n,
		 const std::complex<double> &alpha,const std::complex<double>* A, int lda,
		 std::complex<double>* B, int ldb){
 ztrsm_(&side,&uplo,&transa,&diag,&m,&n,&alpha,A,&lda,B,&ldb);
}
  // ***************************************************************************

inline void GEMV(char c, int M, int N, 
		 const float &alpha, const float *A, int ldA,
		 const float *x, int incX,
		 const float &beta, float *y, int incY) {
  sgemv_(&c,&M,&N,&alpha,A,&ldA,x,&incX,&beta,y,&incY);
}
// ----------------------------------------------------------------------------
inline void GEMV(char c, int M, int N, 
		 const double &alpha, const double *A, int ldA,
		 const double *x, int incX,
		 const double &beta, double *y, int incY) {
  dgemv_(&c,&M,&N,&alpha,A,&ldA,x,&incX,&beta,y,&incY);
}
// ---------------------------------------------------------------------------
inline void GEMV(char c, int M, int N, 
		 const std::complex<float> &alpha, const std::complex<float> *A, int ldA,
		 const std::complex<float> *x, int incX,
		 const std::complex<float> &beta, std::complex<float> *y, int incY) {
  cgemv_(&c,&M,&N,&alpha,A,&ldA,x,&incX,&beta,y,&incY);
}
// ---------------------------------------------------------------------------
inline void GEMV(char c, int M, int N, 
		 const std::complex<double> &alpha, const std::complex<double> *A, 
		 int ldA,const std::complex<double> *x, int incX,
		 const std::complex<double> &beta, std::complex<double> *y, int incY) {
  zgemv_(&c,&M,&N,&alpha,A,&ldA,x,&incX,&beta,y,&incY);
}
// ----------------------------------------------------------------------------
inline void GBMV(char trans,int m,int n,int kl,int ku,const float &alpha,const float *A,
		 int lda,const float *x,int incx,const float &beta,float *y,int incy){
  sgbmv_(&trans,&m,&n,&kl,&ku,&alpha,A,&lda,x,&incx,&beta,y,&incy);
  }
inline void GBMV(char trans,int m,int n,int kl,int ku,const double &alpha,const double *A,
		 int lda,const double *x,int incx,const double &beta,double *y,int incy){
  dgbmv_(&trans,&m,&n,&kl,&ku,&alpha,A,&lda,x,&incx,&beta,y,&incy);
}
inline void GBMV(char trans,int m,int n,int kl,int ku,const std::complex<float> &alpha,
		 const std::complex<float> *A,int lda,const std::complex<float> *x,int incx,
		 const std::complex<float> &beta,std::complex<float> *y,int incy){
  cgbmv_(&trans,&m,&n,&kl,&ku,&alpha,A,&lda,x,&incx,&beta,y,&incy);
}
inline void GBMV(char trans,int m,int n,int kl,int ku, const std::complex<double> &alpha,
		 const std::complex<double> *A,int lda,const std::complex<double> *x,int incx,
		 const std::complex<double> &beta,std::complex<double> *y,int incy){
  zgbmv_(&trans,&m,&n,&kl,&ku,&alpha,A,&lda,x,&incx,&beta,y,&incy);
}
// ****************************************************************************
inline void HEMV(char uplo,int n,const std::complex<float> &alpha,const std::complex<float> *a,
		 int lda,const std::complex<float> *x,int incx,const std::complex<float> &beta,
		 std::complex<float> *y,int incy){
  chemv_(&uplo,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}
inline void HEMV(char uplo,int n,const std::complex<double> &alpha,const std::complex<double> *a,
		 int lda,const std::complex<double> *x,int incx,const std::complex<double> &beta,
		 std::complex<double> *y,int incy){
  zhemv_(&uplo,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}
// **************************************************************************
inline void HBMV(char uplo,int n,int k,const std::complex<float> &alpha,
		 const std::complex<float> *a,int lda,const std::complex<float> *x,
		 int incx,const std::complex<float> &beta,std::complex<float> *y,int incy){
  chbmv_(&uplo,&n,&k,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}
inline void HBMV(char uplo,int n,int k,const std::complex<double> &alpha,
		 const std::complex<double> *a,int lda,const std::complex<double> *x,
		 int incx,const std::complex<double> &beta,std::complex<double> *y,int incy){
  zhbmv_(&uplo,&n,&k,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}
// ***************************************************************************
inline void HPMV(char uplo,int n,const std::complex<float> &alpha,const std::complex<float> *ap,
		 const std::complex<float> *x,int incx,const std::complex<float> &beta,
		 std::complex<float> *y,int incy){
  chpmv_(&uplo,&n,&alpha,ap,x,&incx,&beta,y,&incy);
}
inline void HPMV(char uplo,int n,const std::complex<double> &alpha,const std::complex<double> *ap,
		 const std::complex<double> *x,int incx,const std::complex<double> &beta,
		 std::complex<double> *y,int incy){
  zhpmv_(&uplo,&n,&alpha,ap,x,&incx,&beta,y,&incy);
}
// ***************************************************************************
inline void SYMV(char uplo,int n,const float &alpha,const float *a,int lda,const float *x,
		 int incx,const float &beta,float *y,int incy){
  ssymv_(&uplo,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}
inline void SYMV(char uplo,int n,const double &alpha,const double *a,int lda,const double *x,
		 int incx,const double &beta,double *y,int incy){
  dsymv_(&uplo,&n,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}
// ****************************************************************************
inline void SBMV(char uplo,int n,int k,const float &alpha,const float *a,int lda,
		 const float *x,int incx,const float &beta, float *y,int incy){
  ssbmv_(&uplo,&n,&k,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}
inline void SBMV(char uplo,int n,int k,const double &alpha,const double *a,int lda,
		 const double *x,int incx,const double &beta,double *y,int incy){
  dsbmv_(&uplo,&n,&k,&alpha,a,&lda,x,&incx,&beta,y,&incy);
}
// ****************************************************************************
inline void SPMV(char uplo,int n,const float &alpha,const float *ap,const float *x,
		 int incx,const float &beta,float *y,int incy){
  sspmv_(&uplo,&n,&alpha,ap,x,&incx,&beta,y,&incy);
}
inline void SPMV(char uplo,int n,const double &alpha,const double *ap,const double *x,
		 int incx,const double &beta,double *y,int incy){
  dspmv_(&uplo,&n,&alpha,ap,x,&incx,&beta,y,&incy);
}
// ****************************************************************************
inline void TRMV(char uplo,char trans,char diag,int n,const float *a,int lda, 
		 float *x,int incx){
  strmv_(&uplo,&trans,&diag,&n,a,&lda,x,&incx);
}
inline void TRMV(char uplo,char trans,char diag,int n,const double *a,int lda, 
		 double *x,int incx){
  dtrmv_(&uplo,&trans,&diag,&n,a,&lda,x,&incx);
}
inline void TRMV(char uplo,char trans,char diag,int n,const std::complex<float> *a,int lda, 
		 std::complex<float> *x,int incx){
  ctrmv_(&uplo,&trans,&diag,&n,a,&lda,x,&incx);
}
inline void TRMV(char uplo,char trans,char diag,int n,const std::complex<double> *a,int lda, 
		 std::complex<double> *x,int incx){
  ztrmv_(&uplo,&trans,&diag,&n,a,&lda,x,&incx);
}
// ****************************************************************************
inline void TBMV(char uplo,char trans,char diag,int n,int k,const float *a,int lda,
		 float *x,int incx){
  stbmv_(&uplo,&trans,&diag,&n,&k,a,&lda,x,&incx);
}
inline void TBMV(char uplo,char trans,char diag,int n,int k,const double *a,int lda,
		 double *x,int incx){
  dtbmv_(&uplo,&trans,&diag,&n,&k,a,&lda,x,&incx);
}
inline void TBMV(char uplo,char trans,char diag,int n,int k,const std::complex<float> *a,
		 int lda,std::complex<float> *x,int incx){
  ctbmv_(&uplo,&trans,&diag,&n,&k,a,&lda,x,&incx);
}
inline void TBMV(char uplo,char trans,char diag,int n,int k,const std::complex<double> *a,
		 int lda,std::complex<double> *x,int incx){
  ztbmv_(&uplo,&trans,&diag,&n,&k,a,&lda,x,&incx);
}
// ****************************************************************************
inline void TPMV(char uplo,char trans,char diag,int n,const float *ap,float *x,int incx){
  stpmv_(&uplo,&trans,&diag,&n,ap,x,&incx);
}
inline void TPMV(char uplo,char trans,char diag,int n,const double *ap,double *x,int incx){
  dtpmv_(&uplo,&trans,&diag,&n,ap,x,&incx);
}
inline void TPMV(char uplo,char trans,char diag,int n,const std::complex<float> *ap,
		 std::complex<float> *x,int incx){
  ctpmv_(&uplo,&trans,&diag,&n,ap,x,&incx);
}
inline void TPMV(char uplo,char trans,char diag,int n,const std::complex<double> *ap,
		 std::complex<double> *x,int incx){
  ztpmv_(&uplo,&trans,&diag,&n,ap,x,&incx);
}
// ****************************************************************************
inline void TRSV(char uplo,char trans,char diag,int n,const float *a,int lda, 
		 float *x,int incx){
  strsv_(&uplo,&trans,&diag,&n,a,&lda,x,&incx);
}
inline void TRSV(char uplo,char trans,char diag,int n,const double *a,int lda, 
		 double *x,int incx){
  dtrsv_(&uplo,&trans,&diag,&n,a,&lda,x,&incx);
}
inline void TRSV(char uplo,char trans,char diag,int n,const std::complex<float> *a,int lda, 
		 std::complex<float> *x,int incx){
  ctrsv_(&uplo,&trans,&diag,&n,a,&lda,x,&incx);
}
inline void TRSV(char uplo,char trans,char diag,int n,const std::complex<double> *a,int lda, 
		 std::complex<double> *x,int incx){
  ztrsv_(&uplo,&trans,&diag,&n,a,&lda,x,&incx);
}
// ****************************************************************************
inline void TBSV(char uplo,char trans,char diag,int n,int k,const float *a,int lda,
		 float *x,int incx){
  stbsv_(&uplo,&trans,&diag,&n,&k,a,&lda,x,&incx);
}
inline void TBSV(char uplo,char trans,char diag,int n,int k,const double *a,int lda,
		 double *x,int incx){
  dtbsv_(&uplo,&trans,&diag,&n,&k,a,&lda,x,&incx);
}
inline void TBSV(char uplo,char trans,char diag,int n,int k,const std::complex<float> *a,
		 int lda,std::complex<float> *x,int incx){
  ctbsv_(&uplo,&trans,&diag,&n,&k,a,&lda,x,&incx);
}
inline void TBSV(char uplo,char trans,char diag,int n,int k,const std::complex<double> *a,
		 int lda,std::complex<double> *x,int incx){
  ztbsv_(&uplo,&trans,&diag,&n,&k,a,&lda,x,&incx);
}
// ****************************************************************************
inline void TPSV(char uplo,char trans,char diag,int n,const float *ap,float *x,int incx){
  stpsv_(&uplo,&trans,&diag,&n,ap,x,&incx);
}
inline void TPSV(char uplo,char trans,char diag,int n,const double *ap,double *x,int incx){
  dtpsv_(&uplo,&trans,&diag,&n,ap,x,&incx);
}
inline void TPSV(char uplo,char trans,char diag,int n,const std::complex<float> *ap,
		 std::complex<float> *x,int incx){
  ctpsv_(&uplo,&trans,&diag,&n,ap,x,&incx);
}
inline void TPSV(char uplo,char trans,char diag,int n,const std::complex<double> *ap,
		 std::complex<double> *x,int incx){
  ztpsv_(&uplo,&trans,&diag,&n,ap,x,&incx);
}
// ****************************************************************************
inline void GER(int m,int n,const float &alpha,const float *x,int incx,const float *y,
		int incy,float *a,int lda){
  sger_(&m,&n,&alpha,x,&incx,y,&incy,a,&lda);
}
inline void GER(int m,int n,const double &alpha,const double *x,int incx,const double *y,
		int incy,double *a,int lda){
  dger_(&m,&n,&alpha,x,&incx,y,&incy,a,&lda);
}
// ****************************************************************************
inline void GERU(int m,int n,const std::complex<float> &alpha,
		 const std::complex<float> *x,int incx,
		 const std::complex<float> *y,int incy,
		 std::complex<float> *a,int lda){
  cgeru_(&m,&n,&alpha,x,&incx,y,&incy,a,&lda);
}
inline void GERU(int m,int n,const std::complex<double> &alpha,
		 const std::complex<double> *x,int incx,
		 const std::complex<double> *y,int incy,
		 std::complex<double> *a,int lda){
  zgeru_(&m,&n,&alpha,x,&incx,y,&incy,a,&lda);
}
// ****************************************************************************
inline void GERC(int m,int n,const std::complex<float> &alpha,
		 const std::complex<float> *x,int incx,
		 const std::complex<float> *y,int incy,
		 std::complex<float> *a,int lda){
  cgerc_(&m,&n,&alpha,x,&incx,y,&incy,a,&lda);
}
inline void GERC(int m,int n,const std::complex<double> &alpha,
		 const std::complex<double> *x,int incx,
		 const std::complex<double> *y,int incy,
		 std::complex<double> *a,int lda){
  zgerc_(&m,&n,&alpha,x,&incx,y,&incy,a,&lda);
}
// *****************************************************************************
inline void HER(char uplo,int n,const std::complex<float> &alpha,
		const std::complex<float> *x,int incx,
		std::complex<float> *a,int lda){
  cher_(&uplo,&n,&alpha,x,&incx,a,&lda);
}
inline void HER(char uplo,int n,const std::complex<double> &alpha,
		const std::complex<double> *x,int incx,
		std::complex<double> *a,int lda){
  zher_(&uplo,&n,&alpha,x,&incx,a,&lda);
}
// *****************************************************************************
inline void HPR(char uplo,int n,const std::complex<float> &alpha,
		const std::complex<float> *x,int incx,
		std::complex<float> *ap){
  chpr_(&uplo,&n,&alpha,x,&incx,ap);
}
inline void HPR(char uplo,int n,const std::complex<double> &alpha,
		const std::complex<double> *x,int incx,
		std::complex<float> *ap){
  zhpr_(&uplo,&n,&alpha,x,&incx,ap);  
}
// *****************************************************************************
inline void HER2(char uplo,int n,const std::complex<float> &alpha,
		 const std::complex<float> *x,int incx,
		 const std::complex<float> *y,int incy,
		 std::complex<float> *a,int lda){
  cher2_(&uplo,&n,&alpha,x,&incx,y,&incy,a,&lda);
}
inline void HER2(char uplo,int n,const std::complex<double> &alpha,
		 const std::complex<double> *x,int incx,
		 const std::complex<double> *y,int incy,
		 std::complex<double> *a,int lda){
  zher2_(&uplo,&n,&alpha,x,&incx,y,&incy,a,&lda);
}
// *****************************************************************************
inline void HPR2(char uplo,int n,const std::complex<float> &alpha,
		 const std::complex<float> *x,int incx,
		 const std::complex<float> *y,int incy,
		 std::complex<float> *ap){
  chpr2_(&uplo,&n,&alpha,x,&incx,y,&incy,ap);
}
inline void HPR2(char uplo,int n,const std::complex<double> &alpha,
		 const std::complex<double> *x,int incx,
		 const std::complex<double> *y,int incy,
		 std::complex<double> *ap){
  zhpr2_(&uplo,&n,&alpha,x,&incx,y,&incy,ap);
}
// *****************************************************************************
  inline void SYR(char uplo,int n,const float &alpha,const float *x,int incx,
		  float *a,int lda){
    ssyr_(&uplo,&n,&alpha,x,&incx,a,&lda);
  }
  inline void SYR(char uplo,int n,const double &alpha,const double *x,int incx,
		  double*a,int lda){
    dsyr_(&uplo,&n,&alpha,x,&incx,a,&lda);
  }
// ****************************************************************************
  inline void SPR(char uplo,int n,const float &alpha,const float *x,int incx,float *ap){
    sspr_(&uplo,&n,&alpha,x,&incx,ap);
  }
  inline void SPR(char uplo,int n,const double &alpha,const double *x,int incx,double *ap){
    dspr_(&uplo,&n,&alpha,x,&incx,ap);
  }
// ****************************************************************************
  inline void SYR2(char uplo,int n,const float &alpha,const float *x,
		   int incx,const float *y,int incy,float *a,int lda){
    ssyr2_(&uplo,&n,&alpha,x,&incx,y,&incy,a,&lda);
  }
  inline void SYR2(char uplo,int n,const double &alpha,const double *x,
		   int incx,const double *y,int incy,double *a,int lda){
    dsyr2_(&uplo,&n,&alpha,x,&incx,y,&incy,a,&lda);
  }
// ****************************************************************************
  inline void SPR2(char uplo,int n,const float &alpha,const float *x,
		   int incx,const float *y,int incy,float *ap){
    sspr2_(&uplo,&n,&alpha,x,&incx,y,&incy,ap);
  }
  inline void SPR2(char uplo,int n,const double &alpha,const double *x,
		   int incx,const double *y,int incy,double *ap){
    dspr2_(&uplo,&n,&alpha,x,&incx,y,&incy,ap);
  }

// ****************************************************************************


inline
void AXPY(int size,const float &a,const float* x, int sx, float* y, int sy){
  saxpy_(&size,&a,x,&sx,y,&sy);
}

inline
void AXPY(int size,const double &a,const double* x, int sx, double* y, int sy){
  daxpy_(&size,&a,x,&sx,y,&sy);
}

inline
void AXPY(int size,const std::complex<float> &a,const std::complex<float>* x, int sx, 
	  std::complex<float>* y, int sy){
  caxpy_(&size,&a,x,&sx,y,&sy);
}

inline
void AXPY(int size,const std::complex<double> &a,const std::complex<double>* x, int sx, 
	  std::complex<double>* y, int sy){
  zaxpy_(&size,&a,x,&sx,y,&sy);
}
// ----------------------------------------------------------------------------
inline
void COPY(int size,const float* x,int sx,float* y,int sy) {
  scopy_(&size,x,&sx,y,&sy);
}

inline
void COPY(int size,const double* x,int sx,double* y,int sy) {
  dcopy_(&size,x,&sx,y,&sy);
}

inline
void COPY(int size,const std::complex<float>* x,int sx,std::complex<float>* y,int sy) {
  ccopy_(&size,x,&sx,y,&sy);
}

inline
void COPY(int size,const std::complex<double>* x,int sx,std::complex<double>* y,int sy) {
  zcopy_(&size,x,&sx,y,&sy);
}
// ----------------------------------------------------------------------------
inline
void SCAL(int size,const float &a,float* y,int sy) {
  sscal_(&size,&a,y,&sy);
}

inline
void SCAL(int size,const double &a,double *y,int sy) {
  dscal_(&size,&a,y,&sy);
}

inline
void SCAL(int size,const std::complex<float> &a,std::complex<float>* y,int sy) {
  cscal_(&size,&a,y,&sy);
}

inline
void SCAL(int size,const std::complex<double> &a,std::complex<double>* y,int sy) {
  zscal_(&size,&a,y,&sy);
}


}      /* namespace BLAS */

#endif /* PSIMAG_BLAS */
