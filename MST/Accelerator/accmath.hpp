#ifndef ACCMATH_HPP
#define ACCMATH_HPP

__host__ __device__ static __inline__ cuDoubleComplex cuCexp(cuDoubleComplex z) {
   //   double real_part = z.x; // cuCreal(z);
   // double imag_part = z.y; // cuCimag(z);
   // double e_real = exp(real_part);
   // double cos_imag, sin_imag;
   // sincos(imag_part, &sin_imag, &cos_imag);
   // return make_cuDoubleComplex(e_real * cos_imag, e_real * sin_imag);

   cuDoubleComplex res;
   double t = exp(z.x);       // Calculate e^x using CUDA's double exp function
   double s, c;
   sincos(z.y, &s, &c);       // Calculate sin(y) and cos(y) simultaneously

   res.x = t * c;             // Real part: e^x * cos(y)
   res.y = t * s;             // Imaginary part: e^x * sin(y)
  
   return res;
}

__device__ void computeLegendre(int lmax, double x, double *plm) {
   double zero=0.0;
   double one=1.0;
   double two=2.0;
   double tol=1e-12;

   if (lmax == 0) {
      plm[0]=one;
      return;
   }
   else if (one-fabs(x) <= tol) {
      int jmax=(lmax+1)*(lmax+2)/2;
      for (int j=0; j<jmax; j++) {
         plm[j]=zero;
      }
      if (x < zero) {
         for (int l=0; l<=lmax; l++) {
            int j=(l+1)*l/2; 
            plm[j]=one-two*(l%2);
         }
      } 
      else {
         for (int l=0; l<=lmax; l++) {
            int j=(l+1)*l/2; 
            plm[j]=one;
         }
      }
      return;
   }

   plm[0]=one;
   plm[1]=x;
   double somx2=sqrt((one-x)*(one+x));
   if (lmax == 1) {
      plm[2]=somx2;
   }
   else {
      //  ==========================================================
      //                           m       m
      //  calculate the first two P   and P
      //                           m       m+1
      //  ==========================================================
      double pmm;
      double fact;
      for (int m=1; m<lmax; m++) {
         pmm=one;
         fact=one;
         for (int i=0; i<m; i++) {
            pmm=pmm*fact*somx2;
            fact=fact+two;
         }
         int mm=(m+1)*(m+2)/2;
         plm[mm-1]=pmm;
         plm[mm+m]=x*(2*m+1)*pmm;
      }
      pmm=one;
      fact=one;
      for (int i=0; i<lmax; i++) {
         pmm=pmm*fact*somx2;
         fact=fact+two;
      }
      int mm=(lmax+1)*(lmax+2)/2;
      plm[mm-1]=pmm;
      //  ==========================================================
      //                      m        m
      //  calculate the rest P     to P
      //                      m+2      lmax
      //  ==========================================================
      for (int m=0; m<=lmax; m++) {
         for (int l=m+2; l<=lmax; l++) {
            int ll = (l-1)*l/2+m;
            plm[ll+l]=(x*(2*l-1)*plm[ll]-(l+m-1)*plm[ll-l+1])/double(l-m);
         }
      }
   }
}

__global__ void createUnitMatrixKernel(cuDoubleComplex *matrix_d, int n);
__global__ void copyBlockKernel(cuDoubleComplex *A, int nA, cuDoubleComplex *B, int nB);
__global__ void copySubBlockKernel(cuDoubleComplex *A, int nA, int i, int j, cuDoubleComplex *B, int nB);

void copySubBlockToMatrix(double _Complex *a, int *size_a, int *row, int *col, double _Complex *b, int *size_b);

#endif
