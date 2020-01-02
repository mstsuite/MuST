#ifndef CUDOUBLE_COMPLEX
#define CUDOUBLE_COMPLEX

class cudaDoubleComplex : public double2 {

  public:

    __inline__ __device__ __host__ cudaDoubleComplex() {}
    __inline__ __device__ __host__ cudaDoubleComplex(double x, double y) {
      this->x=x;
      this->y=y;
    }

    __inline__ __device__ __host__ cudaDoubleComplex& operator=(const double &a) {
      x=a;y=0;
      return *this;
    }
    __inline__ __device__ __host__ double real() { return x;}
    __inline__ __device__ __host__ double imag() { return y;}
};

__inline__ __device__ __host__  cudaDoubleComplex operator+(const cudaDoubleComplex &a, const double &b) {
  return cudaDoubleComplex(a.x+b,a.y);
}
__inline__ __device__ __host__ cudaDoubleComplex operator+(const double &b, const cudaDoubleComplex &a) {
  return cudaDoubleComplex(a.x+b,a.y);
}
__inline__ __device__ __host__ cudaDoubleComplex operator+(const cudaDoubleComplex &a, const cudaDoubleComplex &b) {
  return cudaDoubleComplex(a.x+b.x,a.y+b.y);
}

__inline__ __device__ __host__ cudaDoubleComplex operator-(const cudaDoubleComplex &a, const double &b) {
  return cudaDoubleComplex(a.x-b,a.y);
}
__inline__ __device__ __host__ cudaDoubleComplex operator-(const double &b, const cudaDoubleComplex &a) {
  return cudaDoubleComplex(b-a.x,-a.y);
}
__inline__ __device__ __host__ cudaDoubleComplex operator-(const cudaDoubleComplex &a, const cudaDoubleComplex &b) {
  return cudaDoubleComplex(a.x-b.x,a.y-b.y);
}

__inline__ __device__ __host__ cudaDoubleComplex operator*(const cudaDoubleComplex &a, const double &b) {
  return cudaDoubleComplex(a.x*b,a.y*b);
}
__inline__ __device__ __host__ cudaDoubleComplex operator*(const double &b, const cudaDoubleComplex &a) {
  return cudaDoubleComplex(a.x*b,a.y*b);
}
__device__ __host__ __inline__ cudaDoubleComplex operator*(const cudaDoubleComplex a, const cudaDoubleComplex b) {
  return cudaDoubleComplex(a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x);
}

__inline__ __device__ __host__ cudaDoubleComplex operator/(const cudaDoubleComplex &a, const double &b) {
  return cudaDoubleComplex(a.x/b,a.y/b);
}
__inline__ __device__ __host__ cudaDoubleComplex operator/(const double &b, const cudaDoubleComplex &a) {
  return cudaDoubleComplex(a.x/b,a.y/b);
}
    
//scaled implementaiton
__inline__ __device__ __host__ cudaDoubleComplex operator/(const cudaDoubleComplex &a, const cudaDoubleComplex &b) {

  double s = fabs(b.x) + fabs(b.y);
  double oos = 1.0 / s;
  double ars = a.x * oos;
  double ais = a.y * oos;
  double brs = b.x * oos;
  double bis = b.y * oos;
  s = (brs * brs) + (bis * bis);
  oos = 1.0 / s;
  return cudaDoubleComplex(((ars * brs) + (ais * bis)) * oos, ((ais * brs) - (ars * bis)) * oos);

}

__inline__ __device__ __host__ cudaDoubleComplex operator-(const cudaDoubleComplex &a) {
  return cudaDoubleComplex(-a.x,-a.y);
}

__inline__ __device__ __host__ cudaDoubleComplex conj(const cudaDoubleComplex &a) {
  return cudaDoubleComplex(a.x,-a.y);
}

__device__ inline cudaDoubleComplex exp(const cudaDoubleComplex &a) {
  double e=exp(a.x);
  return cudaDoubleComplex(e*cos(a.y),e*sin(a.y));
}

__device__ inline bool operator==(const cudaDoubleComplex &a, const cudaDoubleComplex &b) {
  return a.x==b.x && a.y==b.y;
}
__device__ inline bool operator!=(const cudaDoubleComplex &a, const cudaDoubleComplex &b) {
  return a.x!=b.x || a.y!=b.y;
}
#endif
