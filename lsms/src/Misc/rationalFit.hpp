/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_RATIONAL_FIT_HPP
#define LSMS_RATIONAL_FIT_HPP

#include "Real.hpp"

/// Fit a function as a rational function with a  and quadratic denominator
/// at point \f$ r_i\f$.
///
/// Define \f$ x=r-r_i \f$ and \f$ \Delta=r_{i+1}-r_i \f$
/// make the ansatz:
/// \f[
///   f(x) = f_i - c_0 + c_1 x + \frac{(c_2-c_1) x + c_0}{(1+c_3 x (\Delta - x)}
/// \f]
template<typename T>
class RationalFit {
public:
  T c[4];
  T r0, delta, f0;
  const T tol=1.0e-6;

  void setWithStride(T *r, T *f, int i0, int n, int stride)
  {
    int ip1=i0+1;
    r0=r[i0];
    f0=f[i0*stride];
    delta=r[i0+1]-r0;
    T dr=r[ip1]-r[i0];
    T drv=f[ip1*stride]-f[i0*stride];
    c[2]=drv/dr;

    if(n==2)
    {
      c[0]=f[i0*stride];
      c[1]=c[3]=0.0;
      return;
    }

    int i1=std::min(3,n-1);
    if(i0>0) i1=i0-1;

    int i2=std::max(1,n-3);
    if(i0+2<n) i2=i0+2;

    if(i1!=i2)
    {
      c[1]=(f[i2*stride]-f[ip1*stride])/(r[i2]-r[ip1]) - (f[i1*stride]-f[i0*stride])/(r[i1]-r[i0]);
      if(c[1]*(f[i2*stride]-f[i1*stride])*(r[i2]-r[i1]) > 0.0) c[1]=-c[1];
      T h1=r[i1]-r[i0];
      T c0=f[i1*stride]-f[i0*stride];
      T c1=c0-c[2]*h1;
      T c2=c0-c[1]*h1;
      h1=h1*(r[ip1]-r[i1]);
      T eqn12= -h1*c2;
      T h2=r[i2]-r[i0];
      c0=f[i2*stride]-f[i0*stride];
      T c3=c0-c[2]*h2;
      T c4=c0-c[1]*h2;
      h2=h2*(r[ip1]-r[i2]);
      T eqn22=-h2*c4;
      T gj=(c4-c2)*h1*h2;
      if(gj==0.0)
      {
        c[0]=f[i0*stride];
        c[3]=0.0;
      } else {
        c[0]=(c1*eqn22 - c3*eqn12)/gj;
        c[3]=(c1*h2-c3*h1)/gj;
      }

      gj=c[3]*dr*dr;
      if((gj > -4.0) && (fabs(gj)>1.0e-14) && (gj<1.0e+14))
      {
        c[0]=c[0]/c[3];
      } else {
        c[0]=f[i0*stride];
        c[3]=0.0;
      }
    } else { // i1==i2
      T fj=0.5*(f[i0*stride]-f[ip1*stride]);
      c[0]=f[i0*stride]-fj;
      c[1]=0.0;
      T rj=r[i0]+0.5*dr;
      if(fabs(r[i1]-rj) > fabs(r[i2]-rj)) i1=i2;
      T h1=r[i1]-r[i0];
      T c1=f[i1*stride]-f[i0*stride]-c[2]*h1;
      h1=h1*(r[ip1]-r[i1]);
      if(f[i1*stride]!=fj)
      {
        c[3]=-c1/((f[i1*stride]-fj)*h1);
        T gj=c[3]*dr*dr;
        if((gj > -4.0) && (fabs(gj)>1.0e-14) && (gj<1.0e+14))
        {
          c[0]=f[i0*stride];
          c[3]=0.0;
        }
      } else {
        c[0]=f[i0*stride];
        c[3]=0.0;
      }
    }
  }

  void set(T *r, T *f, int i0, int n)
  {
    int ip1=i0+1;
    r0=r[i0];
    f0=f[i0];
    delta=r[i0+1]-r0;
    T dr=r[ip1]-r[i0];
    T drv=f[ip1]-f[i0];
    c[2]=drv/dr;

    if(n==2)
    {
      c[0]=f[i0];
      c[1]=c[3]=0.0;
      return;
    }

    int i1=std::min(3,n-1);
    if(i0>0) i1=i0-1;

    int i2=std::max(1,n-3);
    if(i0+2<n) i2=i0+2;

    if(i1!=i2)
    {
      c[1]=(f[i2]-f[ip1])/(r[i2]-r[ip1]) - (f[i1]-f[i0])/(r[i1]-r[i0]);
      if(c[1]*(f[i2]-f[i1])*(r[i2]-r[i1]) > 0.0) c[1]=-c[1];
      T h1=r[i1]-r[i0];
      T c0=f[i1]-f[i0];
      T c1=c0-c[2]*h1;
      T c2=c0-c[1]*h1;
      h1=h1*(r[ip1]-r[i1]);
      T eqn12= -h1*c2;
      T h2=r[i2]-r[i0];
      c0=f[i2]-f[i0];
      T c3=c0-c[2]*h2;
      T c4=c0-c[1]*h2;
      h2=h2*(r[ip1]-r[i2]);
      T eqn22=-h2*c4;
      T gj=(c4-c2)*h1*h2;
      if(gj==0.0)
      {
        c[0]=f[i0];
        c[3]=0.0;
      } else {
        c[0]=(c1*eqn22 - c3*eqn12)/gj;
        c[3]=(c1*h2-c3*h1)/gj;
      }

      gj=c[3]*dr*dr;
      if((gj > -4.0) && (fabs(gj)>1.0e-14) && (gj<1.0e+14))
      {
        c[0]=c[0]/c[3];
      } else {
        c[0]=f[i0];
        c[3]=0.0;
      }
    } else { // i1==i2
      T fj=0.5*(f[i0]-f[ip1]);
      c[0]=f[i0]-fj;
      c[1]=0.0;
      T rj=r[i0]+0.5*dr;
      if(fabs(r[i1]-rj) > fabs(r[i2]-rj)) i1=i2;
      T h1=r[i1]-r[i0];
      T c1=f[i1]-f[i0]-c[2]*h1;
      h1=h1*(r[ip1]-r[i1]);
      if(f[i1]!=fj)
      {
        c[3]=-c1/((f[i1]-fj)*h1);
        T gj=c[3]*dr*dr;
        if((gj > -4.0) && (fabs(gj)>1.0e-14) && (gj<1.0e+14))
        {
          c[0]=f[i0];
          c[3]=0.0;
        }
      } else {
        c[0]=f[i0];
        c[3]=0.0;
      }
    }
  }
  
  void set(std::vector<T> &r, std::vector<T> &f, int i0)
  {
    set(&r[0], &f[0], i0, f.size());
  }
// this allows us to provide a function to be applied to f
// fn: T, T -> T; fn(r,f)
  template<typename Func>
  void set(T *r, T *f, Func fn, int i0, int n)
  {
    std::vector<T> ff;
    ff.resize(n);
    for(int i=0; i<n; i++)
      ff[i]=fn(r[i],f[i]);

    set(r, &ff[0], i0, n);
  }

  template<typename Func>
  inline void set(std::vector<T> &r, std::vector<T> &f, Func fn, int i0)
  {
    set(&r[0], &f[0], fn, i0, f.size());
  }

  // calculate the coefficients of the rational fit for function fn(r[i]) at i0 
  template<typename Func>
  void set(std::vector<T> &r, Func fn, int i0)
  {
    std::vector<T> f;
    f.resize(r.size());
    for(int i=0; i<r.size(); i++)
      f[i]=fn(r[i],i);

    set(&r[0], &f[0], i0, r.size());
  }

  // obtain the value of r
  T operator()(T r)
  {
    T x=r-r0;
    return f0-c[0]+c[1]*x+((c[2]-c[1])*x+c[0])/(1.0+c[3]*x*(delta-x));
  }

  // calculate the derivative at r
  T derivative(T r)
  {
    T x=r-r0;
    T den=1.0+c[3]*x*(delta-x);
    T fac=c[2]*delta-c[0];
    return c[1] + ((c[2]-c[1])*den - ((c[2]-c[1])*x+c[0])*c[3]*(delta-2.0*x))/(den*den);
  }

  // calculate the integral of the fit from r1 to r2
  T integral(T r1, T r2)
  {
    T x1=r1-r0;
    T x2=r2-r0;
    // polynomial terms:
    T i1=(f0-c[0])*(x2-x1) + 0.5*c[1]*(x2*x2-x1*x1);
    // c_0/(1+c_3\Delta x -c_3 x^2) term:
    T i2=0.0;
    if(c[3]==0.0)
    {
      i2=c[0]*(x2-x1) + 0.5*(c[2]-c[1])*(x2*x2-x1*x1);
    } else {
      // see Gradshteyn / Ryzhik p 57, 2.103. 5.
      T m=c[2]-c[1];
      T n=c[0];
      // T a=1.0;
      T b=0.5*c[3]*delta;
      T cc=-c[3];
      if(cc>b*b) // ac>b^2
      {
        i2 = ((m/(2.0*cc))*std::log(std::fabs(1.0+2.0*b*x2+cc*x2*x2))
              + ((n*cc-m*b)/(cc*std::sqrt(cc-b*b)))*std::atan((cc*x2+b)/std::sqrt(cc-b*b)))
          - ((m/(2.0*cc))*std::log(std::fabs(1.0+2.0*b*x1+cc*x1*x1))
             + ((n*cc-m*b)/(cc*std::sqrt(cc-b*b)))*std::atan((cc*x1+b)/std::sqrt(cc-b*b)));
      } else { // ac<b^2
        T b2ac=std::sqrt(b*b-cc);
        i2 = ((m/(2.0*cc))*std::log(std::fabs(1.0+2.0*b*x2+cc*x2*x2))
              +((n*cc-m*b)/(2.0*cc*b2ac))*std::log(std::fabs((cc*x2+b-b2ac)/(cc*x2+b+b2ac))))
          -((m/(2.0*cc))*std::log(std::fabs(1.0+2.0*b*x1+cc*x1*x1))
            +((n*cc-m*b)/(2.0*cc*b2ac))*std::log(std::fabs((cc*x1+b-b2ac)/(cc*x1+b+b2ac))));
      }
    }
    
    return i1+i2;
  }
};

/// Given a table of function values r[i] -> f(r[i])
/// find the interpolated value of f(x)
template<typename T>
T interpolate(std::vector<T> &r, std::vector<T> &f, T x)
{
// find index i such that r[i]<x<r[i+1]
  if(r[0]>x) return f[0];
  int n=r.size();
  int i0=0, i1=n;
  while(i1-i0>1)
  {
    int d=i0+(i1-i0)/2;
    if(r[d]>x) i1=d; else i0=d;
  }
  RationalFit<T> fit;
  fit.set(r,f,i0);
  return fit(x);
}

// Interpolate from a table of points on a grid xOrigin to a different grid xTarget
template<typename T>
void interpolateTable(std::vector<T> &xOrigin, std::vector<T> &fOrigin, std::vector<T> &xTarget, std::vector<T> &fTarget)
{
  fTarget.resize(xTarget.size());
  for(int i=0; i<xTarget.size(); i++)
  {
    fTarget[i] = interpolate(xOrigin, fOrigin, xTarget[i]);
  }
}

template<typename T>
void calculateDerivative(T *r, T *f, T *df, size_t n)
{
  RationalFit<T> fit;
  for(size_t i=0; i<n; i++)
  {
    fit.set(r,f,i,n);
    df[i]=fit.derivative(r[i]);
  }
}

template<typename T>
void calculateDerivative(T *r, T *f, T *df, size_t n, size_t strideF, size_t strideDf)
{
  RationalFit<T> fit;
  for(size_t i=0; i<n; i++)
  {
    fit.setWithStride(r,f,i,n, strideF);
    df[i*strideDf]=fit.derivative(r[i]);
  }
}

#endif
