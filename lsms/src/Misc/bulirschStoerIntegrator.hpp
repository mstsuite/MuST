// bulischStoerIntegrator.hpp
// integrate the initial value problem of a system of coupled ODEs
// using the Bulirsch-Stoer method
// R. Bulirsch and J. Stoer, Numerische Matematik 8, 1-13 (1966)
//

#include <functional>
#include <vector>
#include <cmath>

template<typename Rx, typename Ry>
void modifiedMidpoint(Rx x0, Rx x1, Ry *y0,
                      Ry *y1, int n, std::function<void(Rx x, Ry* y, Ry* dy)> rhs,
                      int steps)
{
  Ry y2[n], dy[n];
  Rx h = (x1-x0)/Rx(steps);
  Rx h2 = Rx(2)*h;

  Rx x= x0;

  rhs(x,y0,dy);
  for(int i=0; i<n; i++)
  {
    y1[i]=y0[i];
    y2[i]=y0[i] + h*dy[i];
  }
  for(int j=1; j<steps; j++)
  {
    x += h;
    rhs(x,y2,dy);
    for(int i=0; i<n; i++)
    {
      Ry temp = y1[i] + h2*dy[i];
      y1[i]=y2[i];
      y2[i]=temp;
    }
  }
  x += h;
  rhs(x,y2,dy);
  for(int i=0; i<n; i++)
  {
    y1[i]=Rx(0.5)*(y1[i]+y2[i]+h*dy[i]);
  }
}

template<typename Rh, typename Rf>
class ExtrapolatorT0 {
public:
  inline int idx(int d, int i, int j) { return d+dim*(i*(i+1)/2 + j); }
  std::vector<Rf> T;
  std::vector<Rh> h;
  int kNum,dim;
  ExtrapolatorT0(int d, int kMax) {T.resize(d*(kMax+1)*(kMax+2)/2); h.resize(kMax+1); kNum=0; dim=d;}
  void reset() {kNum=0;}
  void add(Rh hNew, Rf *fNew)
  {
    h[kNum]=hNew;
    int i=kNum;
    kNum++;
    for(int d=0; d<dim; d++)
    {
      T[idx(d,i,0)] = fNew[d];
      // printf("hNew=%g  h[i]=%g  ",hNew,h[i]);
      for(int k=1; k<kNum; k++)
      {
        Rh r = h[i-k]/h[i];
        Rh den = Rh(1)/(r*r - Rh(1));
        // printf("den=%g\n",den);
        T[idx(d,i,k)] = T[idx(d,i,k-1)] + (T[idx(d,i,k-1)] - T[idx(d,i-1,k-1)])*den;
      }
    }
  }
  Rh error()
  {
    if(kNum<2) return Rh(dim);
    Rh err=Rh(0);
    for(int d=0; d<dim; d++)
      err+=std::abs(T[idx(d,kNum-1,kNum-1)] - T[idx(d,kNum-2,kNum-2)]);
    return err;
  }
  Rh relativeError()
  {
    Rh mag=Rh(0);
    for(int d=0; d<dim; d++)
      mag+=std::abs(T[idx(d,kNum-1,kNum-1)]);
    if(mag==0.0) mag=1.0e-15;
    return error()/mag;
  }
  void estimate(Rf *f)
  {
    for(int d=0; d<dim; d++)
      f[d]=T[idx(d,kNum-1,kNum-1)];
  }
};

// integrate a system of n coupled ordinary differential equations
// from x0 to x1 with initial values in y0
// return values at x1 in y1
// rhs is a function that can evaluate the right hand side of the system of ODEs
// anywhere in the interval [x0,x1]
// eps is the desired max error
// returns an error code: 0 -> success
// templated on the type of real numbers
template <typename Rx, typename Ry>
int bulirschStoerIntegrator(Rx x0, Rx x1, Ry *y0,
                            Ry *y1, int n, std::function<void(Rx x, Ry* y, Ry* dy)> rhs,
                            Rx eps=1.0e-12)
{
  const int nSteps=11;
  const int stepTarget=7;
  const Rx stepGrowth=1.25;
  const Rx stepShrink=0.8;
  // Bulirsch Stoer sequence (n_j = 2 n_{j-2})
  const int stepSequence[]={2,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024};
  // Deuflhard sequence (n_j = 2 j)
  // const int stepSequence[]={2,4,6,8,10,12,14,16,18,20,22,24,26,28,30};
  // Romberg sequence (n_j = 2^j)
  // const int stepSequence[]={2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536};
  // Harmonic sequence (n_j = j)
  // const int stepSequence[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14}

  Rx xEpsilon = 1.0e-12 * std::abs(x0);
  Rx err;
  Ry y[n];

  Rx x = x0;
  Rx h0 = x1-x0;
  Rx h=h0;
  for(int i=0; i<n; i++) y[i]=y0[i];
  
  ExtrapolatorT0<Rx,Ry> extrapolation(n,nSteps+1);

  do
  {
    // fill extrapolation Tableau
    int step=0;
    Rx err=HUGE_VAL;
    extrapolation.reset();
    while(step < nSteps && err > eps)
    {
      modifiedMidpoint<Rx,Ry>(x,x+h,y,y1,n,rhs,stepSequence[step]);
      extrapolation.add(h/Rx(stepSequence[step]), y1);
      err=extrapolation.relativeError();
      step++;
    }
    // printf("step %d err %g\n",step,err);
    if(err < eps)
    {
      x = x + h;
      if(step<stepTarget)
        h *= stepGrowth;
      else
        h *= Rx(stepSequence[stepTarget])/Rx(stepSequence[step]);
      h= std::copysign(std::min(std::abs(h), std::abs(x1-x)), h0);
      extrapolation.estimate(y);
    } else {
      // the error is to large -> reduce the interval
      h = Rx(0.25) * h;
      if(h+x==x)
      {
        printf("Bulirsch-Stoer didn't converge to %g (error=%g)\n",eps,err);
        extrapolation.estimate(y1);
        return 1;
      }
    }
  } while (std::abs(x-x1)>xEpsilon);
  

  extrapolation.estimate(y1);
  return 0;
}


template <typename Rx, typename Ry>
int generalizedBulirschStoerIntegrator(Rx x0, Rx x1, Ry *y0,
                            Ry *y1, int n, std::function<void(Rx x, Ry* y, Ry* dy)> rhs,
                                       std::function<void(Rx x0, Rx x1, Ry *y0, Ry *y1, int n,
                                                          std::function<void(Rx x, Ry* y, Ry* dy)> rhs,
                                                          int steps)> method,
                                       Rx eps=1.0e-12)
{
  const int nSteps=11;
  const int stepTarget=7;
  const Rx stepGrowth=1.25;
  const Rx stepShrink=0.8;
  // Bulirsch Stoer sequence (n_j = 2 n_{j-2})
  const int stepSequence[]={2,4,6,8,12,16,24,32,48,64,96,128,192,256,384,512,768,1024};
  // Deuflhard sequence (n_j = 2 j)
  // const int stepSequence[]={2,4,6,8,10,12,14,16,18,20,22,24,26,28,30};
  // Romberg sequence (n_j = 2^j)
  // const int stepSequence[]={2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536};
  // Harmonic sequence (n_j = j)
  // const int stepSequence[]={1,2,3,4,5,6,7,8,9,10,11,12,13,14}

  Rx xEpsilon = 1.0e-12 * std::abs(x0);
  Rx err;
  Ry y[n];

  Rx x = x0;
  Rx h0 = x1-x0;
  Rx h=h0;
  for(int i=0; i<n; i++) y[i]=y0[i];
  
  ExtrapolatorT0<Rx,Ry> extrapolation(n,nSteps+1);

  do
  {
    // fill extrapolation Tableau
    int step=0;
    Rx err=HUGE_VAL;
    extrapolation.reset();
    while(step < nSteps && err > eps)
    {
      method(x,x+h,y,y1,n,rhs,stepSequence[step]);
      extrapolation.add(h/Rx(stepSequence[step]), y1);
      err=extrapolation.relativeError();
      step++;
    }
    // printf("step %d err %g\n",step,err);
    if(err < eps)
    {
      x = x + h;
      if(step<stepTarget)
        h *= stepGrowth;
      else
        h *= Rx(stepSequence[stepTarget])/Rx(stepSequence[step]);
      h= std::copysign(std::min(std::abs(h), std::abs(x1-x)), h0);
      extrapolation.estimate(y);
    } else {
      // the error is to large -> reduce the interval
      h = Rx(0.25) * h;
      if(h+x==x)
      {
        printf("Bulirsch-Stoer didn't converge to %g (error=%g)\n",eps,err);
        extrapolation.estimate(y1);
        return 1;
      }
    }
  } while (std::abs(x-x1)>xEpsilon);
  

  extrapolation.estimate(y1);
  return 0;
}
