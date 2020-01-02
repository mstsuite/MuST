#include "Real.hpp"
#include "rationalFit.hpp"
#include <vector>


// #define DEFAULT_STEP_SUBDIVISION 3
#define DEFAULT_STEP_SUBDIVISION 0

template <class Fn>
Real integrateStep(Real x0, Real x1, Fn &fn)
{
  const Real sixth=1.0/6.0;
  const Real nintieth=1.0/90.0;
  Real h = x1-x0;
// trapezoid:
//  return 0.5*h*(fn(x0)+fn(x1));

// Simpson:
  return sixth*h*(fn(x0)+4.0*fn(0.5*(x0+x1))+fn(x1));

// Newton-Cotes of 4-th order:
//  return nintieth*h*(7.0*fn(x0)+32.0*fn(x0+0.25*h)+12.0*fn(0.5*(x0+x1))+32.0*fn(x0+0.75*h)+7.0*fn(x1));
}

// calculate the 1-d integral of an integrand defined on the grid from 
template <size_t stepSubdivision=DEFAULT_STEP_SUBDIVISION>
void integrateOneDim(Real *grid, Real *integrand, Real *integral, size_t n)
{
  RationalFit<Real> fit;
  const Real invStepSubdivision=1.0/Real(stepSubdivision);
  integral[0]=0.0;
  if(stepSubdivision==0) // use fit.integral for integration steps
    for(size_t i=0; i<n-1; i++)
    {
      fit.set(grid,integrand,i,n);
      integral[i+1]=integral[i]+fit.integral(grid[i],grid[i+1]);
    }
  else if(stepSubdivision==1)
    for(size_t i=0; i<n-1; i++)
    {
      fit.set(grid,integrand,i,n);
      integral[i+1]=integral[i]+integrateStep(grid[i],grid[i+1],fit);
    }
  else
    for(size_t i=0; i<n-1; i++)
    {
      fit.set(grid,integrand,i,n);
      Real h=(grid[i+1]-grid[i])*invStepSubdivision;
      Real x0=grid[i]; Real x1=x0+h;
      Real stepIntegral=0.0;
      for(size_t j=0; j<stepSubdivision; j++)
      {
        stepIntegral+=integrateStep(x0,x1,fit);
        x0=x1; x1=x0+h;
      }
      integral[i+1]=integral[i]+stepIntegral;
    }
}

template <typename Func, size_t stepSubdivision=DEFAULT_STEP_SUBDIVISION>
void integrateOneDim(std::vector<Real> &grid, Func fn, std::vector<Real> &integral)
{
  RationalFit<Real> fit;
  const Real invStepSubdivision=1.0/Real(stepSubdivision);
  integral[0]=0.0;
  if(stepSubdivision==0) // use fit.integral for integration steps
    for(size_t i=0; i<grid.size()-1; i++)
    {
      fit.set(grid,fn,i);
      integral[i+1]=integral[i]+fit.integral(grid[i],grid[i+1]);
    }
  else if(stepSubdivision==1)
    for(size_t i=0; i<grid.size()-1; i++)
    {
      fit.set(grid,fn,i);
      integral[i+1]=integral[i]+integrateStep(grid[i],grid[i+1],fit);
    }
    else
      for(size_t i=0; i<grid.size()-1; i++)
      {
        fit.set(grid,fn,i);
        Real h=(grid[i+1]-grid[i])*invStepSubdivision;
        Real x0=grid[i]; Real x1=x0+h;
        Real stepIntegral=0.0;
        for(size_t j=0; j<stepSubdivision; j++)
        {
          stepIntegral+=integrateStep(x0,x1,fit);
          x0=x1; x1=x0+h;
        }
        integral[i+1]=integral[i]+stepIntegral;
      }
}

template <size_t stepSubdivision=DEFAULT_STEP_SUBDIVISION>
inline void integrateOneDim(std::vector<Real> &grid, std::vector<Real> &integrand, std::vector<Real> &integral)
{
  integrateOneDim<stepSubdivision>(&grid[0], &integrand[0], &integral[0], integrand.size());
}

template <size_t stepSubdivision=DEFAULT_STEP_SUBDIVISION>
void integrateOneDimSpherical(Real *grid, Real *integrand, Real *integral, size_t n)
{
  RationalFit<Real> fit;
  const Real invStepSubdivision=1.0/Real(stepSubdivision);
  integral[0]=0.0;
  
  if(stepSubdivision==0) // use fit.integral for integration steps
    for(size_t i=0; i<n-1; i++)
    {
      fit.set(grid,integrand,[](Real r, Real f){return 4.0*M_PI*r*r*f;}, i, n);
      integral[i+1]=integral[i]+fit.integral(grid[i],grid[i+1]);
    }
  else if(stepSubdivision==1)
    for(size_t i=0; i<n-1; i++)
    {
      fit.set(grid,integrand,[](Real r, Real f){return 4.0*M_PI*r*r*f;}, i, n);
      integral[i+1]=integral[i]+integrateStep(grid[i],grid[i+1],fit);
    }
  else
    for(size_t i=0; i<n-1; i++)
    {
      fit.set(grid,integrand,[](Real r, Real f){return 4.0*M_PI*r*r*f;}, i, n);
      Real h=(grid[i+1]-grid[i])*invStepSubdivision;
      Real x0=grid[i]; Real x1=x0+h;
      Real stepIntegral=0.0;
      for(size_t j=0; j<stepSubdivision; j++)
      {
        stepIntegral+=integrateStep(x0,x1,fit);
        x0=x1; x1=x0+h;
      }
      integral[i+1]=integral[i]+stepIntegral;
    }
}

template <size_t stepSubdivision=DEFAULT_STEP_SUBDIVISION>
inline void integrateOneDimSpherical(std::vector<Real> &grid, std::vector<Real> &integrand, std::vector<Real> &integral)
{
  integrateOneDimSpherical<stepSubdivision>(&grid[0], &integrand[0], &integral[0], integrand.size());
}

template <size_t stepSubdivision=DEFAULT_STEP_SUBDIVISION>
Real integrateOneDim(std::vector<Real> &grid, std::vector<Real> &integrand, std::vector<Real> &integral, Real r)
{
  integrateOneDim<stepSubdivision>(grid, integrand, integral);
  return interpolate(grid, integral, r);
}

template <size_t stepSubdivision=DEFAULT_STEP_SUBDIVISION>
Real integrateOneDimSpherical(std::vector<Real> &grid, std::vector<Real> &integrand, std::vector<Real> &integral, Real r)
{
  integrateOneDimSpherical<stepSubdivision>(grid, integrand, integral);
  return interpolate(grid, integral, r);
}

template <size_t stepSubdivision=DEFAULT_STEP_SUBDIVISION>
void integrateOneDimRPower(Real *grid, Real *integrand, Real *integral, size_t n, int p)
{
  RationalFit<Real> fit;
  const Real invStepSubdivision=1.0/Real(stepSubdivision);
  integral[0]=0.0;
  
  if(stepSubdivision==0) // use fit.integral for integration steps
    for(size_t i=0; i<n-1; i++)
    {
      fit.set(grid,integrand,[p](Real r, Real f){return std::pow(r,p)*f;}, i, n);
      integral[i+1]=integral[i]+fit.integral(grid[i],grid[i+1]);
    }
  else if(stepSubdivision==1)
    for(size_t i=0; i<n-1; i++)
    {
      fit.set(grid,integrand,[p](Real r, Real f){return std::pow(r,p)*f;}, i, n);
      integral[i+1]=integral[i]+integrateStep(grid[i],grid[i+1],fit);
    }
  else
    for(size_t i=0; i<n-1; i++)
    {
      fit.set(grid,integrand,[p](Real r, Real f){return std::pow(r,p)*f;}, i, n);
      Real h=(grid[i+1]-grid[i])*invStepSubdivision;
      Real x0=grid[i]; Real x1=x0+h;
      Real stepIntegral=0.0;
      for(size_t j=0; j<stepSubdivision; j++)
      {
        stepIntegral+=integrateStep(x0,x1,fit);
        x0=x1; x1=x0+h;
      }
      integral[i+1]=integral[i]+stepIntegral;
    }
}

template <size_t stepSubdivision=DEFAULT_STEP_SUBDIVISION>
inline void integrateOneDimRPower(std::vector<Real> &grid, std::vector<Real> &integrand, std::vector<Real> &integral, int p)
{
  integrateOneDimRPower<stepSubdivision>(&grid[0], &integrand[0], &integral[0], integrand.size(), p);
}

template <size_t stepSubdivision=DEFAULT_STEP_SUBDIVISION>
Real integrateOneDimRPower(std::vector<Real> &grid, std::vector<Real> &integrand, std::vector<Real> &integral, Real r, int p)
{
  integrateOneDimRPower<stepSubdivision>(grid, integrand, integral, p);
  return interpolate(grid, integral, r);
}
