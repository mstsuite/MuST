#ifndef RANDOM_EVEC_H
#define RANDOM_EVEC_H

#include <cmath>
#include <random>

void inline random_evec(double ev[3])
{
  std::mt19937 rng;
  do {
    ev[0] = rng();
    ev[1] = rng();
  } while(ev[0]*ev[0]+ev[1]*ev[1]>1.0);
  ev[2] = rng();
  double r = std::sqrt((1.0-ev[2]*ev[2])/(ev[0]*ev[0]+ev[1]*ev[1]));
  ev[0] *= r;
  ev[1] *= r;
  if(rand()%2) ev[0] = -ev[0];
  if(rand()%2) ev[1] = -ev[1];
  if(rand()%2) ev[2] = -ev[2];
}

#ifdef XYZZ_MEIS

#ifdef ISING

    void inline random_evec_1(double ev[3])
  {
    ev[0]=ev[1]=0.0;
    ev[2]=1.0;
    if(rng()%2 == 0) ev[2]=-ev[2];
  }
#else
  void inline random_evec_1(double ev[3])
  {
    double x,y,z;
    do {
      x = rnd(rng);
      y = rnd(rng);
    } while(x*x+y*y>1);
    z = rnd(rng);
    double r = sqrt((1-z*z)/(x*x+y*y));
    x *= r;
    y *= r;
    if (rng() % 2 == 0) x = -x;
    if (rng() % 2 == 0) y = -y;
    if (rng() % 2 == 0) z = -z;
    r=1.0/sqrt(x*x+y*y+z*z);
    ev[0]=x*r; ev[1]=y*r; ev[2]=z*r;
  }
#endif

#ifdef ISING    
  void inline random_evec(double ev[3])
  {
    ev[2]=-ev[2];
  }
#else
  void inline random_evec(double ev[3])
  {
    double x, y, z;
    do {
      x = rnd(rng); y = rnd(rng);
    } while(x*x+y*y>1); 
    z = rnd(rng);
    double r = sqrt((1-z*z)/(x*x+y*y));
    x *= r; y*= r;
    if (rng() % 2 == 0) x = -x;
    if (rng() % 2 == 0) y = -y;
    if (rng() % 2 == 0) z = -z; 
    // Project out the parallel component;
    r = x*ev[0] + y*ev[1] + z*ev[2];
    x -= r*ev[0]; y -= r*ev[1]; z -= r*ev[2];
    r = x*x + y*y + z*z;
    double t = 1-0.3*rnd(rng);
    ev[0] *= t; ev[1] *= t; ev[2] *= t;
    r = sqrt((1-t*t)/r);
    ev[0] += x*r; ev[1] += y*r; ev[2] += z*r;
    r=1.0/sqrt(ev[0]*ev[0]+ev[1]*ev[1]+ev[2]*ev[2]);
    ev[0]*=r; ev[1]*=r; ev[2]*=r;
    
    /*  
    ev[2]=1.0-2.0*rnd(rng);
    // ev[2]=rnd11(rng);
    double phi=2.0*M_PI*rnd(rng);
    // double phi=rnd0pi(rng);
    double cos_theta=sqrt(1-ev[2]*ev[2]);
    ev[0]=cos_theta*cos(phi);
    ev[1]=cos_theta*sin(phi);
    */  
  }
#endif
#endif

#endif
