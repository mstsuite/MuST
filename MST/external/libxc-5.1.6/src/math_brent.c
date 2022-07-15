/*
 Copyright (C) 2019 M.A.L. Marques

 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

/*******************************************************************************
 based on C++ code by Grant Williams, that was by turn based on wikipedia's
 pseudocode
********************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "util.h"

GPU_FUNCTION
double xc_math_brent
(xc_brent_f f,
 double lower_bound, double upper_bound,
 double TOL, double MAX_ITER,
 void *f_params)
{
  double a, b, fa, fb, fs;
  double tmp;
  double c, fc, s, d;
  int mflag, iter;

  a = lower_bound;
  b = upper_bound;
  fa = f(a, f_params);
  fb = f(b, f_params);
  fs = 0;

  if (fa * fb > 0){
#ifndef HAVE_CUDA
    fprintf(stderr, "Brent: bracketing error [%lf,%lf]\n", a, b);
    exit(1);
#endif
  }

  if (fabs(fa) < fabs(fb)){
    tmp = a;   a =  b;  b = tmp; /* swap  a and  b */
    tmp = fa; fa = fb; fb = tmp; /* swap fa and fb */
  }

  c = a;     /* c now equals the largest magnitude of the lower and upper bounds */
  fc = fa;   /* precompute function evalutation for point c by assigning it the same value as fa */
  mflag = 1; /* boolean flag used to evaluate if statement later on */
  s = 0;     /* Our root that will be returned */
  d = 0;     /* Only used if mflag is unset (mflag == false) */

  for (iter=1; iter<MAX_ITER; ++iter){
    /* stop if converged or error is less than tolerance */
    if (fabs(b - a) < TOL)
      return (b + a)/2.0;

    if (fa != fc && fb != fc){
      /* use inverse quadratic interopolation */
      s = ( a * fb * fc / ((fa - fb) * (fa - fc)) )
        + ( b * fa * fc / ((fb - fa) * (fb - fc)) )
        + ( c * fa * fb / ((fc - fa) * (fc - fb)) );
    }else{
      /* secant method */
      s = b - fb * (b - a) / (fb - fa);
    }

    /*
      (condition 1) s is not between  (3a+b)/4  and b or
      (condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
      (condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
      (condition 4) (mflag is set and |b−c| < |TOL|) or
      (condition 5) (mflag is false and |c−d| < |TOL|)
        */
    if(((s < (3 * a + b) * 0.25) || (s > b) ) ||
       (  mflag && (fabs(s-b) >= (fabs(b-c) * 0.5))) ||
       ( !mflag && (fabs(s-b) >= (fabs(c-d) * 0.5))) ||
       (  mflag && (fabs(b-c) < TOL)) ||
       ( !mflag && (fabs(c-d) < TOL))){
      /* bisection */
      s = (a+b)*0.5;
      mflag = 1;
    }else
      mflag = 0;

    fs = f(s, f_params);
    d = c;
    c = b;
    fc = fb;

    if (fa*fs < 0){
      b = s;
      fb = fs;
    }else{
      a = s;
      fa = fs;
    }

    if (fabs(fa) < fabs(fb)){
      tmp = a;   a =  b;  b = tmp; /* swap  a and  b */
      tmp = fa; fa = fb; fb = tmp; /* swap fa and fb */
    }

  }

#ifndef HAVE_CUDA
  fprintf(stderr, "Warning: Convergence not reached in brent\n");
#endif

  return (b + a)/2.0;

}
