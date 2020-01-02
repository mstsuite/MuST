#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP

// try gmp:
#include <gmp.h>

#include <cmath>


namespace NewGaunt {
//
//        ( n )    ( n ) (n+1)-k
// return (   ) == (   ) -------
//        ( k )    (k-1)    k
//
  void binomialCoefficient(mpz_t b, long n, long k)
  {
    mpz_set_ui(b,0);
    if(k>n || k<0) return;
    mpz_set_si(b,n);
    mpz_bin_ui(b,b,k);
  }

//helper function iabs
// provided as
// void mpz_abs (mpz_t rop, mpz_t op)
template<typename Int>
Int iabs(Int i)
{
  if(i<Int(0)) return -i;
  return i;
}

void w3j_intterm(mpz_t sum, long j1, long j2, long j3, long m1, long m2, long m3)
{
  mpz_t term,h;
  mpz_init(term);
  mpz_init(h);

  long I1 = 0;
  if(j1-j3+m2>I1) I1=j1-j3+m2;
  if(j2-j3-m1>I1) I1=j2-j3-m1;
  long I2 = j1+j2-j3;
  if(j1-m1<I2) I2=j1-m1;
  if(j2+m2<I2) I2=j2+m2;

  mpz_set_ui(sum,0);

  long sgn=1;
  if(iabs(I1)%2==1) sgn=-1;

  for(long k=I1; k<=I2; k++)
  {
//    sum+=sgn*binomialCoefficient(j1+j2-j3,k)*binomialCoefficient(j1-j2+j3,j1-m1-k)
//            *binomialCoefficient(-j1+j2+j3,j2+m2-k);
    binomialCoefficient(term, j1+j2-j3, k);
    binomialCoefficient(h,j1-j2+j3,j1-m1-k);
    mpz_mul(term,term,h);
    binomialCoefficient(h,-j1+j2+j3,j2+m2-k);
    mpz_mul(term,term,h);
    mpz_mul_si(term,term,sgn);
    mpz_add(sum,sum,term);
    sgn*= -1;
  }
  mpz_clear(h);
  mpz_clear(term);
}

void w3j_sqrt_sq(mpq_t r, long j1, long j2, long j3, long m1, long m2, long m3)
{
  mpz_t n,d,h;
  mpz_init(n);
  mpz_init(d);
  mpz_init(h);

  binomialCoefficient(n, 2*j1,j1-j2+j3);
  binomialCoefficient(h, 2*j2,j1+j2-j3);
  mpz_mul(n,n,h);
  binomialCoefficient(h, 2*j3,-j1+j2+j3);
  mpz_mul(n,n,h);;

  mpq_set_z(r,n);

  binomialCoefficient(d, 2*j1,j1+m1);
  binomialCoefficient(h, 2*j2,j2+m2);
  mpz_mul(d,d,h);
  binomialCoefficient(h, 2*j3,j3+m3);
  mpz_mul(d,d,h);

  mpq_set_den(r,d);
  mpq_canonicalize(r);

  mpz_clear(d); mpz_clear(h); mpz_clear(n);
}

void w3j_Delta_sq(mpq_t r, long j1, long j2, long j3)
{
  // n = (j1+j2-j3)! (j1-j2+j3)! (-j1+j2+j3)!
  // d = (j1+j2+j3+1)!

  mpz_t h;
  mpz_init(h);
  mpq_set_ui(r,1,1);

  mpz_fac_ui(h,j1+j2-j3);
  mpz_mul(mpq_numref(r),mpq_numref(r),h);
  mpz_fac_ui(h,j1+j2+j3+1);
  mpz_mul(mpq_denref(r),mpq_denref(r),h);
  mpq_canonicalize(r);
  mpz_fac_ui(h,j1-j2+j3);
  mpz_mul(mpq_numref(r),mpq_numref(r),h);
  mpz_fac_ui(h,j2+j3-j1);
  mpz_mul(mpq_numref(r),mpq_numref(r),h);

  mpq_canonicalize(r);

  mpz_clear(h);
}

void w3j(mpf_t w, long j1, long j2, long j3, long m1, long m2, long m3)
{
  mpq_t delta_sq,r;
  mpz_t i;
  mpf_t h;

  mpq_init(delta_sq);
  mpq_init(r);
  mpz_init(i);
  mpf_init(h);
  mpq_set_si(r,0,1);

  if(m1+m2+m3!=0) return;
  if((iabs(m1)>j1) || (iabs(m2)>j2) || (iabs(m3)>j3)) return;
  if((j3<iabs(j1-j2)) || ((j1+j2)<j3)) return;

  w3j_Delta_sq(delta_sq, j1, j2, j3);
  w3j_intterm(i, j1, j2, j3, m1, m2, m3);
  if(iabs(j1-j2-m3)%2 == 1) mpz_neg(i,i);

  w3j_sqrt_sq(r, j1, j2, j3, m1, m2, m3);

  mpq_mul(r,r,delta_sq);
  mpf_set_q(w,r);
  mpf_sqrt(w,w);
  mpf_set_z(h,i);
  mpf_mul(w,w,h);

  mpf_clear(h);
  mpz_clear(i);
  mpq_clear(r);
  mpq_clear(delta_sq);
}

} // namepsape NewGaunt

template<typename R, typename Int>
R gaunt(Int lp, Int l1, Int l2, Int mp, Int m1, Int m2)
{
  R gg;
  mpf_t g,h;

  if((lp+l1+l2)%Int(2)==Int(1)) return R(0);
  if(NewGaunt::iabs(mp)>lp || NewGaunt::iabs(m1)>l1 || NewGaunt::iabs(m2)>l2) return R(0);
  mpf_init(g);
  mpf_init(h);
  NewGaunt::w3j(g,lp,l1,l2,0,0,0);
  NewGaunt::w3j(h,lp,l1,l2,-mp,m1,m2);
  mpf_mul(g,g,h);
  mpf_set_si(h,(2*lp+1)*(2*l1+1)*(2*l2+1));
  mpf_sqrt(h,h);
  mpf_mul(g,g,h);

  gg=mpf_get_d(g)/sqrt(4.0*M_PI);
  if(NewGaunt::iabs(mp)%Int(2)==Int(1)) gg=-gg;
  mpf_clear(g);
  mpf_clear(h);
  return gg;
}
#endif
