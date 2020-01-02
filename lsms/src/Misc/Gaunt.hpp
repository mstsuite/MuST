#ifndef BINOMIAL_HPP
#define BINOMIAL_HPP

#include <cmath>


namespace NewGaunt {
//
//        ( n )    ( n ) (n+1)-k
// return (   ) == (   ) -------
//        ( k )    (k-1)    k
//
template<typename Int>
Int binomialCoefficient(Int n, Int k)
{
  if(k>n || k<Int(0)) return Int(0);
  Int b=Int(1);
  for(Int i=Int(1); i<=k; i++)
    b=(b*((n+1)-i))/i;

  return b;
}

template<typename Int>
Int factorial(Int n)
{
  if(n<0) return Int(0);
  Int f=1;
  for(Int i=1; i<=n; i++)
    f=f*i;
  return f;
}

// helper function: gcd
template<typename Int>
Int gcd(Int i, Int j)
{
  Int a=i;
  Int b=j;
  Int h;

  while(b!=Int(0))
  {
    h=a%b;
    a=b;
    b=h;
  }

  return a;
}

//helper function iabs
template<typename Int>
Int iabs(Int i)
{
  if(i<Int(0)) return -i;
  return i;
}

template<typename Int>
Int w3j_intterm(Int j1, Int j2, Int j3, Int m1, Int m2, Int m3)
{
  Int sum=Int(0);
  Int I1 = Int(0);
  if(j1-j3+m2>I1) I1=j1-j3+m2;
  if(j2-j3-m1>I1) I1=j2-j3-m1;
  Int I2 = j1+j2-j3;
  if(j1-m1<I2) I2=j1-m1;
  if(j2+m2<I2) I2=j2+m2;

  Int sgn=Int(1);
  if(iabs(I1)%Int(2)==Int(1)) sgn=Int(-1);

  for(Int k=I1; k<=I2; k++)
  {
    sum+=sgn*binomialCoefficient(j1+j2-j3,k)*binomialCoefficient(j1-j2+j3,j1-m1-k)
            *binomialCoefficient(-j1+j2+j3,j2+m2-k);
    sgn*=Int(-1);
  }
  return sum;
}

template <typename Int>
Int w3j_sqrt_numerator(Int j1, Int j2, Int j3, Int m1, Int m2, Int m3)
{
  return binomialCoefficient(Int(2)*j1,j1-j2+j3)*binomialCoefficient(Int(2)*j2,j1+j2-j3)
        *binomialCoefficient(Int(2)*j3,-j1+j2+j3);
}

template<typename Int>
Int w3j_sqrt_denominator(Int j1, Int j2, Int j3, Int m1, Int m2, Int m3)
{
  return binomialCoefficient(Int(2)*j1,j1+m1)*binomialCoefficient(Int(2)*j2,j2+m2)
        *binomialCoefficient(Int(2)*j3,j3+m3);
}

template<typename R, typename Int>
R w3j_sqrt(Int j1, Int j2, Int j3, Int m1, Int m2, Int m3)
{
  Int n,d;
  Int g;
  n=binomialCoefficient(Int(2)*j1,j1-j2+j3);
  d=binomialCoefficient(Int(2)*j1,j1+m1);
  g=gcd(n,d);
  d=d/g;
  n=n/g;
  n=n*binomialCoefficient(Int(2)*j2,j1+j2-j3);
  d=d*binomialCoefficient(Int(2)*j2,j2+m2);
  g=gcd(n,d);
  d=d/g;
  n=n/g;
  n=n*binomialCoefficient(Int(2)*j3,-j1+j2+j3);
  d=d*binomialCoefficient(Int(2)*j3,j3+m3);
  g=gcd(n,d);
  d=d/g;
  n=n/g;

  return sqrt(R(n)/R(d));
}

template<typename R, typename Int>
R w3j_Delta(Int j1, Int j2, Int j3)
{
  Int n,d; // numerator and denominator
  Int g; // gcd used in intermediate calculations
  R dsqr; // Delta square
  // n = (j1+j2-j3)! (j1-j2+j3)! (-j1+j2+j3)!
  // d = (j1+j2+j3+1)!

  d=factorial(j1+j2+j3+Int(1));
  n=factorial(j1+j2-j3);
  g=gcd(n,d);
  d=d/g;
  n=n/g;
  n=n*factorial(j1-j2+j3);
  g=gcd(n,d);
  d=d/g;
  n=n/g;
  n=n*factorial(j2+j3-j1);
  g=gcd(n,d);
  d=d/g;
  n=n/g;

  dsqr=R(n)/R(d);
  return sqrt(dsqr);
}

template<typename R, typename Int>
R w3j(Int j1, Int j2, Int j3, Int m1, Int m2, Int m3)
{
  Int n,d,s,i,g;
  R delta, w3j_f;
  // printf("w3j(%d,%d,%d, %d,%d,%d)\n",j1,j2,j3,m1,m2,m3);
  if(m1+m2+m3!=Int(0)) return R(0);
  if((iabs(m1)>j1) || (iabs(m2)>j2) || (iabs(m3)>j3)) return R(0);
  if((j3<iabs(j1-j2)) || ((j1+j2)<j3)) return R(0);
  
  delta=w3j_Delta<R,Int>(j1,j2,j3);
/*
  n=w3j_sqrt_numerator<Int>(j1, j2, j3, m1, m2, m3);
  d=w3j_sqrt_denominator<Int>(j1, j2, j3, m1, m2, m3);
  // printf("d=%d n=%d i=%d s=%d delta=%lf\n",d,n,i,s,delta);
  g=gcd(n,d);
  d=d/g;
  n=n/g;
*/
  i=w3j_intterm<Int>(j1, j2, j3, m1, m2, m3);
  s=Int(1);
  if(iabs(j1-j2-m3)%Int(2) == Int(1)) s=Int(-1);

  // printf("g=%d -> d=%d n=%d i=%d s=%d delta=%lf\n",g,d,n,i,s,delta);

  return delta*w3j_sqrt<R,Int>(j1, j2, j3, m1, m2, m3)*R(s*i);
}

} // namepsape NewGaunt

template<typename R, typename Int>
R gaunt(Int lp, Int l1, Int l2, Int mp, Int m1, Int m2)
{
  if((lp+l1+l2)%Int(2)==Int(1)) return R(0);
  if(NewGaunt::iabs(mp)>lp || NewGaunt::iabs(m1)>l1 || NewGaunt::iabs(m2)>l2) return R(0);
  R g=sqrt(R((2*lp+1)*(2*l1+1)*(2*l2+1))/(4.0*M_PI))
      *NewGaunt::w3j<R,Int>(lp,l1,l2,0,0,0)*NewGaunt::w3j<R,Int>(lp,l1,l2,-mp,m1,m2);
  if(NewGaunt::iabs(mp)%Int(2)==Int(1)) return -g;
  return g;
}
#endif
