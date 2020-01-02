#include "Complex.hpp"

void calculateTrPxTau(Complex *tau00, int kkrsz, Complex *tr_pxtau)
{
  Array3<Complex> p_vec(kkrsz,kkrsz,3);
  Array3<Complex> tau_vec(kkrsz,kkrsz,3);

  tr_pxtau[0]=tr_pxtau[1]=tr_pxtau[2]=Complex(0.0);


}
