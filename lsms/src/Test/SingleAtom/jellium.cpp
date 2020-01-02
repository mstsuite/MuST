#include <cmath>
#include <stdio.h>
#include <vector>

typedef double Real;

class Atom {
public:
  int Z;
  int numCore, numSemicore, numValence;

  Real rws;
};

Real calculateWignerSeitzRadiusBCC(Real a)
{
  // Real atomicVolume = 0.5 * a*a*a;
  // return std::cbrt((3.0*atomicVolume)/(4.0*M_PI));
  return std::cbrt(3.0/(8.0*M_PI)) * a;
}

Real calculateWignerSeitzRadiusSC(Real a)
{
  return std::cbrt(3.0/(4.0*M_PI)) * a;
}

Real calculateWignerSeitzRadiusFCC(Real a)
{
  return std::cbrt(3.0/(16.0*M_PI)) * a;
}

int main()
{
  Real latticeConstant = 6.831; //

  Atom a;
  a.Z = 29; a.numCore = 10; a.numSemicore = 8; a.numValence = 11; // Copper
  a.rws =  calculateWignerSeitzRadiusFCC(latticeConstant);

  printf("Atom in Jellium\n");
  printf("Z = %d\nfcc: latticeConstant = %lf\n",a.Z,latticeConstant);
  printf("Wigner-Seitz radius  = %lf\n", a.rws);
  printf("Z/rws = %lf\n",Real(a.Z)/a.rws);

  return 0;
}
