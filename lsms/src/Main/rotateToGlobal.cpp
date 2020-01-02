// rotate results from the relativistic Green's function calculation
// green_function_rel.f to the global frame of reference
// this is from LSMS_1.9 gettau_c.f after the call to green_function_rel.

// this is rotr from LSMS_1.9
// ! rotation matrix for 3D vectors 
// ! Simon L. Altmann: Rotations,Quaternions,..., p.75, Eq.(3.3.11)
// ! input: 
// !        tvec   normal vector of axis
// !        phi    angle of rotation
// ! output: 
// !        drot   matrix of rotation
// !

#include "Complex.hpp"
#include <cmath>
#include "Matrix.hpp"
#include "Array3d.hpp"
#include "../SingleSite/AtomData.hpp"

void calculateRotationMatrix(Matrix<Real> &drot, Real *tvec, Real phi)
{
  Real sp=sin(phi);
  Real sp2=sin(0.5*phi);
  Real tx=tvec[0]*sp;
  Real ty=tvec[1]*sp;
  Real tz=tvec[2]*sp;
  Real tx2=tvec[0]*sp2;
  Real ty2=tvec[1]*sp2;
  Real tz2=tvec[2]*sp2;

  drot(0,0) = 1.0 - 2.0 * (ty2*ty2 + tz2*tz2);
  drot(1,1) = 1.0 - 2.0 * (tx2*tx2 + tz2*tz2);
  drot(2,2) = 1.0 - 2.0 * (tx2*tx2 + ty2*ty2);
  drot(0,1) = -tz + 2.0*tx2*ty2;
  drot(1,0) = tz + 2.0*tx2*ty2;
  drot(0,2) = ty + 2.0*tx2*tz2;
  drot(2,0) = -ty + 2.0*tx2*tz2;
  drot(1,2) = -tx + 2.0*ty2*tz2;
  drot(2,1) = tx + 2.0*ty2*tz2;
}

void rotateToGlobal(AtomData &atom, Matrix<Complex> &dos, Matrix<Complex> &dosck,
                    Matrix<Complex> &dos_orb, Matrix<Complex> &dosck_orb,
                    Array3d<Complex> &green, Array3d<Complex> &dens_orb, int i)
{
  Real axis[3];
  Matrix<Real> rot(3,3);
  
  Real a=std::sqrt(atom.evec[0]*atom.evec[0] + atom.evec[1]*atom.evec[1]);
  Real phi = std::acos(atom.evec[2]);

  if(a == 0.0)
  {
    axis[0]=axis[1]=0.0;
    axis[2]=copysign(1.0,atom.evec[2]);
    a=1.0;
  } else {
    axis[0] = -atom.evec[1];
    axis[1] = atom.evec[0];
    axis[2] = 0.0;
    a=1.0/a;
  }
  axis[0]*=a;
  axis[1]*=a;
  axis[2]*=a;

  calculateRotationMatrix(rot, axis, phi);
  Complex t1 = rot(0,0)*dos(1,i)+rot(0,1)*dos(2,i)+rot(0,2)*dos(3,i);
  Complex t2 = rot(1,0)*dos(1,i)+rot(1,1)*dos(2,i)+rot(1,2)*dos(3,i);
  Complex t3 = rot(2,0)*dos(1,i)+rot(2,1)*dos(2,i)+rot(2,2)*dos(3,i);
  dos(1,i) = t1;
  dos(2,i) = t2;
  dos(3,i) = t3;
  
  t1 = rot(0,0)*dosck(1,i)+rot(0,1)*dosck(2,i)+rot(0,2)*dosck(3,i);
  t2 = rot(1,0)*dosck(1,i)+rot(1,1)*dosck(2,i)+rot(1,2)*dosck(3,i);
  t3 = rot(2,0)*dosck(1,i)+rot(2,1)*dosck(2,i)+rot(2,2)*dosck(3,i);
  dosck(1,i) = t1;
  dosck(2,i) = t2;
  dosck(3,i) = t3;
  
  t1 = rot(0,0)*dos_orb(0,i)+rot(0,1)*dos_orb(1,i)+rot(0,2)*dos_orb(2,i);
  t2 = rot(1,0)*dos_orb(0,i)+rot(1,1)*dos_orb(1,i)+rot(1,2)*dos_orb(2,i);
  t3 = rot(2,0)*dos_orb(0,i)+rot(2,1)*dos_orb(1,i)+rot(2,2)*dos_orb(2,i);
  dos_orb(0,i) = t1;
  dos_orb(1,i) = t2;
  dos_orb(2,i) = t3;
  
  t1 = rot(0,0)*dosck_orb(0,i)+rot(0,1)*dosck_orb(1,i)+rot(0,2)*dosck_orb(2,i);
  t2 = rot(1,0)*dosck_orb(0,i)+rot(1,1)*dosck_orb(1,i)+rot(1,2)*dosck_orb(2,i);
  t3 = rot(2,0)*dosck_orb(0,i)+rot(2,1)*dosck_orb(1,i)+rot(2,2)*dosck_orb(2,i);
  dosck_orb(0,i) = t1;
  dosck_orb(1,i) = t2;
  dosck_orb(2,i) = t3;

  for(int j=0; j<green.l_dim1(); j++)
  {
    // do j=1,jws
    t1 = rot(0,0)*green(j,1,i)+rot(0,1)*green(j,2,i)+rot(0,2)*green(j,3,i);
    t2 = rot(1,0)*green(j,1,i)+rot(1,1)*green(j,2,i)+rot(1,2)*green(j,3,i);
    t3 = rot(2,0)*green(j,1,i)+rot(2,1)*green(j,2,i)+rot(2,2)*green(j,3,i);
    green(j,1,i) = t1;
    green(j,2,i) = t2;
    green(j,3,i) = t3;
    t1 = rot(0,0)*dens_orb(j,0,i)+rot(0,1)*dens_orb(j,1,i)+rot(0,2)*dens_orb(j,2,i);
    t2 = rot(1,0)*dens_orb(j,0,i)+rot(1,1)*dens_orb(j,1,i)+rot(1,2)*dens_orb(j,2,i);
    t3 = rot(2,0)*dens_orb(j,0,i)+rot(2,1)*dens_orb(j,1,i)+rot(2,2)*dens_orb(j,2,i);
    dens_orb(j,0,i) = t1;
    dens_orb(j,1,i) = t2;
    dens_orb(j,2,i) = t3;
  }
}
