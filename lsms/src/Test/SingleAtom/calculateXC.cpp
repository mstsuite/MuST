/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
// calculate the exchange correlation potential for a charge density

#include "Real.hpp"
#include <cmath>
#include <vector>

// calculate the Wigner-Seitz desnity parameter r_s
// r_s = (4 pi rho / 3)^{-1/3}
Real densityParameterRs(Real rho)
{
  return std::pow( (4.0*M_PI/3.0) * rho, -1.0/3.0);
}

// v_xv = e_xc + rho * d e_xc/d rho

/*
      function alpha2(rs,dz,sp,iexch,exchg)

c     von barth-hedin  exch-corr potential
c     j. phys. c5,1629(1972)
c
      data     ccp,rp,ccf,rf/0.0450d+00,21.0d+00,
     >                       0.02250d+00,52.9166820d+00/
c
c
  10  continue
      fm=2.0d+00**(4.0d+00/3.0d+00)-2.0d+00
      fdz = ((1.0d+00+dz)**(4.0d+00/3.0d+00)
     >     +(1.0d+00-dz)**(4.0d+00/3.0d+00)-2.0d+00)/fm
      ex=-0.916330d+00/rs
      exf=ex*2.0d+00**0.333333330d+00
      xp=rs/rp
      xf=rs/rf
      gp = (1.0d+00+xp**3)*log(1.0d+00+1.0d+00/xp)
     >    -xp*xp +xp/2.0d+00 - 0.333333330d+00
      gf = (1.0d+00+xf**3)*log(1.0d+00+1.0d+00/xf)
     >    -xf*xf +xf/2.0d+00 - 0.333333330d+00
      exc = ex-ccp*gp
      excf=exf-ccf*gf
      dedz= (4.0d+00/3.0d+00)*(excf-exc)
     >     *((1.0d+00+dz)**(1.0d+00/3.0d+00)
     >     -(1.0d+00-dz)**(1.0d+00/3.0d+00))/fm
      gpp = 3.0d+00*xp*xp*log(1.0d+00+1.0d+00/xp)-1.0d+00/xp
     >     +1.50d+00-3.0d+00*xp
      gfp = 3.0d+00*xf*xf*log(1.0d+00+1.0d+00/xf)-1.0d+00/xf
     >     +1.50d+00-3.0d+00*xf
      depd=-ex/rs-ccp/rp*gpp
      defd=-exf/rs-ccf/rf*gfp
      decd=depd+(defd-depd)*fdz
c     exchange-correlation energy
      exchg= exc + (excf-exc)*fdz
c     exchange-correlation potential
      alpha2 = exc+(excf-exc)*fdz-rs*decd/3.0d+00
     >        +sp*(1.0d+00-sp*dz)*dedz
      return
c
*/
Real alpha2(Real rs, Real dz, Real sp, Real &eXC)
{
//    von barth-hedin  exch-corr potential
//    j. phys. c5,1629(1972)
  const Real ccp = 0.0450;
  const Real rp = 21.0;
  const Real ccf = 0.02250;
  const Real rf = 52.9166820;
  
  Real fm=std::pow(2.0, (4.0/3.0))-2.0;
  Real fdz = (std::pow((1.0+dz),(4.0/3.0))
              + std::pow((1.0-dz),(4.0d+00/3.0d+00)) - 2.0)/fm;
  Real ex = -0.916330/rs;
  Real exf = ex * std::pow(2.0,1.0/3.0);
  Real xp = rs/rp;
  Real xf = rs/rf;
  Real gp = std::pow(1.0+xp,3) * std::log(1.0 + 1.0/xp)
    - xp*xp + xp/2.0 - 0.333333330;
  Real gf = std::pow(1.0+xf,3) * std::log(1.0 + 1.0/xf)
    - xf*xf + xf/2.0 - 0.333333330;
  Real exc = ex - ccp*gp;
  Real excf = exf - ccf*gf;
  Real dedz= (4.0/3.0) * (excf-exc)
    * (std::pow(1.0 + dz, 1.0/3.0)
       - std::pow(1.0 - dz, 1.0/3.0))/fm;
  Real gpp = 3.0 * xp*xp * std::log(1.0 + 1.0/xp) - 1.0/xp
    + 1.50 - 3.0*xp;
  Real gfp = 3.0 * xf*xf * std::log(1.0 + 1.0/xf) - 1.0/xf
    + 1.50 - 3.0*xf;
  Real depd=-ex/rs-ccp/rp*gpp;
  Real defd=-exf/rs-ccf/rf*gfp;
  Real decd=depd+(defd-depd)*fdz;
//     exchange-correlation energy
  eXC = exc + (excf-exc)*fdz;
//     exchange-correlation potential
  return exc + (excf-exc) * fdz - rs*decd/3.0 + sp*(1.0 - sp*dz)*dedz;
}

// spin polarization weighting function from Cachiyo2016
Real fChachiyo2016(Real zeta)
{
  const Real c = 1.0/(2.0*(std::cbrt(2) - 1.0));
  return c*(std::pow(1.0+zeta, 4.0/3.0) + std::pow(1.0-zeta, 4.0/3.0) -2.0);
}

Real exchangeLDA(Real rho, Real &ex)
{
  ex = -6.0*std::pow(rho * 3.0/(64.0*M_PI), 1.0/3.0);
  return ex + ex/3.0;
}

// correlation energy from T. Chachiyo, J. Chem. Phys. 145, 021101 (2016).
Real chachiyo2016Correlation(Real rho, Real &ec)
{
  // total XC energy = int rho*exc dr
  // exchange energy density:
  //
  
  // correlation energy density:
  // e_c = a ln(1 + b/r_s + b/r_s^2)
  // constatnts from eq. 3 in T. Chachiyo, J. Chem. Phys. 145, 021101 (2016)
  Real const a = (std::log(2) - 1.0)/(M_PI * M_PI); // converted from Hartree in Chachiyo to Rydberg
  Real const b = 20.4562557;
  // r_s = (4 pi rho / 3) ^ -1/3 ->
  // 1/rs = (4 pi rho / 3) ^ 1/3
  Real rsInv = std::pow(4.0 * M_PI * rho / 3.0, 1.0/3.0);
  Real ec = a*std::log(1.0 + b*rsInv + b*rsInv*rsInv);
  
  // exchange correlation potential:
  // v_xc(r) = e_xc(rho(r)) + rho(r) * d e_xc(rho)/d rho
 
  Real rho_dec = a*b*(rsInv + 2.0*rsInv*rsInv)/(3.0*(1.0 + b*rsInv + b*rsInv*rsInv));

  return ec + rho_dec;
}

// correlation energy from T. Chachiyo, J. Chem. Phys. 145, 021101 (2016).
Real chachiyo2016XC(Real rho, Real &exc)
{
  // total XC energy = int rho*exc dr
  // exchange energy density:
  //
  //Real ex = 0;
  Real ex = -6.0*std::pow(rho * 3.0/(64.0*M_PI), 1.0/3.0);
  // correlation energy density:
  // e_c = a ln(1 + b/r_s + b/r_s^2)
  // constatnts from eq. 3 in T. Chachiyo, J. Chem. Phys. 145, 021101 (2016)
  Real const a = (std::log(2) - 1.0)/(M_PI * M_PI); // converted from Hartree in Chachiyo to Rydberg
  Real const b = 20.4562557;
  // r_s = (4 pi rho / 3) ^ -1/3 ->
  // 1/rs = (4 pi rho / 3) ^ 1/3
  Real rsInv = std::pow(4.0 * M_PI * rho / 3.0, 1.0/3.0);
  Real ec = a*std::log(1.0 + b*rsInv + b*rsInv*rsInv);
  exc = ex + ec;
  
  // exchange correlation potential:
  // v_xc(r) = e_xc(rho(r)) + rho(r) * d e_xc(rho)/d rho
  // rho_dex = rho * d e_x / d rho = 1/3 e_x
  Real rho_dex = ex/3.0;
  //
  Real rho_dec = a*b*(rsInv + 2.0*rsInv*rsInv)/(3.0*(1.0 + b*rsInv + b*rsInv*rsInv));

  return exc + rho_dex + rho_dec;
}

// correlation energy from T. Chachiyo, J. Chem. Phys. 145, 021101 (2016).
Real chachiyo2016PolarizedCorrelation(Real rho, Real &ec)
{
  // total XC energy = int rho*exc dr
  // 
  // correlation energy density:
  // e_c = a ln(1 + b/r_s + b/r_s^2)
  // constatnts from eq. 3 in T. Chachiyo, J. Chem. Phys. 145, 021101 (2016)
  Real const a = (std::log(2) - 1.0)/(2.0*M_PI * M_PI); // converted from Hartree in Chachiyo to Rydberg
  Real const b = 27.4203609;
  // r_s = (4 pi rho / 3) ^ -1/3 ->
  // 1/rs = (4 pi rho / 3) ^ 1/3
  Real rsInv = std::pow(4.0 * M_PI * rho / 3.0, 1.0/3.0);
  Real ec = a*std::log(1.0 + b*rsInv + b*rsInv*rsInv);
  
  // exchange correlation potential:
  // v_xc(r) = e_xc(rho(r)) + rho(r) * d e_xc(rho)/d rho
  
  Real rho_dec = a*b*(rsInv + 2.0*rsInv*rsInv)/(3.0*(1.0 + b*rsInv + b*rsInv*rsInv));

  return ec + rho_dec;
}

Real exchangeCorrelationPotentialLDAPolarized(Real rhoUp, Real rhoDown, Real r)
{
  Real eXC, eC_0, eC_1, vC_0, vC_1, rho, zeta;
  Real vXC, vX, vXup, vXdown, eX, eXup, eXdown;
  rho = rhoUp+rhoDown;
  zeta = (rhoUp-rhoDown)/rho;
  /*
  if(rho < 1.0e-9)
    return 0.0;
  // return alpha2(std::pow(3.0/rho , 1.0/3.0), 0.0, 1.0, eXC);
  return alpha2(std::pow(3.0*r*r/rho , 1.0/3.0), 0.0, 1.0, eXC);
  //*/
  vXup = exchangeLDA(2.0*rhoUp/(4.0*M_PI*r*r), eXup);
  vXdown = exchangeLDA(2.0*rhoDown/(4.0*M_PI*r*r), eXdown);
  eX = 0.5*(eXup + eXdown);
  vX = 0.5*(vXup + vXdown);

  vXC_0 = chachiyo2016(rho/(4.0*M_PI*r*r), eXC_0);
  vXC_1 = chachiyo2016Polarized(rho/(4.0*M_PI*r*r), eXC_1);
  eXC = eX + eC_0 + (eC_1 - eC_0)*fChachiyo2016(zeta);
  // vXCup, vXCdown = eXC + rho (d eXC/d rho) +- (1 -+ zeta) (d eXC/d zeta)
  // vXC = vX + vC_0 + (vC_1 - vC_0)*fChachiyo2016(zeta);
  vXCup = eXC + rho * decx + (1.0 - zeta) * decDzeta;
  vXCdown = eXC + rho * decx - (1.0 + zeta) * decDzeta;
}

Real exchangeCorrelationPotentialLDA(Real rho, Real r)
{
  Real eXC;
  /*
  if(rho < 1.0e-9)
    return 0.0;
  // return alpha2(std::pow(3.0/rho , 1.0/3.0), 0.0, 1.0, eXC);
  return alpha2(std::pow(3.0*r*r/rho , 1.0/3.0), 0.0, 1.0, eXC);
  //*/
  return chachiyo2016(rho/(4.0*M_PI*r*r), eXC);
}

void exchangeCorrelationPotentialLDA(std::vector<Real> &rho, std::vector<Real> &r_mesh, std::vector<Real> &vXC)
{
  for(int i=0; i<rho.size(); i++)
  {
    vXC[i] = exchangeCorrelationPotentialLDA(rho[i], r_mesh[i]);
  }
}

