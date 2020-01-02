/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#include "Real.hpp"
#include "Complex.hpp"
#include <vector>
#include "Matrix.hpp"
#include "Array3d.hpp"

#include "calculateDensities.hpp"

extern "C"
{
  void interp_(Real *r, Real *f, int *nr, Real *rs, Real *ps, Real *dps, int *deriv);
  void newint_(int *nr, Real *r, Real *f, Real *g, int *ip0);
}

void calculateDensities(LSMSSystemParameters &lsms, int iatom, int is, int ie, int nume, Complex energy, Complex dele1,
                        Matrix<Complex> &dos,Matrix<Complex> &dosck,Array3d<Complex> &green,
                        Array3d<Complex> &dipole,
                        AtomData &atom)
{
// from LSMS_19. zplanint after gettau
  Complex ede=energy*dele1;

  for(int isp=0; isp<lsms.n_spin_cant*lsms.n_spin_cant; isp++)
  {
    int isq=isp;
    if(lsms.n_spin_cant<2) isq=is;
    atom.dos_real(ie,isq)=std::imag(dos(isq,iatom));
    atom.evalsum[isq]+=std::imag(ede*dos(isq,iatom));
    atom.dosint[isq]+=std::imag(dos(isq,iatom)*dele1);
    atom.dosckint[isq]+=std::imag(dosck(isq,iatom)*dele1);
    if(lsms.relativity!=full)
    {
      for(int ir=0; ir<green.l_dim1(); ir++)
        atom.greenint(ir,isq)+=std::imag(green(ir,isq,iatom)*dele1);
    } else {
      // the fully relativistic densities at this point are actually
      // 4 pi r^2 rho(r), whereas the non relativistic densities are 4 pi rho(r)
      // here we correct for this:
      for(int ir=0; ir<green.l_dim1(); ir++)
      {
        Real r2inv=1.0/(atom.r_mesh[ir]*atom.r_mesh[ir]);
        atom.greenint(ir,isq)+=std::imag(green(ir,isq,iatom)*dele1)*r2inv;
      }
    }
  }
  for(int m=0; m<6; m++)
  {
    atom.dip[m]=std::imag(dipole(m,0,iatom)*dele1);
    if(lsms.n_spin_cant==2)
    {
      for(int isp=1; isp<4; isp++)
        atom.dip[m]+=std::imag(dipole(m,isp,iatom)*dele1*atom.evec[isp-1]);
      atom.dip[m]*=0.5;
    }
  }
  if(ie==nume-1)
  {
    for(int isp=0; isp<lsms.n_spin_cant*lsms.n_spin_cant; isp++)
    {
      int isq=isp;
      if(lsms.n_spin_cant<2) isq=is;
      atom.doslast[isq]=std::imag(dos(isq,iatom));
      atom.doscklast[isq]=std::imag(dosck(isq,iatom));
      if(lsms.relativity!=full)
      {
        for(int ir=0; ir<green.l_dim1(); ir++)
          atom.greenlast(ir,isq)=std::imag(green(ir,isq,iatom));
      } else {
        // the fully relativistic densities at this point are actually
      // 4 pi r^2 rho(r), whereas the non relativistic densities are 4 pi rho(r)
      // here we correct for this:
        for(int ir=0; ir<green.l_dim1(); ir++)
        {
          Real r2inv=1.0/(atom.r_mesh[ir]*atom.r_mesh[ir]);
          atom.greenlast(ir,isq)=std::imag(green(ir,isq,iatom));
        }
      }
    }
  }
}



// LSMS_1.9 getchg_c.f
void calculateChargeDensity(LSMSSystemParameters &lsms, AtomData &atom, Real edote,
                            Matrix<Real> &rhonew, Real &qvalmt, Real *qrms)
{
  Real rtmp[atom.jws+2];
  Real rhotmp[atom.jws+2];
  Real w1[atom.jws+3], w2[atom.jws+3];
  Real r_sph=atom.rInscribed;
  int ir_sph;

// need these constants to call fortran routines:
  int four=4;
  int five=5;
  int seven=7;
  Real dzero=0.0;

  if(lsms.mtasa>0) r_sph=atom.rws;

// construct rho: rho(ir,0) = up density; rho(ir,1) = down density

  rtmp[0]=0.0;
  for(int i=0; i<atom.jws; i++)
  {
    rhonew(i,0)=rhonew(i,1)=0.0;
    rtmp[i+1]=std::sqrt(atom.r_mesh[i]);
    if(atom.r_mesh[i] < r_sph) ir_sph = i+2;  // i+1?
  }
  if(ir_sph > atom.jws) ir_sph = atom.jws;

  if(lsms.n_spin_cant==2)
  {
// calculate valence charge density and store it in w1
    w1[0]=w2[0]=0.0;
    for(int ir=0; ir<atom.jws; ir++)
    {
      w1[ir+1]=-atom.r_mesh[ir]*atom.r_mesh[ir]*atom.greenint(ir,0)/M_PI;
      // calculate valence moment density by
      // projecting magnetization onto evec
      w2[ir+1]=-atom.r_mesh[ir]*atom.r_mesh[ir]*(atom.greenint(ir,1)*atom.evec[0]+
                                                 atom.greenint(ir,2)*atom.evec[1]+
                                                 atom.greenint(ir,3)*atom.evec[2])/M_PI;
      // protect against moment being larger than total charge because this
      // doesnt make sense and because it causes undefined exchange energy.
      // ! meis: do this only in the non relativistic case:
      if(lsms.relativity != full) // non relativistic or scalar relativistic calculation
      {
        if(w2[ir+1]>0.0 && w2[ir+1] >  w1[ir+1]) w2[ir+1] =  w1[ir+1];
        if(w2[ir+1]<0.0 && w2[ir+1] < -w1[ir+1]) w2[ir+1] = -w1[ir+1];
      }
    }

    for(int ir=0; ir<atom.jws; ir++)
    {
      rhonew(ir,0)=0.5*(w1[ir+1]+w2[ir+1]);
      rhonew(ir,1)=0.5*(w1[ir+1]-w2[ir+1]);
    }
  } else { // n_spin_cant != 2
    for(int ir=0; ir<atom.jws; ir++)
    {
      rhonew(ir,0)=-atom.r_mesh[ir]*atom.r_mesh[ir]*atom.greenint(ir,0)/M_PI;
      rhonew(ir,1)=-atom.r_mesh[ir]*atom.r_mesh[ir]*atom.greenint(ir,1)/M_PI;
    }
  }

#if 0
  FILE *of;
  of=fopen("rhonew.dat","w");
  fprintf(of,"#r  rho m  rhoup rhodown\n");
  for(int ir=0; ir<atom.jws; ir++)
  {
    fprintf(of,"%lf  %lf %lf  %lf %lf\n",atom.r_mesh[ir],w1[ir+1],w2[ir+1],rhonew(ir,0),rhonew(ir,1));
  }
  fclose(of);
  exit(1);
#endif

  // if forceZeroMoment != force the up and down densities to be identical
  if(atom.forceZeroMoment && (lsms.n_spin_pola!=1))
  {
    Real diff = 0.0;
    for(int ir=0; ir<atom.jws; ir++)
    {
      diff += (rhonew(ir,0)-rhonew(ir,1))*(rhonew(ir,0)-rhonew(ir,1));
      // rhonew(ir,1) = rhonew(ir,0) = 0.5*(rhonew(ir,0)+rhonew(ir,1));
      // atom.corden(ir,1) = atom.corden(ir,0) = 0.5*(atom.corden(ir,0) + atom.corden(ir,1));
      // atom.semcor(ir,1) = atom.semcor(ir,0) = 0.5*(atom.semcor(ir,0) + atom.semcor(ir,1));
    }

    // printf("forceZeroAtom: diff(rhonew up, down)=%lg\n",std::sqrt(diff));

    // atom.averageSpins();
  }

  // calculate valence charge
  qvalmt=0.0;
  for(int is=0; is<lsms.n_spin_pola; is++)
  {
    Real dummy, w1_r_sph;
    Real sqrt_r_sph = std::sqrt(r_sph);
    int ir_sph_p1 = ir_sph+1;
    int f = 0;

    for(int ir=0; ir<atom.jws; ir++)
      rhotmp[ir+1] = 2.0*rhonew(ir,is) / (atom.r_mesh[ir]*atom.r_mesh[ir]);

    interp_(&rtmp[1],&rhotmp[1],&four,&dzero,&rhotmp[0],&dummy,&f);
    newint_(&ir_sph_p1,rtmp,rhotmp,w1,&five);
    interp_(&rtmp[0],w1,&ir_sph_p1,&sqrt_r_sph,&w1_r_sph,&dummy,&f);

    qvalmt+=w1_r_sph;
  }

  if(lsms.global.iprint>=0) printf("Total valence charge = %25.20f\n",qvalmt);

  // add velence electron density and core and semi-core densities
  if (lsms.n_spin_pola != 2)
  {
    for(int ir=0; ir<atom.jws; ir++)
      w1[ir+1] = w2[ir+1] = atom.corden(ir,0) + atom.semcor(ir,0);
  }
  else {
    for (int ir=0; ir<atom.jws; ir++)
    {
      w1[ir+1] = atom.corden(ir,0) + atom.semcor(ir,0);
      w2[ir+1] = atom.corden(ir,1) + atom.semcor(ir,1);
    }
  }
  for (int is=0; is<lsms.n_spin_pola; is++)
  {
    Real fac = (1.0-2.0*is)*edote;
    for(int ir=0; ir<atom.jws; ir++) 
      rhonew(ir,is)+=0.5*(w1[ir+1]+w2[ir+1]+fac*(w1[ir+1]-w2[ir+1]));
  }

  // if(atom.forceZeroMoment && (lsms.n_spin_pola!=1))
  // {
  //   for(int ir=0; ir<atom.jws; ir++)
  //   { 
  //     rhonew(ir,1) = rhonew(ir,0) = 0.5*(rhonew(ir,0)+rhonew(ir,1)); 
  //     // atom.corden(ir,1) = atom.corden(ir,0) = 0.5*(atom.corden(ir,0) + atom.corden(ir,1));
  //     // atom.semcor(ir,1) = atom.semcor(ir,0) = 0.5*(atom.semcor(ir,0) + atom.semcor(ir,1));
  //   }
  // }


// calculate the rms between old and new charge densities
  for(int is=0; is<lsms.n_spin_pola; is++)
  {
    Real dummy;
    Real sqrt_r_sph = std::sqrt(r_sph);
    int ir_sph_p1 = ir_sph + 1;
    int f = 0;
    int isold = is;
    if (atom.spinFlipped)
      isold = lsms.n_spin_pola - is - 1;       // 0->1 and 1->0 and if n_spin_pola==1: 0->0 !
    for(int ir=0; ir<atom.jws; ir++)
      w2[ir+1] = 2.0 * std::pow((rhonew(ir,is)-atom.rhotot(ir,is)),2) / std::pow(atom.r_mesh[ir],4);
    interp_(&atom.r_mesh[0],&w2[1],&four,&dzero,&w2[0],&dummy,&f);
    newint_(&ir_sph_p1,&rtmp[0],&w2[0],w1,&seven);
    interp_(&rtmp[0],&w1[0],&ir_sph_p1,&sqrt_r_sph,&qrms[is],&dummy,&f);
    qrms[is] = std::sqrt(qrms[is] / (4.0*M_PI*atom.omegaMT));
    // qrms[is] = qrms[is] / (4.0*M_PI*atom.omegaMT);
  }
}

void calculateAllLocalChargeDensities(LSMSSystemParameters &lsms, LocalTypeInfo &local)
{
#pragma omp parallel for default(none) shared(lsms,local)
  for(int i=0; i<local.num_local; i++)
  {
    Real edote = local.atom[i].evec[0]*local.atom[i].evecNew[0]+
                 local.atom[i].evec[1]*local.atom[i].evecNew[1]+
                 local.atom[i].evec[2]*local.atom[i].evecNew[2];
    calculateChargeDensity(lsms, local.atom[i], edote,
                           local.atom[i].rhoNew, local.atom[i].qvalmt, &local.atom[i].qrms[0]); // &local_qrms[2*i]);
    // checkCharge(lsms, local.atom[i]);
  }

  int n=0;
  local.qrms[0]=local.qrms[1]=0.0;
  for(int i=0; i<local.num_local; i++)
  {
    n+=local.n_per_type[i];
    local.qrms[0]+=local.atom[i].qrms[0]*Real(local.n_per_type[i]);
    local.qrms[1]+=local.atom[i].qrms[1]*Real(local.n_per_type[i]);
  }
  local.qrms[0] = local.qrms[0]/Real(n);
  local.qrms[1] = local.qrms[1]/Real(n);
}

void checkAllLocalCharges(LSMSSystemParameters &lsms, LocalTypeInfo &local)
{
#pragma omp parallel for default(none) shared(lsms,local)
  for(int i=0; i<local.num_local; i++)
  {
    checkCharge(lsms, local.atom[i]);
  }
}
