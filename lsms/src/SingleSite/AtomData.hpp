/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
#ifndef LSMS_ATOMDATA_H
#define LSMS_ATOMDATA_H

#include <vector>
#include <algorithm>

#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"
#include "Array3d.hpp"
#include "VORPOL/VORPOL.hpp"

extern "C"
{
void spin_trafo_(Real *evec, Complex * u, Complex *ud);
}


#ifdef BUILDKKRMATRIX_GPU
void * allocateDConst(void);
void freeDConst(void *);
#endif


class AtomData {
public:

  AtomData() {reset();}

  ~AtomData() {}

  void reset(void)
  {
    b_con[0] = b_con[1] = b_con[2] = 0.0;
    for(int i=0; i<9; i++)
      b_basis[i] = 0.0;
    b_basis[0] = b_basis[4] = b_basis[8] = 1.0;
    mConstraint = -1.0;
    spinFlipped = false;
  }

  void reset_b_basis(void)
  {
    for(int i=0; i<9; i++) b_basis[i] = 0.0;
    b_basis[0] = b_basis[4] = b_basis[8] = 1.0;
  }

  void get_b_basis(void)
  {
    const Real tol=1.0e-8;
    Real norm2=evec[0]*evec[0]+evec[1]*evec[1]+evec[2]*evec[2];
    Real norm=std::sqrt(norm2);
    evec[0]=evec[0]/norm; evec[1]=evec[1]/norm; evec[2]=evec[2]/norm;
    Real cost=evec[2];
    Real sint = std::sqrt(1.0-cost*cost);
    Real cosp=1.0;
    Real sinp=0.0;
    if(std::abs(sint)>tol)
    {
      cosp=evec[0]/sint; sinp=evec[1]/sint;
    }

    b_basis[6+0]=evec[0]; b_basis[6+1]=evec[1]; b_basis[6+2]=evec[2];

    b_basis[0]=cost*cosp; b_basis[1]=cost*sinp; b_basis[2]=-sint;
    norm2=b_basis[0]*b_basis[0]+b_basis[1]*b_basis[1]+b_basis[2]*b_basis[2];
    norm=std::sqrt(norm2);
    b_basis[0]=b_basis[0]/norm; b_basis[1]=b_basis[1]/norm; b_basis[2]=b_basis[2]/norm;


    b_basis[3+0] = b_basis[6+1]*b_basis[0+2] - b_basis[6+2]*b_basis[0+1];
    b_basis[3+1] = -(b_basis[6+0]*b_basis[0+2] - b_basis[6+2]*b_basis[0+0]);
    b_basis[3+2] = b_basis[6+0]*b_basis[0+1] - b_basis[6+1]*b_basis[0+0];
  }

  void newConstraint(void)
  {
    Real b_con_mag=0.0;
    Real b_con_g[3];
    Real b_con_o[3];
    Real b_con_in[3];
    Real moment_dir[3];

    for(int i=0; i<3; i++)
    {
      b_con_g[i]=b_con[0]*b_basis[0+i] + b_con[1]*b_basis[3+i] + b_con[2]*b_basis[6+i];
      b_con_o[i]=b_con_g[i];
      b_con_in[i]=b_con[i];
      b_con_mag+=b_con_g[i]*b_con_g[i];
    }
    b_con_mag=std::sqrt(b_con_mag);
    // Real moment_mag=std::sqrt(moment[0]*moment[0]+moment[1]*moment[1]+moment[2]*moment[2]);
    // moment_dir[0]=moment[0]/moment_mag; moment_dir[1]=moment[1]/moment_mag;
    // moment_dir[2]=moment[2]/moment_mag;
    Real moment_mag=std::sqrt(evecOut[0]*evecOut[0]+evecOut[1]*evecOut[1]+evecOut[2]*evecOut[2]);
    moment_dir[0]=evecOut[0]/moment_mag; moment_dir[1]=evecOut[1]/moment_mag;
    moment_dir[2]=evecOut[2]/moment_mag;

    Real mproj=moment_dir[0]*evec[0]+moment_dir[1]*evec[1]+moment_dir[2]*evec[2];
    Real bproj=moment_dir[0]*b_con_g[0]+moment_dir[1]*b_con_g[1]+moment_dir[2]*b_con_g[2];

    //Real bb1=0.0;
    for(int i=0; i<3; i++)
    {
      b_con_g[i]=b_con_g[i]+evec[i]-(mproj+bproj)*moment_dir[i];
      //bb1+=b_con_g[i]*b_con_g[i];
    }
    //bb1=std::sqrt(bb1);

//     ================================================================
//     project components in b_basis frame ...................
//     ================================================================

    for(int j=0; j<3; j++)
    {
      b_con[j]=0.0;
      for(int i=0; i<3; i++)
        b_con[j]+=b_con_g[i]*b_basis[3*j+i];
    }

// meis: as a Test set b_con[2]=0.0
//    b_con[2]=0.0;

  }

  void resizePotential(int npts)
  {
    int n = npts;

    vr.resize(n,2);
    vr = 0.0;
    vSpinShift = 0.0;

    rhotot.resize(n,2);
    rhotot = 0.0;

    corden.resize(n,2);
    corden = 0.0;

    semcor.resize(n,2);
    r_mesh.resize(n);
    x_mesh.resize(n);

    vrNew.resize(n,2);
    vrNew = 0.0;

    rhoNew.resize(n,2);
    rhoNew = 0.0;

    dos_real.resize(40,4);
    greenint.resize(n,4);
    greenlast.resize(n,4);

    // Check if they should be put here...
    exchangeCorrelationPotential.resize(n,2);
    exchangeCorrelationPotential = 0.0;
    exchangeCorrelationEnergy.resize(n,2);
    exchangeCorrelationEnergy = 0.0;

  }


  void resizeCore(int ncs)
  {
    ec.resize(ncs,2);
    nc.resize(ncs,2);
    lc.resize(ncs,2);
    kc.resize(ncs,2);
  }

  void changeNspin(int nspinNew)
  {
    if(nspin == nspinNew) return;
    if(nspinNew == 2) // extend from non spin polarized to spin polarized
    {
      xvalws[1] = xvalws[0];
      xvalwsNew[1] = xvalwsNew[0];
      xvalmt[1] = xvalmt[0];
      for(int i=0; i<vr.l_dim(); i++)
	vr(i,1) = vr(i,0);
      for(int i=0; i<rhotot.l_dim(); i++)
	rhotot(i,1) = rhotot(i,0);
      for(int i=0; i<vrNew.l_dim(); i++)
	vrNew(i,1) = vrNew(i,0);
      for(int i=0; i<rhoNew.l_dim(); i++)
	rhoNew(i,1) = rhoNew(i,0);
      
      for(int i=0; i<exchangeCorrelationPotential.l_dim(); i++)
	exchangeCorrelationPotential(i,1) = exchangeCorrelationPotential(i,0);
      for(int i=0; i<exchangeCorrelationEnergy.l_dim(); i++)
	exchangeCorrelationEnergy(i,1) = exchangeCorrelationEnergy(i,0);

      exchangeCorrelationV[1] = exchangeCorrelationV[0];

      for(int i=0; i<ec.l_dim(); i++)
	ec(i,1) = ec(i,0);
      for(int i=0; i<nc.l_dim(); i++)
	nc(i,1) = nc(i,0);
      for(int i=0; i<lc.l_dim(); i++)
	lc(i,1) = lc(i,0);
      for(int i=0; i<kc.l_dim(); i++)
	kc(i,1) = kc(i,0);

      ecorv[1] = ecorv[0]; esemv[1] = esemv[0];
      for(int i=0; i<corden.l_dim(); i++)
	corden(i,1) = corden(i,0);
      for(int i=0; i<semcor.l_dim(); i++)
	semcor(i,1) = semcor(i,0);
      
      nspin = 2;
    }
    if(nspinNew == 1) // extend from non spin polarized to spin polarized
    {
      xvalws[0] = 0.5*(xvalws[0]+xvalws[1]);
      xvalwsNew[0] = 0.5*(xvalwsNew[0]+xvalwsNew[1]);
      xvalmt[0] = 0.5*(xvalmt[0]+xvalmt[1]);
      for(int i=0; i<vr.l_dim(); i++)
	vr(i,0) = 0.5*(vr(i,0)+vr(i,1));;
      for(int i=0; i<rhotot.l_dim(); i++)
	rhotot(i,0) = 0.5*(rhotot(i,0)+rhotot(i,1));
      for(int i=0; i<vrNew.l_dim(); i++)
	vrNew(i,0) = 0.5*(vrNew(i,0)+vrNew(i,1));
      for(int i=0; i<rhoNew.l_dim(); i++)
	rhoNew(i,0) = 0.5*(rhoNew(i,0)+rhoNew(i,1));
      
      for(int i=0; i<exchangeCorrelationPotential.l_dim(); i++)
	exchangeCorrelationPotential(i,0) = 0.5*(exchangeCorrelationPotential(i,0)+exchangeCorrelationPotential(i,1));
      for(int i=0; i<exchangeCorrelationEnergy.l_dim(); i++)
	exchangeCorrelationEnergy(i,0) = exchangeCorrelationEnergy(i,0)+exchangeCorrelationEnergy(i,1);

      exchangeCorrelationV[0] = 0.5*(exchangeCorrelationV[0]+exchangeCorrelationV[1]);

      /*
      for(int i=0; i<ec.l_dim(); i++)
	ec(i,1) = ec(i,0);
      for(int i=0; i<nc.l_dim(); i++)
	nc(i,1) = nc(i,0);
      for(int i=0; i<lc.l_dim(); i++)
	lc(i,1) = lc(i,0);
      for(int i=0; i<kc.l_dim(); i++)
	kc(i,1) = kc(i,0);
      */

      ecorv[0] = ecorv[0]+ecorv[1]; esemv[0] = esemv[0]+esemv[1];
      for(int i=0; i<corden.l_dim(); i++)
	corden(i,0) = 0.5*(corden(i,0)+corden(i,1));
      for(int i=0; i<semcor.l_dim(); i++)
	semcor(i,0) = 0.5*(semcor(i,0)+semcor(i,1));
      
      nspin = 1;
    }
  }

  void averageSpins(void)
  {
    if(nspin == 1) return;
    
    xvalws[0] = 0.5*(xvalws[0]+xvalws[1]);
    xvalws[1] = xvalws[0];
    xvalwsNew[0] = 0.5*(xvalwsNew[0]+xvalwsNew[1]);
    xvalwsNew[1] = xvalwsNew[0];
    xvalmt[0] = 0.5*(xvalmt[0]+xvalmt[1]);
    xvalmt[1] = xvalmt[0];
    for(int i=0; i<vr.l_dim(); i++)
    {
      vr(i,0) = 0.5*(vr(i,0)+vr(i,1));
      vr(i,1) = vr(i,0);
    }
    for(int i=0; i<rhotot.l_dim(); i++)
    {
      rhotot(i,0) = 0.5*(rhotot(i,0)+rhotot(i,1));
      rhotot(i,1) = rhotot(i,0);
    }
    for(int i=0; i<vrNew.l_dim(); i++)
    {
      vrNew(i,0) = 0.5*(vrNew(i,0)+vrNew(i,1));
      vrNew(i,1) = vrNew(i,0);
    }
    for(int i=0; i<rhoNew.l_dim(); i++)
    {
      rhoNew(i,0) = 0.5*(rhoNew(i,0)+rhoNew(i,1));
      rhoNew(i,1) = rhoNew(i,0);
    }
      
    for(int i=0; i<exchangeCorrelationPotential.l_dim(); i++)
    {
      exchangeCorrelationPotential(i,0) = 0.5*(exchangeCorrelationPotential(i,0)+exchangeCorrelationPotential(i,1));
      exchangeCorrelationPotential(i,1) = exchangeCorrelationPotential(i,0);
    }
       
    exchangeCorrelationV[0] = 0.5*(exchangeCorrelationV[0]+exchangeCorrelationV[1]);
    exchangeCorrelationV[1] = exchangeCorrelationV[0];

    for(int i=0; i<corden.l_dim(); i++)
    {
      corden(i,0) = 0.5*(corden(i,0)+corden(i,1));
      corden(i,1) = corden(i,0);
    }
    for(int i=0; i<semcor.l_dim(); i++)
    {
      semcor(i,0) = 0.5*(semcor(i,0)+semcor(i,1));
      semcor(i,1) = semcor(i,0);
    }
  }
  
  AtomData &operator=(const AtomData &a)
  {
    jmt = a.jmt;
    jws = a.jws;
    xstart = a.xstart;
    rmt = a.rmt;
    h = a.h;
    r_mesh = a.r_mesh;
    x_mesh = a.x_mesh;

    alat = a.alat;
    efermi = a.efermi;
    vdif = a.vdif;
    ztotss = a.ztotss;
    zcorss = a.zcorss;
    zsemss = a.zsemss;
    zvalss = a.zvalss;

    nspin = a.nspin;
    forceZeroMoment = a.forceZeroMoment;
    numc = a.numc;
    spinFlipped = a.spinFlipped;
    localEnergy = a.localEnergy;

    evec[0] = a.evec[0];
    evec[1] = a.evec[1];
    evec[2] = a.evec[2];
    evecNew[0] = a.evecNew[0];
    evecNew[1] = a.evecNew[1];
    evecNew[2] = a.evecNew[2];
    evecOut[0] = a.evecOut[0];
    evecOut[1] = a.evecOut[1];
    evecOut[2] = a.evecOut[2];
    xvalws[0] = a.xvalws[0];
    xvalws[1] = a.xvalws[1];
    xvalwsNew[0] = a.xvalwsNew[0];
    xvalwsNew[1] = a.xvalwsNew[1];
    xvalmt[0] = a.xvalmt[0];
    xvalmt[1] = a.xvalmt[1];
    qtotws = a.qtotws;
    mtotws = a.mtotws;
    qtotmt = a.qtotmt;
    mtotmt = a.mtotmt;
    qvalws = a.qvalws;
    mvalws = a.mvalws;
    qvalmt = a.qvalmt;
    mvalmt = a.mvalmt;

    qInt = a.qInt;
    mInt = a.mInt;
    rhoInt = a.rhoInt;
    mIntComponent[0] = a.mIntComponent[0];
    mIntComponent[1] = a.mIntComponent[1];
    mIntComponent[2] = a.mIntComponent[2];

    for(int i=0; i<80; i++) header[i] = a.header[i];

    vr = a.vr;
    vSpinShift=a.vSpinShift;
    rhotot = a.rhotot;
    
    exchangeCorrelationPotential = a.exchangeCorrelationPotential;
    exchangeCorrelationEnergy = a.exchangeCorrelationEnergy;
    exchangeCorrelationV[0] = a.exchangeCorrelationV[0];
    exchangeCorrelationV[1] = a.exchangeCorrelationV[1];
    exchangeCorrelationE = a.exchangeCorrelationE;

    ec = a.ec;
    nc = a.nc;
    lc = a.lc;
    kc = a.kc;

    ecorv[0] = a.ecorv[0];
    ecorv[1] = a.ecorv[1];
    esemv[0] = a.esemv[0];
    esemv[1] = a.esemv[1];

    corden = a.corden;
    semcor = a.semcor;
    qcpsc_mt = a.qcpsc_mt;
    qcpsc_ws = a.qcpsc_ws;
    mcpsc_mt = a.mcpsc_mt;
    mcpsc_ws = a.mcpsc_ws;

    b_con[0] = a.b_con[0];
    b_con[1] = a.b_con[1];
    b_con[2] = a.b_con[2];

    for(int i=0; i<9; i++) b_basis[i] = a.b_basis[i];

    vrms[0]=a.vrms[0]; vrms[1]=a.vrms[1]; qrms[0]=a.qrms[0]; qrms[1]=a.qrms[1];

    return *this;
  }


  void generateRadialMesh(void)
  {
    int N = std::max((int)r_mesh.size(), jws);
    N = std::max(N, jmt);
    if (N != r_mesh.size()) r_mesh.resize(N);
    Real xmt = std::log(rmt);
    h = (xmt-xstart) / (jmt-1);
    for(int j=0; j<N; j++)
    {
      x_mesh[j] = xstart + (Real)j*h;
      r_mesh[j] = std::exp(x_mesh[j]);
      //printf("j = %5d x_mesh = %45.40e r_mesh = %45.40e\n", j, x_mesh[j], r_mesh[j]);
    }
    generateNewMesh = false;
  }


  void setEvec(Real x, Real y, Real z)
  {
    evec[0] = x;
    evec[1] = y;
    evec[2] = z;
    spin_trafo_(evec, ubr, ubrd);
  }


// Local Interaction Zone
  int numLIZ;
  std::vector<int> LIZGlobalIdx, LIZStoreIdx, LIZlmax;
  int nrmat;                          // sum (LIZlmax+1)^2
  std::vector<Real> LIZDist;
  Matrix<Real> LIZPos;

// Mesh Data:
  int jmt,jws;
  Real xstart,rmt,h;
  Real rInscribed; // LSMS_1.9: rins
  std::vector<Real> r_mesh, x_mesh;
  bool generateNewMesh;

// General Data
  char header[80];
  int lmax, kkrsz;
  Real alat, efermi;
  Real ztotss, zcorss, zsemss, zvalss;
  Real vdif, vdifNew;
  Real evec[3], evecNew[3], evecOut[3];
  Complex ubr[4], ubrd[4];            // Spin transformation matrices
  Matrix<Complex> dmat, dmatp;        // Spin rotation matrices for the relativistic case
  Complex wx[4], wy[4], wz[4];
  Real xvalws[2], xvalwsNew[2];
  Real xvalmt[2];
  Real qtotws, mtotws;
  Real qtotmt, mtotmt;
  Real qvalws, mvalws;
  Real qvalmt, mvalmt;
  Real qInt, mInt, rhoInt;            // Interstitial charge, moment and charge density
  Real mIntComponent[3];              // Interstitial moment components
  int nspin;                          // Number of spin direction-related settings
                                      //  (determines n_spin_cant & n_spin_pola)
  int forceZeroMoment;                // if != 0, average spin up and spin down densities
                                      // to force zero magnetic moment
  int numc;                           // Number of core states
  bool spinFlipped;                   // Flag for antiferromagnetic condition

// local energy
  Real localEnergy;
  Real localMadelungEnergy;
  
// Alloy Class
  int alloy_class;

// Volumes:
  Real omegaMT;                       // Muffin-Tin volume
  Real omegaWS;                       // Wigner-Seitz volume
  Real rws;                           // Wigner-Seitz radius

// omegaInt - interstitial volume is in voronoi.omegaInt

// Madelung matrix
  std::vector<Real> madelungMatrix;

// Potential and charge density
  Real vSpinShift; // relativ shift of the spin up and spin down potentials
                   // for use in WL-LSMS. Applied using the PotentialShifter class
  Matrix<Real> vr, rhotot;

// Storage for newly calculated potential and chage densities before mixing
  Matrix<Real> vrNew,rhoNew;

// Exchange-correlation parameters
  Matrix<Real> exchangeCorrelationPotential;     // Exchange-correlation potential
  Matrix<Real> exchangeCorrelationEnergy;        // Exchange-correlation energy
  Real exchangeCorrelationE;                     // Exchange-correlation energy
  Real exchangeCorrelationV[2];                  // Exchange-correlation potential for spin up/down

// Core state info
  Matrix<Real> ec;
  Matrix<int> nc, lc, kc;
  Real ecorv[2], esemv[2];
  Matrix<Real> corden, semcor;
  Real qcpsc_mt, qcpsc_ws, mcpsc_mt, mcpsc_ws;

// Constraint data
  enum {None, Direction, Moment} constraintType;
  Real b_con[3];
  Real b_basis[9];
  Real mConstraint;

// vector for the energy points in eGroup
  std::vector<Matrix<Complex> > pmat_m;

  VoronoiPolyhedra voronoi;

// local densities
  Matrix<Real> dos_real;
  Real doslast[4];
  Real doscklast[4];
  Real evalsum[4];
  Real dosint[4];
  Real dosckint[4];
  Matrix<Real> greenint;
  Matrix<Real> greenlast;
  Real dip[6];

// rms changes between iterations:
  Real vrms[2];
  Real qrms[2];

  void resetLocalDensities(void)
  {
    dos_real=0.0;
    greenint=0.0;
    greenlast=0.0;
    doslast[0]=doslast[1]=doslast[2]=doslast[3]=0.0;
    doscklast[0]=doscklast[1]=doscklast[2]=doscklast[3]=0.0;
    evalsum[0]=evalsum[1]=evalsum[2]=evalsum[3]=0.0;
    dosint[0]=dosint[1]=dosint[2]=dosint[3]=0.0;
    dosckint[0]=dosckint[1]=dosckint[2]=dosckint[3]=0.0;
    dip[0]=dip[1]=dip[2]=dip[3]=dip[4]=dip[5]=0.0;
  }
};

#endif
