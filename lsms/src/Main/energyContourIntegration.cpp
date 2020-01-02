/* -*- c-file-style: "bsd"; c-basic-offset: 2; indent-tabs-mode: nil -*- */
// this replaces zplanint from LSMS_1.9

#include <vector>
#include <mpi.h>
#include <complex>
#include "Complex.hpp"

#include "Communication/LSMSCommunication.hpp"
#include "SingleSite/SingleSiteScattering.hpp"
#include "MultipleScattering/MultipleScattering.hpp"
#include "EnergyContourIntegration.hpp"
#include "Misc/Coeficients.hpp"
#include "calculateDensities.hpp"
// #include <omp.h>
#ifdef USE_NVTX
#include <nvToolsExt.h>
#endif

#if defined(ACCELERATOR_CULA) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C)
void copyTmatStoreToDevice(LocalTypeInfo &local);
#ifdef BUILDKKRMATRIX_GPU
#include "Accelerator/buildKKRMatrix_gpu.hpp"
extern std::vector<DeviceConstants> deviceConstants;
// extern std::vector<void *> deviceConstants;
// extern void * deviceStorage;
#endif
//#if defined(ACCELERATOR_CULA) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C)
#include "Accelerator/DeviceStorage.hpp"
extern DeviceStorage *deviceStorage;
#endif
void solveSingleScatterers(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                           std::vector<Matrix<Real> > &vr, Complex energy,
                           std::vector<NonRelativisticSingleScattererSolution> &solution,int iie);
void solveSingleScatterers(LSMSSystemParameters &lsms, LocalTypeInfo &local,
                           std::vector<Matrix<Real> > &vr, Complex energy,
                           std::vector<RelativisticSingleScattererSolution> &solution,int iie);

void rotateToGlobal(AtomData &atom, Matrix<Complex> &dos, Matrix<Complex> &dosck,
                    Matrix<Complex> &dos_orb, Matrix<Complex> &dosck_rob,
                    Array3d<Complex> &green, Array3d<Complex> &dens_orb, int i);

extern "C"
{
//     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
void constraint_(int *jmt,Real *rmt,int *n_spin_pola,
                 Real*vr,Real *r_mesh,Real *pi4,
                 Real *evec,Real *evec_r,Real *b_con,Real *b_basis,
                 int *i_vdif,Real *h_app_para_mag,Real *h_app_perp_mag,
                 int *iprpts,
                 int *iprint,char *istop,int len_sitop);
//     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
void u_sigma_u_(Complex *ubr,Complex *ubrd,
               Complex *wx,Complex *wy,Complex *wz);
//     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
void green_function_(int *mtasa,int *n_spin_pola,int *n_spin_cant,
                     int *lmax,int *kkrsz,Complex *wx, Complex *wy, Complex *wz,
                     Real *rins,Real *r_sph,Real *r_mesh,int *jmt,int *jws,
                     Complex *pnrel,Complex *tau00_l,Complex *matom,Complex *zlr,Complex *jlr,
                     int *nprpts, int* nplmax,
                     int *ngaussr,
                     Real *cgnt, int *lmax_cg,
                     Complex *dos,Complex *dosck,Complex *green,Complex *dipole,
                     int *ncrit,Real *grwylm,Real *gwwylm,Complex *wylm,
                     int *iprint,char *istop,int len_sitop);

void green_function_rel_(int *mtasa,
                         int *lmax, int *kkrsz,
                         Complex *wx, Complex *wy, Complex *wz,
                         Real *rins, Real *r_mesh, int *jmt, int *jws, Real *h,
                         Complex *pnrel, Complex *tau00_l, Complex *matom,
                         Complex *gz, Complex *fz,
                         Complex *gj, Complex *fj,
                         int *nuz, int *indz,
                         int *nprpts,
                         int *ngaussr, Real *cgnt, int *lmax_cg,
                         Complex *dos, Complex *dosck, Complex *green, Complex *dipole,
                         Complex *dos_orb, Complex *dosck_orb, Complex *dens_orb,
                         int *iprint, char *istop, int len_istop);
//     ================================================================
void clebsch_(void);
void gfill_(int *iplmax);
void gafill_(int *iplmax);
void matrot1_(Real * rGlobal, Real *evec_r, int *lmax,
              Complex *dmat, Complex *dmatp);
}

void buildEnergyContour(int igrid,Real ebot,Real etop,Real eibot, Real eitop,
                        std::vector<Complex> &egrd, std::vector<Complex> &dele1, int npts, int &nume,
                        int iprint, char *istop)
{
  Real pi=2.0*std::asin(1.0);
  int ipepts;
  if(iprint>=0) printf("Energy Contour Parameters: grid=%d npts=%d, ebot=%lf etop=%lf eibot=%lf eitop=%lf\n",
                      igrid,npts,ebot,etop,eibot,eitop);
  switch (igrid)
  {
  case 0: // single energy point
    egrd.resize(1); dele1.resize(1);
    egrd[0]=std::complex<Real>(ebot,eibot); dele1[0]=0.0; nume=1;
    break;
  case 2: // Gaussian Contour
    egrd.resize(npts+2); dele1.resize(npts+2);
    ipepts=egrd.size();
    congauss_(&ebot,&etop,&eibot,&egrd[0],&dele1[0],&npts,&nume,&pi,&ipepts,&iprint,istop,32);
    break;
  default:
    fprintf(stderr,"Unknown energy grid type %d in 'buildEnergyContour'\n",igrid);
    exit(1);
  }
}

void energyContourIntegration(LSMSCommunication &comm,LSMSSystemParameters &lsms, LocalTypeInfo &local)
{
  double timeEnergyContourIntegration_1=MPI_Wtime();

  if(lsms.global.iprint>=0) printf("** Energy Contour Integration **\n");

// calculate coefficients and matrices for spherical relativistic calculations
  if(lsms.relativity==full)
  {
    clebsch_();
    gfill_(&lsms.maxlmax);
    gafill_(&lsms.maxlmax);
  }

// energy grid info
  std::vector<Complex> egrd,dele1;
  int nume;
// constrained potentials:
  std::vector<Matrix<Real > > vr_con;
  Matrix<Real> evec_r;

  int i_vdif=0;

  Real pi4=4.0*2.0*std::asin(1.0);

  vr_con.resize(local.num_local);
  evec_r.resize(3,local.num_local);
#pragma omp parallel for default(none) shared(local,lsms,vr_con,evec_r)
  for(int i=0; i<local.num_local; i++)
  {
    Real pi4=4.0*2.0*std::asin(1.0);
    int i_vdif=0;
//     ================================================================
//     set up the spin space stransformation matrix....................
//     ================================================================
// check evec
    Real evec_norm=std::sqrt(local.atom[i].evec[0]*local.atom[i].evec[0]
                             +local.atom[i].evec[1]*local.atom[i].evec[1]
                             +local.atom[i].evec[2]*local.atom[i].evec[2]);
    if(std::abs(evec_norm-1.0)>1.0e-5) printf("|atom[%d].evec|=%lf\n",i,evec_norm);
//
    spin_trafo_(&local.atom[i].evec[0],&local.atom[i].ubr[0],&local.atom[i].ubrd[0]);
//     ================================================================
//     set up Constraint ..............................................
//     copy vr into vr_con which contains the B-field constraint.......
//     calls to gettau_c etc. require vr_con...........................
//     ================================================================
    vr_con[i]=local.atom[i].vr;
    if(lsms.n_spin_cant==2)
    {
      Real h_app_para_mag=0.0;
      Real h_app_perp_mag=0.0;
      // int iprpts=local.atom[i].r_mesh.size();
      int iprpts=vr_con[i].l_dim();
      Real rmt = local.atom[i].rmt;
      if(lsms.mtasa>0) rmt = local.atom[i].rws;
      int jmt = local.atom[i].jmt;
      if(lsms.mtasa==1) jmt = local.atom[i].jws;
// here I leave out the i_vdif<0 case!
      constraint_(&jmt,&rmt,&lsms.n_spin_pola,
                  &(vr_con[i])(0,0),&local.atom[i].r_mesh[0],&pi4,
                  &local.atom[i].evec[0],&evec_r(0,i),local.atom[i].b_con,
                  local.atom[i].b_basis,&i_vdif,&h_app_para_mag,&h_app_perp_mag,
                  &iprpts,
                  &lsms.global.iprint,lsms.global.istop,32);
      if(lsms.relativity != full)
      {
        spin_trafo_(&evec_r(0,i),&local.atom[i].ubr[0],&local.atom[i].ubrd[0]);
      } else { //. relativistic
        int matrot_size=2*(local.atom[i].lmax+1)*(local.atom[i].lmax+1);
        local.atom[i].dmat.resize(matrot_size,matrot_size);
        local.atom[i].dmatp.resize(matrot_size,matrot_size);
        Real rGlobal[3];
        rGlobal[0]=0.0;
        rGlobal[1]=0.0;
        rGlobal[2]=1.0;

        matrot1_(rGlobal,&evec_r(0,i),&local.atom[i].lmax,
                 &local.atom[i].dmat(0,0),&local.atom[i].dmatp(0,0));
      }
    } else { // n_spin_cant != 2 i.e. non spin polarized
// call zcopy(4,u,1,ubr,1)
// call zcopy(4,ud,1,ubrd,1)
    }
    u_sigma_u_(&local.atom[i].ubr[0],&local.atom[i].ubrd[0],
               &local.atom[i].wx[0],&local.atom[i].wy[0],&local.atom[i].wz[0]);
  }
  // Real e_top;
  // e_top=lsms.energyContour.etop;
  // if(lsms.energyContour.etop==0.0) etop=lsms.chempot;
  buildEnergyContour(lsms.energyContour.grid, lsms.energyContour.ebot, lsms.chempot,
                     lsms.energyContour.eibot, lsms.energyContour.eitop, egrd, dele1,
                     lsms.energyContour.npts, nume, lsms.global.iprint, lsms.global.istop);

  for(int i=0; i<local.num_local; i++)
  {
    if(local.atom[i].dos_real.l_dim()<nume) 
      local.atom[i].dos_real.resize(nume,4);
  }

  std::vector<std::vector<NonRelativisticSingleScattererSolution> >solutionNonRel;
  std::vector<std::vector<RelativisticSingleScattererSolution> >solutionRel;

  if(lsms.relativity!=full)
  {
    solutionNonRel.resize(lsms.energyContour.groupSize());
    for(int ie=0; ie<lsms.energyContour.groupSize(); ie++) solutionNonRel[ie].resize(local.num_local);
  } else {
    solutionRel.resize(lsms.energyContour.groupSize());
    for(int ie=0; ie<lsms.energyContour.groupSize(); ie++) solutionRel[ie].resize(local.num_local);
  }
  
  int maxkkrsz=(lsms.maxlmax+1)*(lsms.maxlmax+1);
  int maxkkrsz_ns=lsms.n_spin_cant*maxkkrsz;
  // for non spin canted, but spin polarized we store tau00_l for local site i for spin up and down as
  // tau00_l(*,i) and tau00_l(*,i + local.num_local)
  // i.e. tau00_l hase size (maxkkrsz_ns*maxkkrsz_ns,2*local.num_local)
  // n.b. n_spin_pola/n_spin_cant == 1 if non polarized or spin canted; == 2 iff collinear
  Matrix<Complex> tau00_l(maxkkrsz_ns*maxkkrsz_ns,local.num_local*lsms.n_spin_pola/lsms.n_spin_cant); // This would be cleaner as a std::vector<Matrix<Complex>>
  Matrix<Complex> dos(4,local.num_local);
  // dos=0.0;
  Matrix<Complex> dosck(4,local.num_local);
  // dosck=0.0;
  Array3d<Complex> dipole(6,4,local.num_local);
  // dipole=0.0;
  Array3d<Complex> green(local.maxjws(),4,local.num_local); // would be better as std::vector<Matrix<Complex>> to avoid problems with different jws.
  // Array3d<Complex> green(local.atom[0].jws,4,local.num_local); // would be better as std::vector<Matrix<Complex>> to avoid problems with different jws.
                                                               // or declare as lsms.global.iprpts instead!
  // green=0.0;
  // orbital dos and densities
  Matrix<Complex> dos_orb(3,local.num_local);
  Matrix<Complex> dosck_orb(3,local.num_local);
  Array3d<Complex> dens_orb(local.maxjws(),3,local.num_local);
  // Array3d<Complex> dens_orb(local.atom[0].jws,3,local.num_local);

// setup Device constant on GPU
  int maxNumLIZ=0;
#ifdef BUILDKKRMATRIX_GPU
  #pragma omp parallel for default(none) shared(lsms,local,deviceConstants)
  for(int i=0; i<local.num_local; i++)
  {
    setupForBuildKKRMatrix_gpu(lsms,local.atom[i],deviceConstants[i]);
    // setupForBuildKKRMatrix_gpu_opaque(lsms,local.atom[i],deviceConstants[i]);
  }
#endif
  for(int i=0; i<local.num_local; i++)
  {
    if(local.atom[i].numLIZ>maxNumLIZ) maxNumLIZ=local.atom[i].numLIZ;
  }
#if defined(ACCELERATOR_CULA) || defined(ACCELERATOR_LIBSCI) || defined(ACCELERATOR_CUDA_C)
  // initDStore(deviceStorage,maxkkrsz,lsms.n_spin_cant,maxNumLIZ,lsms.global.GPUThreads);
  deviceStorage->allocate(maxkkrsz,lsms.n_spin_cant,maxNumLIZ,lsms.global.GPUThreads);
#endif

// inside an omp for to ensure first touch
#pragma omp parallel for default(none) shared(local,dos,dosck,dipole,green)
  for(int i=0; i<local.num_local; i++)
  {
    for(int j=0; j<4; j++)
    {
      dos(j,i)=dosck(j,i)=0.0;
      for(int k=0; k<6; k++) dipole(k,j,i)=0.0;
      for(int k=0; k<local.atom[i].jws; k++) green(k,j,i)=0.0;
    }
    local.atom[i].resetLocalDensities();
  }

  timeEnergyContourIntegration_1=MPI_Wtime()-timeEnergyContourIntegration_1;

  double timeEnergyContourIntegration_2=MPI_Wtime();
  double timeCalculateAllTauMatrices=0.0;

// energy groups:
  int eGroupRemainder=nume%lsms.energyContour.groupSize();
  int numEGroups=nume/lsms.energyContour.groupSize()+std::min(1,eGroupRemainder);
  std::vector<int> eGroupIdx(numEGroups+1);
  for(int ig=0; ig<numEGroups; ig++) eGroupIdx[ig]=ig*lsms.energyContour.groupSize();
  eGroupIdx[numEGroups]=nume;

  for(int ig=0; ig<numEGroups; ig++)
  {
// solve single site problem
#ifdef USE_NVTX
    nvtxEventAttributes_t eventAttrib = {0};
    eventAttrib.version = NVTX_VERSION;
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    eventAttrib.colorType = NVTX_COLOR_ARGB;
    eventAttrib.color = 0x00ff00ff;
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
    eventAttrib.message.ascii = "singleScatter";
    nvtxRangePushEx(&eventAttrib);
#endif

    double timeSingleScatterers=MPI_Wtime();
    local.tmatStore=0.0;
    expectTmatCommunication(comm,local);

    if(lsms.global.iprint>=1) printf("calculate single scatterer solutions.\n");

    if(lsms.relativity!=full)
    {
#pragma omp parallel for default(none) shared(local,lsms,eGroupIdx,ig,egrd,solutionNonRel,vr_con)
      for(int ie=eGroupIdx[ig]; ie<eGroupIdx[ig+1]; ie++)
      {
        int iie=ie-eGroupIdx[ig];
        Complex energy=egrd[ie];
        Complex pnrel=std::sqrt(energy);
        
        solveSingleScatterers(lsms,local,vr_con,energy,solutionNonRel[iie],iie);
      }
    } else {
#pragma omp parallel for default(none) shared(local,lsms,eGroupIdx,ig,egrd,solutionRel,vr_con)
      for(int ie=eGroupIdx[ig]; ie<eGroupIdx[ig+1]; ie++)
      {
        int iie=ie-eGroupIdx[ig];
        Complex energy=egrd[ie];
        
        solveSingleScatterers(lsms,local,vr_con,energy,solutionRel[iie],iie);
      }
    }

    if(lsms.global.iprint>=2) printf("About to send t matrices\n");
    sendTmats(comm,local);
    if(lsms.global.iprint>=2) printf("About to finalize t matrices communication\n");
    finalizeTmatCommunication(comm);
    if(lsms.global.iprint>=2) printf("Recieved all t matricies\n");
    timeSingleScatterers=MPI_Wtime()-timeSingleScatterers;
#ifdef USE_NVTX
    nvtxRangePop();
#endif
    if(lsms.global.iprint>=0) printf("timeSingleScatteres = %lf sec\n",timeSingleScatterers);

#ifdef BUILDKKRMATRIX_GPU
  copyTmatStoreToDevice(local);
#endif

  for(int ie=eGroupIdx[ig]; ie<eGroupIdx[ig+1]; ie++)
  {
    int iie=ie-eGroupIdx[ig];
    Complex energy=egrd[ie];
    Complex pnrel=std::sqrt(energy);
    if(lsms.global.iprint>=1) printf("Energy #%d (%lf,%lf)\n",ie,real(energy),imag(energy));

    double timeCATM=MPI_Wtime();
    calculateAllTauMatrices(comm, lsms, local, vr_con, energy, iie, tau00_l);

    timeCalculateAllTauMatrices+=MPI_Wtime()-timeCATM;
    // if(!lsms.global.checkIstop("buildKKRMatrix"))
    {
#ifdef USE_NVTX
    nvtxEventAttributes_t eventAttrib = {0};
    eventAttrib.version = NVTX_VERSION;
    eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
    eventAttrib.colorType = NVTX_COLOR_ARGB;
    eventAttrib.color = 0x00ffff00;
    eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
    eventAttrib.message.ascii = "calculateDensities";
    nvtxRangePushEx(&eventAttrib);
#endif
    double timeCalcDensities=MPI_Wtime();
    if(lsms.relativity != full)
    {
// openMP here
#pragma omp parallel for default(none) \
        shared(local,lsms,dos,dosck,green,dipole,solutionNonRel,gauntCoeficients,dele1,tau00_l) \
        firstprivate(ie,iie,pnrel,energy,nume)
      for(int i=0; i<local.num_local; i++)
      {
        //Real r_sph=local.atom[i].r_mesh[local.atom[i].jws];
        //if (lsms.mtasa==0) r_sph=local.atom[i].r_mesh[local.atom[i].jmt];
        Real r_sph=local.atom[i].rInscribed;
        if(lsms.mtasa>0) r_sph=local.atom[i].rws;
        Real rins=local.atom[i].rmt;
	int jmt = local.atom[i].jmt;
	if(lsms.mtasa==1) jmt = local.atom[i].jws;
//        int nprpts=solutionNonRel[iie][i].zlr.l_dim1();
        int nprpts=local.atom[i].r_mesh.size();
//        int nplmax=solutionNonRel[iie][i].zlr.l_dim2()-1;
        int nplmax=local.atom[i].lmax;
        green_function_(&lsms.mtasa,&lsms.n_spin_pola,&lsms.n_spin_cant,
                        &local.atom[i].lmax, &local.atom[i].kkrsz,
                        &local.atom[i].wx[0],&local.atom[i].wy[0],&local.atom[i].wz[0],
                        &rins,&r_sph,&local.atom[i].r_mesh[0],&jmt,&local.atom[i].jws,
                        &pnrel,&tau00_l(0,i),&solutionNonRel[iie][i].matom(0,0),
                        &solutionNonRel[iie][i].zlr(0,0,0),&solutionNonRel[iie][i].jlr(0,0,0),
                        &nprpts,&nplmax,
                        &lsms.ngaussr, &gauntCoeficients.cgnt(0,0,0), &gauntCoeficients.lmax,
                        &dos(0,i),&dosck(0,i),&green(0,0,i),&dipole(0,0,i),
                        &local.atom[i].voronoi.ncrit,&local.atom[i].voronoi.grwylm(0,0),
                        &local.atom[i].voronoi.gwwylm(0,0),&local.atom[i].voronoi.wylm(0,0,0),
                        &lsms.global.iprint,lsms.global.istop,32);
	if((lsms.n_spin_pola == 2) && (lsms.n_spin_cant == 1)) // spin polarized, collinear case
	{
	  green_function_(&lsms.mtasa,&lsms.n_spin_pola,&lsms.n_spin_cant,
                        &local.atom[i].lmax, &local.atom[i].kkrsz,
                        &local.atom[i].wx[0],&local.atom[i].wy[0],&local.atom[i].wz[0],
                        &rins,&r_sph,&local.atom[i].r_mesh[0],&jmt,&local.atom[i].jws,
                        &pnrel,&tau00_l(0,i+local.num_local),&solutionNonRel[iie][i].matom(0,1),
                        &solutionNonRel[iie][i].zlr(0,0,1),&solutionNonRel[iie][i].jlr(0,0,1),
                        &nprpts,&nplmax,
                        &lsms.ngaussr, &gauntCoeficients.cgnt(0,0,0), &gauntCoeficients.lmax,
                        &dos(1,i),&dosck(1,i),&green(0,1,i),&dipole(0,0,i),
                        &local.atom[i].voronoi.ncrit,&local.atom[i].voronoi.grwylm(0,0),
                        &local.atom[i].voronoi.gwwylm(0,0),&local.atom[i].voronoi.wylm(0,0,0),
                        &lsms.global.iprint,lsms.global.istop,32);
	}

        if(local.atom[i].forceZeroMoment &&(lsms.n_spin_pola>1))
        {
          if(lsms.n_spin_cant>1) // spin canted case
          {
            for(int ir=0; ir<green.l_dim1(); ir++)
            {
              // green(ir,0,i) is the charge density part and green(ir,1:3,i) are the magnetic moment part
              green(ir,1,i) = green(ir,2,i) = green(ir,3,0) = 0.0;
            }
	    dos(1,i) = dos(2,i) = dos(3,i) = dosck(1,i) = dosck(2,i) = dos(3,i) = 0.0;
          } else { // spin polarized collinear case
            for(int ir=0; ir<green.l_dim1(); ir++)
            {
              // green(ir,0,i) is the spin up charge density and green(ir,1,i) is spin down part
              green(ir,0,i) = 0.5*(green(ir,0,i) + green(ir,1,i));
              green(ir,1,i) = green(ir,0,i);
            }
          }
        }

        Complex tr_pxtau[3];
        calculateDensities(lsms, i, 0, ie, nume, energy, dele1[ie],
                           dos,dosck,green,
                           dipole,
                           local.atom[i]);
	if((lsms.n_spin_pola == 2) && (lsms.n_spin_cant == 1)) // spin polarized, collinear case
	{
	  calculateDensities(lsms, i, 1, ie, nume, energy, dele1[ie],
			     dos,dosck,green,
			     dipole,
			     local.atom[i]);
	}

      }
    } else { // fully relativistic
      for(int i=0; i<local.num_local; i++)
      {
        //Real r_sph=local.atom[i].r_mesh[local.atom[i].jws];
        //if (lsms.mtasa==0) r_sph=local.atom[i].r_mesh[local.atom[i].jmt];
        Real r_sph=local.atom[i].rInscribed;
        if(lsms.mtasa>0) r_sph=local.atom[i].rws;
        Real rins=local.atom[i].rmt;
//        int nprpts=solutionNonRel[iie][i].zlr.l_dim1();
        int nprpts=local.atom[i].r_mesh.size();
//        int nplmax=solutionNonRel[iie][i].zlr.l_dim2()-1;
        int nplmax=local.atom[i].lmax;
        // printf("Relativistic version not implemented yet\n");
        // exit(1);
        
        green_function_rel_(&lsms.mtasa,
                            &local.atom[i].lmax, &local.atom[i].kkrsz,
                            &local.atom[i].wx[0],&local.atom[i].wy[0],&local.atom[i].wz[0],
                            &rins,&local.atom[i].r_mesh[0],
                            &local.atom[i].jws, // originally jmt
                            &local.atom[i].jws,&local.atom[i].h,
                            &pnrel,&tau00_l(0,i),&solutionRel[iie][i].matom(0,0),
                            &solutionRel[iie][i].gz(0,0,0),&solutionRel[iie][i].fz(0,0,0),
                            &solutionRel[iie][i].gj(0,0,0),&solutionRel[iie][i].fj(0,0,0),
                            &solutionRel[iie][i].nuz[0],&solutionRel[iie][i].indz(0,0),
                            &nprpts,
                            &lsms.ngaussr, &gauntCoeficients.cgnt(0,0,0), &gauntCoeficients.lmax,
                            &dos(0,i),&dosck(0,i),&green(0,0,i),&dipole(0,0,i),
                            &dos_orb(0,i),&dosck_orb(0,i),&dens_orb(0,0,i),
                            &lsms.global.iprint,lsms.global.istop,32);

        // rotateToGlobal(local.atom[i], dos, dosck, dos_orb, dosck_orb, green, dens_orb, i);
        
        // this is the non-rel version now
        calculateDensities(lsms, i, 0, ie, nume, energy, dele1[ie],
                           dos,dosck,green,
                           dipole,
                           local.atom[i]);
        
      }
    }
#ifdef USE_NVTX
    nvtxRangePop();
#endif
    timeCalcDensities=MPI_Wtime()-timeCalcDensities;
    if(lsms.global.iprint>=1) printf("timeCalculateDensities = %lf sec\n",timeCalcDensities);
  }
  }
  }
  timeEnergyContourIntegration_2=MPI_Wtime()-timeEnergyContourIntegration_2;
  if(lsms.global.iprint>=0)
  {
    printf("time in energyContourIntegration = %lf sec\n",timeEnergyContourIntegration_1+timeEnergyContourIntegration_2);
    printf("  before energy loop             = %lf sec\n",timeEnergyContourIntegration_1);
    printf("  in energy loop                 = %lf sec\n",timeEnergyContourIntegration_2);
    printf("    in calculateAllTauMatrices   = %lf sec\n",timeCalculateAllTauMatrices);
  }
}
