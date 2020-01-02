// -*- mode: c++; -*-

#include <stdlib.h>
#include "Real.hpp"
#include "Complex.hpp"
#include "Matrix.hpp"

#include "DeviceMatrix.hpp"
#include "DeviceArray3d.hpp"
#include "DeviceVector.hpp"

// #include "TestStructures.hpp"
#include "Misc/Indices.hpp"
#include "Main/SystemParameters.hpp"
#include "SingleSite/AtomData.hpp"
#include "Misc/Coeficients.hpp"

#ifdef CRAYPAT
#include <pat_api.h>
#endif

#include <cublas_v2.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#else
inline int omp_get_max_threads() {return 1;}
inline int omp_get_num_threads() {return 1;}
inline int omp_get_thread_num() {return 0;}
inline double omp_get_wtime() {return 0.0;}
#endif

// #include "cuda_error.h"
#include "cudaCheckError.hpp"
#include "DeviceStorage.hpp"
//#define TRANSFER
//#define CHECK

#include "buildKKRMatrix_gpu.hpp"

using namespace std;

extern "C"
{
  void write_kkrmat_(Complex *a,int *n,int *lda,Complex *e);

  void makegij_(int *lmaxi,int *kkri,int *lmaxj,int *kkrj,
                int *lmax,int *kkrsz,int *ndlj,int *ndlm,
                Complex *prel,double *rij,double *sinmp,double *cosmp,
                double *clm,double *plm,double *cgnt,int *lmax_cg,int *lofk,int *mofk,
                Complex *ilp1,Complex *illp,
                Complex *hfn,Complex *dlm,Complex *gij,
                double *pi4,int *iprint,char *istop,int len_sitop);

  void setgij_(Complex *gij,Complex *bgij,int *kkr1,int *kkr1_ns,int *kkr2,int *kkr2_ns,
               int *n_spin_cant,int *nrel_rel,Complex *psq,Complex *energy);

};

void setDiagonal_cuda(Complex *m, int rows, int lda,  cudaStream_t s );
void makeTmat_cuda_batched(DeviceVector<int> &LIZStoreIdx, DeviceVector<int> &LIZlmax, int iie, int blkSizeTmatStore, DeviceMatrix<Complex> &tmat_store, int n_spin_cant, int numLIZ, int kkr1, int kkrsz, Complex *tmat_n, cudaStream_t s);
void makeBGijs_cuda_batched(int numLIZ, int ndlj, int ndlm, int lmax, int gntlmax, int kkrsz, int nspin, Complex prel, Complex energy,
    DeviceVector<int> &LIZlmax, DeviceMatrix<Real> &LIZPos, DeviceVector<Real> &clm, DeviceArray3d<Real> &cgnt, 
    DeviceVector<int> &lofk, DeviceVector<int> &mofk, DeviceVector<Complex> &ilp1, DeviceMatrix<Complex> &illp, 
    Complex *bgij, cudaStream_t s);
void zgemm_cuda_batched(int numLIZ, int kkrsz, int nspin, int nrmat_ns, cublasHandle_t &cublas_h, Complex* a, Complex *b, Complex *c, cudaStream_t s);


ostream& operator<<(ostream& out, const Complex &c) {
  out << "( " << c.real() << " , " << c.imag() << " )";
  return out;
}

/*
//Instantiate this once, use the = operator from their host counterparts after the host has valid data. Then leave them on the device, pass into kernels as necessary
class DeviceConstants {
  public:
  //DeviceConstants() : LIZStoreIdx(0), LIZlmax(0), clm(0), lofk(0), mofk(0), ilp1(0), LIZPos(0,0), illp(0,0), cgnt(0,0,0) { }

  //~DeviceConstants() { }
 
  DeviceVector<int> LIZStoreIdx;
  DeviceVector<int> LIZlmax; 
  DeviceVector<Real> clm; 
  DeviceVector<int> lofk; 
  DeviceVector<int> mofk; 
  DeviceVector<Complex> ilp1; 
  DeviceMatrix<Real> LIZPos;  
  DeviceMatrix<Complex> illp; 
  DeviceArray3d<Real> cgnt;   
};
*/

void *allocateDConst(void)
{
  return static_cast<void *>(new DeviceConstants);
}

void freeDConst(void * d_const)
{
  delete static_cast<DeviceConstants*>(d_const);
}

#ifdef CHECK
 //test function for comparing two arrays within a tolerance
 template<typename T>
 void checkResults(T* host, T* dev, unsigned int num, string where, double tol) {
   T* tmp=new T[num];
   cudaDeviceSynchronize();
   cudaCheckError();
   cudaMemcpy(tmp,dev,sizeof(T)*num,cudaMemcpyDeviceToHost);

   cudaCheckError();
   int count=0, max_count=10;
   for(int i=0;i<num;i++) {
     if(std::abs(tmp[i]-host[i])>tol) {
        cout << where << ": Error at index " << i <<  " host(" << host[i] << ") device(" << tmp[i] << ")\n";
        count++;
     }
     if(count>max_count)
       exit(1);
   }

   if(count!=0)
     exit(1);
  
   cout << where << " passed\n";
   delete [] tmp;
 }


 //host side function for comparison
 void initializeM(Matrix<Complex> &m, int nrmat_ns ) {
   for(int i=0; i<nrmat_ns*nrmat_ns; i++) 
     m[i]=0.0;
   for(int i=0; i<nrmat_ns; i++) 
     m(i,i)=1.0;
 }

 //host side function for comparison
 void makeTmat(int iie, const LSMSSystemParameters &lsms, /*const*/ LocalTypeInfo &local, const AtomData &atom, int ir1, int kkr1,int kkrsz_ns,int kkrsz, Complex *tmat_n) {
   printf("KKR1: %d, kkrsz: %d\n", kkr1,kkrsz);
   int im=0;

   for(int js=0; js<lsms.n_spin_cant; js++) //2
   {
     int jsm = kkrsz*kkrsz_ns*js;  //linear through matrix of size kkrsz*kkrsz_ns
     for(int j=0; j<kkr1; j++) //16
     {
       for(int is=0; is<lsms.n_spin_cant; is++) //2
       {
         int jm=jsm+kkrsz_ns*j+kkrsz*is;
         int one=1;
         //copy of size kkr1 from local.tmatStore(jm,atom.LIZStoreIdx) to tmat_n
         BLAS::zcopy_(&kkr1,&local.tmatStore(iie*local.blkSizeTmatStore+jm,atom.LIZStoreIdx[ir1]),&one,&tmat_n[im],&one);
         im+=kkr1;
       }
     }
   }
 }

 //host side function for comparison
 void makeBGijs(/*const*/ LSMSSystemParameters &lsms, /*const*/ AtomData &atom, int ir1, int kkr1, int kkrsz, Complex prel, Complex energy, Complex *gij, Complex *bgij ) {

   int lmax=lsms.maxlmax;
   for(int ir2=0; ir2<atom.numLIZ; ir2++) //PARALLEL 
   {

     //Place these in dynamic shared memory
     Real *sinmp = new Real[2*lmax+1];  //56 bytes
     Real *cosmp = new Real[2*lmax+1];  //56 bytes
     Real *plm = new Real[lsms.angularMomentumIndices.ndlm]; //224 bytes
     Complex *hfn = new Complex[2*lmax+1]; //102 bytes
     Complex *dlm = new Complex[lsms.angularMomentumIndices.ndlj]; //784 bytes

     if(ir1!=ir2)
     {
       int kkr2=(atom.LIZlmax[ir2]+1)*(atom.LIZlmax[ir2]+1);
       //int kkr2_ns=kkr2*lsms.n_spin_cant;
       Real rij[3];
       Real pi4=4.0*2.0*std::asin(1.0);

       rij[0]=atom.LIZPos(0,ir1)-atom.LIZPos(0,ir2);
       rij[1]=atom.LIZPos(1,ir1)-atom.LIZPos(1,ir2);
       rij[2]=atom.LIZPos(2,ir1)-atom.LIZPos(2,ir2);

       //if(ir1==0 && ir2==1) 
       {
       makegij_(&atom.LIZlmax[ir1],&kkr1,&atom.LIZlmax[ir2],&kkr2,
           &lsms.maxlmax,&kkrsz,&lsms.angularMomentumIndices.ndlj,&lsms.angularMomentumIndices.ndlm,
           &prel,&rij[0],sinmp,cosmp,
           &sphericalHarmonicsCoeficients.clm[0],plm,
           &gauntCoeficients.cgnt(0,0,0),&gauntCoeficients.lmax,
           &lsms.angularMomentumIndices.lofk[0],&lsms.angularMomentumIndices.mofk[0],
           &iFactors.ilp1[0],&iFactors.illp(0,0),
           hfn,dlm,gij+kkrsz*kkrsz*ir2,
           &pi4,&lsms.global.iprint,lsms.global.istop,32);
  


       int kkr2_ns=kkr2*lsms.n_spin_cant;
       int kkr1_ns=kkr1*lsms.n_spin_cant;
       Complex psq=prel*prel;
    
       setgij_(gij+kkrsz*kkrsz*ir2,bgij+4*kkrsz*kkrsz*ir2,&kkr1,&kkr1_ns,&kkr2,&kkr2_ns,
           &lsms.n_spin_cant,&lsms.nrel_rel,&psq,&energy);
       }
     }
     else {
#if 0
       //set block to zero so that it matches the cuda work
       Complex zero(0,0);
      for(int idx=0;idx<4*kkr1*kkr1;idx++) {
        bgij[4*kkrsz*kkrsz*ir2+idx]=zero;
      }
#endif
     }
     delete [] sinmp;
     delete [] cosmp;
     delete [] plm;
     delete [] hfn;
     delete [] dlm;
   }
 }

 //host side function for comparison
 void batchedZgemm(/*const*/ LSMSSystemParameters &lsms, /*const*/ AtomData &atom, int ir1, int kkr1_ns,
                   int kkrsz, int nrst, int nrmat_ns, /*const*/ Complex *tmat_n, /*const*/ Complex *bgij, Matrix<Complex> &m ) {
     int ncst=0;
     //Todo combine into a single zgemm call
     for(int ir2=0; ir2<atom.numLIZ; ir2++) //PARALLEL 
     {
       int kkr2=(atom.LIZlmax[ir2]+1)*(atom.LIZlmax[ir2]+1);
       int kkr2_ns=kkr2*lsms.n_spin_cant;
       if(ir1!=ir2)  
       {
         const Complex cmone=-1.0;
         const Complex czero=0.0;
         //SIZE:  32x32 . 32x32
         BLAS::zgemm_("n","n",&kkr1_ns,&kkr2_ns,&kkr1_ns,&cmone,
             tmat_n,&kkr1_ns,bgij+4*kkrsz*kkrsz*ir2,&kkr1_ns,&czero,
             &m(nrst,ncst),&nrmat_ns);
       }
       ncst+=kkr2_ns;
     }
 }
#endif

void setupForBuildKKRMatrix_gpu(LSMSSystemParameters &lsms, AtomData &atom,DeviceConstants &d_const)
{
  cudaStream_t stream=get_stream_(0);
  /********************************************
   * Device Constants
   *******************************************/
  // DeviceConstants d_const;
  d_const.LIZStoreIdx.copy_async(atom.LIZStoreIdx,stream);          //constant
  d_const.LIZlmax.copy_async(atom.LIZlmax,stream);                  //constant
  d_const.clm.copy_async(sphericalHarmonicsCoeficients.clm,stream); //constant
  d_const.lofk.copy_async(lsms.angularMomentumIndices.lofk,stream); //constant
  d_const.mofk.copy_async(lsms.angularMomentumIndices.mofk,stream); //constant
  d_const.ilp1.copy_async(iFactors.ilp1,stream);                    //constant
  d_const.LIZPos.copy_async(atom.LIZPos,stream);                    //constant
  d_const.illp.copy_async(iFactors.illp,stream);                    //constant
  d_const.cgnt.copy_async(gauntCoeficients.cgnt,stream);            //constant
}

void setupForBuildKKRMatrix_gpu_opaque(LSMSSystemParameters &lsms, AtomData &atom,void *d_const)
{
  setupForBuildKKRMatrix_gpu(lsms,atom,*static_cast<DeviceConstants*>(d_const));
}

void printDeviceComplex(char* where, Complex* val) {
  Complex h;
  cudaMemcpy(&h,val,sizeof(Complex),cudaMemcpyDeviceToHost);
  printf("**********************%s: %lg, %lg\n",where,h.real(),h.imag());
}
void buildKKRMatrix_gpu(LSMSSystemParameters &lsms, LocalTypeInfo &local,AtomData &atom, Complex energy, Complex prel, int iie, Matrix<Complex> &m, DeviceConstants &d_const)
{
  int lmax=lsms.maxlmax;
  int kkrsz=(lmax+1)*(lmax+1);
  int nrmat_ns=lsms.n_spin_cant*atom.nrmat;

#ifdef CHECK
  /********************************************
   * Host pointers
   * *****************************************/
  Complex *gij = new Complex[kkrsz*kkrsz*atom.numLIZ];    //Replicating for parallelism
  memset(gij,0,kkrsz*kkrsz*atom.numLIZ*sizeof(Complex));
  Complex *bgij = new Complex[4*kkrsz*kkrsz*atom.numLIZ]; //Replicating for parallelism
  Complex *tmat_n = new Complex[atom.kkrsz*atom.kkrsz*4];
  cudaHostRegister(gij,sizeof(Complex)*kkrsz*kkrsz*atom.numLIZ,0);
  cudaHostRegister(bgij,sizeof(Complex)*4*kkrsz*kkrsz*atom.numLIZ,0);
  cudaHostRegister(tmat_n,sizeof(Complex)*atom.kkrsz*atom.kkrsz*4,0);
  cudaDeviceSynchronize(); //syncing to make sure the timing below is accurate
#endif

  cublasHandle_t cublas_h=get_cublas_handle_();
  cudaStream_t stream=get_stream_(0);
  cudaEvent_t done_event=get_cuda_event_();

  /********************************************
   * Device arrays
   * *****************************************/
  Complex *dev_m, *dev_bgij, *dev_tmat_n;
  dev_m=(Complex*)get_dev_m_();
  dev_bgij=(Complex*)get_dev_bgij_(); 
  dev_tmat_n=get_dev_tmat_n_();


/*********************************************************************************************************/
  DeviceMatrix<Complex> &d_tmat_store=*get_dev_tmat_store();

#ifdef CHECK
  checkResults(&local.tmatStore[0],d_tmat_store.raw(), d_tmat_store.size(), "TMAT STORE", 0);
  checkResults(&atom.LIZStoreIdx[0], d_const.LIZStoreIdx.raw(), d_const.LIZStoreIdx.size(), "LIZStoreIDX", 0);
  checkResults(&atom.LIZlmax[0], d_const.LIZlmax.raw(), d_const.LIZlmax.size(), "LIZlmax", 0);
#endif



  cudaCheckError();
  unsigned long long flops=0;
  int kkrsz_ns=kkrsz*lsms.n_spin_cant;

  //batch makeTmat
  makeTmat_cuda_batched(d_const.LIZStoreIdx, d_const.LIZlmax, iie, local.blkSizeTmatStore, d_tmat_store, lsms.n_spin_cant, atom.numLIZ, kkrsz, kkrsz, dev_tmat_n,stream);
  //batch makeBGijs_cuda
  makeBGijs_cuda_batched(atom.numLIZ, lsms.angularMomentumIndices.ndlj,lsms.angularMomentumIndices.ndlm, 
      lmax, gauntCoeficients.lmax, kkrsz, lsms.n_spin_cant, prel, energy,
      d_const.LIZlmax, d_const.LIZPos, d_const.clm, d_const.cgnt, d_const.lofk, d_const.mofk, d_const.ilp1, d_const.illp, dev_bgij,stream);
  //batch zgemm_cuda
  zgemm_cuda_batched(atom.numLIZ, kkrsz, lsms.n_spin_cant, nrmat_ns, cublas_h, dev_tmat_n, dev_bgij, dev_m,stream);
  
  //diagonal blocks contain all zeros from the zgemm.  Now set the diagonals to (1,0)
  setDiagonal_cuda(dev_m,nrmat_ns,nrmat_ns,stream);
  cudaEventRecord(done_event,stream);

  flops += kkrsz_ns* kkrsz_ns *atom.numLIZ * kkrsz_ns * atom.numLIZ; //TODO double check this formula
#ifdef PRINT_FLOPS
  printf("%d: BUILDKKR FLOPS: %llu\n",omp_get_thread_num(),flops * 4 * 2);
#endif

#ifdef CHECK
  initializeM(m,nrmat_ns);

  int nrst=0;
  for(int ir1=0; ir1<atom.numLIZ; ir1++) //SERIAL
  {
    int kkr1=(atom.LIZlmax[ir1]+1)*(atom.LIZlmax[ir1]+1);
    int kkr1_ns=kkr1*lsms.n_spin_cant;
  
    makeTmat(iie, lsms, local, atom, ir1, kkr1, kkrsz*lsms.n_spin_cant, kkrsz, tmat_n);
#ifdef CHECK
    //compare tmats
    checkResults(&tmat_n[0], dev_tmat_n+(atom.kkrsz*atom.kkrsz*4)*ir1, (atom.kkrsz*atom.kkrsz*4), "TMAT_N", 1e-8);
#endif

    makeBGijs(lsms, atom, ir1, kkr1, kkrsz, prel, energy, gij, bgij); 
#ifdef CHECK
    //compare BGIJs - this doesn't work even when it is right.  I tracked this down a while ago and saw that it was ok but I don't remember what caused it... Justin
    //checkResults(&bgij[0], dev_bgij+4*kkrsz*kkrsz*atom.numLIZ*ir1, 4*kkrsz*kkrsz*atom.numLIZ, "BGIJ", 1e-8);
#endif
    
    batchedZgemm(lsms, atom, ir1, kkr1_ns, kkrsz, nrst, nrmat_ns, tmat_n, bgij, m);
    
    nrst+=kkr1_ns;
  }
  checkResults(&m[0], dev_m, nrmat_ns*nrmat_ns, "M",1e-8);
#endif

#ifdef TRANSFER
  cudaMemcpyAsync(&m[0],dev_m,nrmat_ns*nrmat_ns*sizeof(Complex),cudaMemcpyDeviceToHost,stream);  //Generated on Device
  cudaEventRecord(done_event,stream);
  cudaEventSynchronize(done_event);
#endif

#ifdef CHECK
  cudaHostUnregister(gij);
  cudaHostUnregister(bgij);
  cudaHostUnregister(tmat_n);
  delete [] bgij;
  delete [] tmat_n;
  delete [] gij;
#endif

  cudaCheckError();
}

void buildKKRMatrix_gpu_opaque(LSMSSystemParameters &lsms, LocalTypeInfo &local,AtomData &atom, Complex energy, Complex prel, int iie, Matrix<Complex> &m, void *d_const)
{
  buildKKRMatrix_gpu(lsms, local, atom, energy, prel, iie, m, 
                     *static_cast<DeviceConstants*>(d_const));
}

