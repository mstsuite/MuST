#ifndef LSMS_MULTIPLESCATTERING_H
#define LSMS_MULTIPLESCATTERING_H

#include "Complex.hpp"
#include "Matrix.hpp"
#include "Main/SystemParameters.hpp"
#include "Communication/LSMSCommunication.hpp"

void calculateAllTauMatrices(LSMSCommunication &comm,LSMSSystemParameters &lsms, LocalTypeInfo &local,
                             std::vector<Matrix<Real> > &vr, Complex energy,
                             int iie,
                             // std::vector<NonRelativisticSingleScattererSolution> &solution,
                             Matrix<Complex> &tau00_l);

extern "C"
{
  void makegij_(int *lmaxi,int *kkri,int *lmaxj,int *kkrj,
                int *lmax,int *kkrsz,int *ndlj,int *ndlm,
                Complex *prel,double *rij,double *sinmp,double *cosmp,
                double *clm,double *plm,double *cgnt,int *lmax_cg,int *lofk,int *mofk,
                Complex *ilp1,Complex *illp,
                Complex *hfn,Complex *dlm,Complex *gij,
                double *pi4,int *iprint,char *istop,int len_sitop);

  void setgij_(Complex *gij,Complex *bgij,int *kkr1,int *kkr1_ns,int *kkr2,int *kkr2_ns,
               int *n_spin_cant,int *nrel_rel,Complex *psq,Complex *energy);

  void block_inv_(Complex *a, Complex *vecs, int *lda, int *na, int *mp, int *ipvt, int *blk_sz, int *nblk, Complex *delta,
                  int *iwork, double *rwork, Complex *work1, int *alg,
                  int *idcol, int *iprint);

  void tau_inv_postproc_nrel_(int *kkrs_ns,int *n_spin_cant,
                              Complex *wbig,Complex *delta, Complex *tmat,int *ipvt,Complex *tau00,
                              Complex *ubr, Complex *ubrd,
                              Complex *tau00_tmp);

 void tau_inv_postproc_rel_(int *kkrs_ns,
                            Complex *wbig,Complex *delta, Complex *tmat,int *ipvt,Complex *tau00,
                            Complex *dmatp, Complex *dmat,
                            Complex *tau00_tmp);
}
#endif
