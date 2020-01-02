#include "mixing.hpp"
#include "Communication/LSMSCommunication.hpp"

// the modified Broyden method follows D. D. Johnson, PRB 38, 12807
template <typename T>
class BroydenMixing {

  int vectorSize;       // size of vectors to be mixed
  int currentIteration; // currently used number of iterations
  int iterationReset;   // number of iterations after which the Broyden mixing is reset
  int maxBroydenLength; // maximum number of iterations that are used
  std::vector<T> vOld, F, dF;
  std::vector<Real> w, cm;
  std::vector<std::vector<T> > u,vt;
  std::vector<int> ipiv;
  Real alpha, w0;
  Matrix<Real> a,b;

public:
  void save(std::vector<T> &fOld, std::vector<T> &fNew, Real wtmp)
  {
    if(currentIteration<maxBroydenLength+1)
    {
      for(int i=0; i<vectorSize; i++)
      {
        u[currentIteration-1][i]=fNew[i];
        vt[currentIteration-1][i]=fOld[i];
      }
      w[currentIteration-1]=wtmp;
    } else {
      for(int j=0; j<maxBroydenLength-1; j++)
      {
        for(int i=0; i<vectorSize; i++)
        {
          u[j][i]=u[j+1][i];
          vt[j][i]=vt[j+1][i];
        }
        w[j]=w[j+1];
      }
      for(int i=0; i<vectorSize; i++)
      {
        u[maxBroydenLength-1][i]=fNew[i];
        vt[maxBroydenLength-1][i]=fOld[i];
      }
      w[maxBroydenLength-1]=wtmp;
    }
  }



// a <- a^-1
  void invert(Matrix<Real> &A,int nn)
  {
    int LDIM=A.l_dim();
    int LWORK = LDIM*LDIM;
    double *WORK = new double[LWORK];
    int INFO;

    LAPACK::dgetrf_(&nn,&nn,&A(0,0),&LDIM,&ipiv[0],&INFO);
    LAPACK::dgetri_(&nn,&A(0,0),&LDIM,&ipiv[0],WORK,&LWORK,&INFO);

    delete[] WORK;
  }

  void init(Real a_, int vs, int mbl=10, int itres=25)
  {
    alpha=a_; vectorSize=vs; currentIteration=0; maxBroydenLength=mbl; iterationReset=itres;
    vOld.resize(vectorSize); F.resize(vectorSize); dF.resize(vectorSize);
    w.resize(maxBroydenLength); cm.resize(maxBroydenLength);
    a.resize(mbl,mbl); b.resize(mbl,mbl); ipiv.resize(mbl+1);
    vt.resize(maxBroydenLength); u.resize(maxBroydenLength);
    for(int i=0; i<maxBroydenLength; i++)
    {
      vt[i].resize(vectorSize); u[i].resize(vectorSize);
    }
    w0=0.01;
  }

  void mix(LSMSCommunication &comm, std::vector<T> &fOld, std::vector<T> &fNew, Real rms)
  {
    if(currentIteration==0)
    {
      // first iteration: perform linear mixing, set up internal storage
      for(int i=0; i<vectorSize; i++)
      {
        F[i]=fNew[i]-fOld[i];
        vOld[i]=fOld[i];
        fNew[i]=fOld[i]+alpha*F[i];
      }
    } else {
      int nn=std::min(currentIteration,maxBroydenLength);
      Real dFnorm=0.0;
      Real Fnorm=0.0;
      for(int i=0; i<vectorSize; i++)
      {
        dF[i]=fNew[i] - fOld[i] -F[i];
        F[i]=fNew[i]-fOld[i];
        dFnorm+=dF[i]*dF[i];
        Fnorm+=F[i]*F[i];
      }
///* // Try without communication!
      globalSum(comm, dFnorm);
      globalSum(comm, Fnorm);
//*/
      dFnorm=std::sqrt(dFnorm);
      Fnorm=std::sqrt(Fnorm);
      Real fac2=1.0/dFnorm;
      Real fac1=alpha*fac2;
      for(int i=0; i<vectorSize; i++)
      {
        fNew[i]=fac1*dF[i]+fac2*(fOld[i]-vOld[i]);
        vOld[i]=fOld[i];
        fOld[i]=fac2*dF[i];
      }

      Real wtmp=0.0;
      if(rms<1.0e-9) wtmp=2.0*std::sqrt(0.01/rms);
      if(wtmp<1.0) wtmp=1.0;

      save(fOld,fNew,wtmp);

      // off diagonal part of a(i,j)
      for(int j=0; j<nn-1; j++)
        for(int i=j+1; i<nn; i++)
        {
          Real aij=0.0;
          for(int k=0; k<vectorSize; k++)
            aij+=(vt[j][k])*(vt[i][k]);
          a(i,j)=a(j,i)=aij;
        }
      // diagonal elements a(i,i) and cm(i)
      for(int i=0; i<nn; i++)
      {
        Real aij=0.0;
        Real cmj=0.0;
        for(int k=0; k<vectorSize; k++)
        {
          cmj+=vt[i][k]*F[k];
          aij+=vt[i][k]*vt[i][k];
        }
        a(i,i)=aij;
        cm[i]=cmj;
      }
      // sum over all sites
      globalSum(comm, &a(0,0),maxBroydenLength*maxBroydenLength);
      globalSum(comm, &cm[0],maxBroydenLength);

// now calculate the b-matrix
//  b = [w(0)*w(0)*delta(i,j) + w(i)*w(j)*a(i,j)]^-1
//
      for(int i=0; i<nn; i++)
      {
        for(int j=0; j<nn; j++)
        {
          b(i,j)=a(i,j)*w[j]*w[i];
        }
        b(i,i)=w0*w0 + a(i,i)*w[i]*w[i];
      }
      invert(b,nn);

// mix vectors
      for(int k=0; k<vectorSize; k++)
        fNew[k]=vOld[k]+alpha*F[k];

      for(int i=0; i<nn; i++)
      {
        Real gmi=0.0;
        for(int j=0; j<nn; j++)
          gmi+=cm[j]*b(j,i)*w[j];
        for(int k=0; k<vectorSize; k++)
          fNew[k]=fNew[k]-gmi*u[i][k]*w[i];
      }

    }
    currentIteration++;
    if(iterationReset > 0 && currentIteration>iterationReset)
      currentIteration=0;
  }
};

Mixing::~Mixing() {}

class NoMixing : public Mixing {
public:
  // void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a) {}
  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}
  // void updatePotential(LSMSSystemParameters &lsms, AtomData &a) {}
  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}
  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}
};



class FrozenPotential : public Mixing {
 
public:

  // void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a)
  // {
  //   a.rhotot = a.rhoNew;
  //   a.xvalws[0] = a.xvalwsNew[0];
  //   a.xvalws[1] = a.xvalwsNew[1];
  // }

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    for(int i = 0; i < as.size(); i++) {
      as[i].rhotot = as[i].rhoNew;
      as[i].xvalws[0] = as[i].xvalwsNew[0];
      as[i].xvalws[1] = as[i].xvalwsNew[1];
    }
  }

  // void updatePotential(LSMSSystemParameters &lsms, AtomData &a) {}

  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}

  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}

};



class SimpleChargeDensityMixing : public Mixing {

public:
  Real alpha;

  SimpleChargeDensityMixing(Real _alpha) : alpha(_alpha) {}

  // void updateChargeDensity(LSMSSystemParameters &lsms, AtomData &a)
  // { 
  //   simpleMixing( &a.rhotot(0,0), &a.rhoNew(0,0), a.rhotot.size(), alpha);
  //   simpleMixing( &a.xvalws[0], &a.xvalwsNew[0], 2, alpha);
  // }

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  { 
    for(int i = 0; i < as.size(); i++)
    { 
      simpleMixing( &as[i].rhotot(0,0), &as[i].rhoNew(0,0), as[i].rhotot.size(), alpha);
      simpleMixing( &as[i].xvalws[0], &as[i].xvalwsNew[0], 2, alpha);
    }
  }

  // void updatePotential(LSMSSystemParameters &lsms, AtomData &a)
  // { 
  //   a.vr = a.vrNew;
  // }

  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  { 
    for(int i = 0; i < as.size(); i++) {
      as[i].vr   = as[i].vrNew;
      as[i].vdif = as[i].vdifNew;
    }
  }

  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}

};



class SimplePotentialMixing : public Mixing {

public:
  Real alpha;
  
  SimplePotentialMixing(Real _alpha) : alpha(_alpha) {}
  
  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  { 
    for(int i = 0; i < as.size(); i++) {
      as[i].rhotot = as[i].rhoNew;
      as[i].xvalws[0] = as[i].xvalwsNew[0];
      as[i].xvalws[1] = as[i].xvalwsNew[1];
    }
  }   
  
  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  { 
    for(int i = 0; i < as.size(); i++)
    {
      simpleMixing( &as[i].vr(0,0), &as[i].vrNew(0,0), as[i].vr.size(), alpha);
      simpleMixing( &as[i].vdif, &as[i].vdifNew, 1, alpha);
    }
  }
  
  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}

};


class BroydenChargeDensityMixing : public Mixing {
  BroydenMixing<Real> mixer;
  std::vector<Real> fNew, fOld;
  int vSize;
  std::vector<size_t> vStarts;
  Real alpha;

public:
  BroydenChargeDensityMixing(Real _alpha) : alpha(_alpha) {}

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    Real rms = 0.0;

    // first: copy potentials into fNew vector before mixing
    for (int i=0; i<as.size(); i++)
    {
      for (int j=0; j<as[i].rhotot.n_row(); j++)
      {
        fNew[vStarts[i]+j]       = as[i].rhoNew(j,0);
        fNew[vStarts[i]+j+vSize] = as[i].rhoNew(j,1);
        fOld[vStarts[i]+j]       = as[i].rhotot(j,0);
        fOld[vStarts[i]+j+vSize] = as[i].rhotot(j,1);

      }
      fNew[2*vSize]   = as[i].xvalwsNew[0];
      fNew[2*vSize+1] = as[i].xvalwsNew[1];
      fOld[2*vSize]   = as[i].xvalws[0];
      fOld[2*vSize+1] = as[i].xvalws[1];

      rms += as[i].qrms[0] + as[i].qrms[1];
    }
    rms = rms / (2.0 * as.size());
    globalSum(comm,rms);
    rms /= comm.size;

    // Broyden mixing
    mixer.mix(comm, fOld, fNew, rms);

    // copy mixed results back
    for(int i=0; i<as.size(); i++)
    {
      for(int j=0; j<as[i].rhotot.n_row(); j++)
      {
        as[i].rhotot(j,0) = fNew[vStarts[i]+j];
        as[i].rhotot(j,1) = fNew[vStarts[i]+j+vSize];
      }
      as[i].xvalws[0] = fNew[2*vSize];
      as[i].xvalws[1] = fNew[2*vSize+1];
    }
  }
 
  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    for(int i = 0; i < as.size(); i++)
    {
      as[i].vr   = as[i].vrNew;
      as[i].vdif = as[i].vdifNew;
    }  
  }
 
  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    vSize = 0;
    vStarts.resize(as.size());

    for(int i=0; i<as.size(); i++)
    { 
      vStarts[i] = vSize;
      vSize += as[i].rhotot.n_row();
    }
    mixer.init(alpha, 2*vSize + 2);
    fNew.resize(2*vSize + 2);
    fOld.resize(2*vSize + 2);
  }

};

class BroydenPotentialMixing : public Mixing {
  BroydenMixing<Real> mixer;
  std::vector<Real> fNew, fOld;
  size_t vSize;
  std::vector<size_t> vStarts;
  Real alpha;

public:
  BroydenPotentialMixing(Real _alpha) : alpha(_alpha) {}

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  { 
    for (int i=0; i<as.size(); i++) {
      as[i].rhotot    = as[i].rhoNew;
      as[i].xvalws[0] = as[i].xvalwsNew[0];
      as[i].xvalws[1] = as[i].xvalwsNew[1];
    }
  }

  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    Real rms = 0.0;

    // first: copy potentials into fNew vector before mixing
    for (int i=0; i<as.size(); i++)
    {
      for (int j=0; j<as[i].vr.n_row(); j++)
      {
        fNew[vStarts[i]+j]       = as[i].vrNew(j,0);
        fNew[vStarts[i]+j+vSize] = as[i].vrNew(j,1);
        fOld[vStarts[i]+j]       = as[i].vr(j,0);
        fOld[vStarts[i]+j+vSize] = as[i].vr(j,1);
      }
      fNew[2*vSize] = as[i].vdifNew;
      fOld[2*vSize] = as[i].vdif;
      
      rms += as[i].vrms[0] + as[i].vrms[1];
    }
    rms = rms / (2.0 * as.size());
    globalSum(comm,rms);
    rms /= comm.size;

    // Broyden mixing
    mixer.mix(comm, fOld, fNew, rms);

    // copy mixed results back
    for (int i=0; i<as.size(); i++)
    {
      for (int j=0; j<as[i].vr.n_row(); j++)
      {
        as[i].vr(j,0) = fNew[vStarts[i]+j];
        as[i].vr(j,1) = fNew[vStarts[i]+j+vSize];
      }
      as[i].vdif = fNew[2*vSize];
    }
  }
  
  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    vSize = 0;
    vStarts.resize(as.size());

    for (int i=0; i<as.size(); i++)
    {
      vStarts[i] = vSize;
      vSize += as[i].vr.n_row();
    }
    mixer.init(alpha, 2*vSize + 1);
    fNew.resize(2*vSize + 1);
    fOld.resize(2*vSize + 1);
  }

};


class EfMixing : public Mixing {
  Real efOld, alpha;

public:

  EfMixing(Real _alpha) : alpha(_alpha) {}

  void updateChargeDensity(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    lsms.chempot = alpha * lsms.chempot + (1.0-alpha) * efOld;
  }

  void updatePotential(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as) {}

  void prepare(LSMSCommunication &comm, LSMSSystemParameters &lsms, std::vector<AtomData> &as)
  {
    efOld = lsms.chempot;
  }

};



void setupMixing(MixingParameters &mix, Mixing* &mixing, int iprint)
{

  if(iprint >= 0)
    printf("\n");
  mixing = NULL;

  // frozen potential by default
  if (!mix.quantity[MixingParameters::no_mixing] &&
      !mix.quantity[MixingParameters::charge] && 
      !mix.quantity[MixingParameters::potential] && 
      !mix.quantity[MixingParameters::moment_magnitude] && 
      !mix.quantity[MixingParameters::moment_direction])
  {
    mixing = new FrozenPotential;
    if(iprint >= 0)
      printf("Mixing method     : frozen potential (default)\n");
  }
  // no mixing
  else if (mix.quantity[MixingParameters::no_mixing])
  {
    mixing = new NoMixing;
    if(iprint >= 0)
      printf("Mixing method     : no mixing\n");
  }
  // charge mixing
  else if (!mix.quantity[MixingParameters::no_mixing] &&
            mix.quantity[MixingParameters::charge] && 
           !mix.quantity[MixingParameters::potential] && 
           !mix.quantity[MixingParameters::moment_magnitude] && 
           !mix.quantity[MixingParameters::moment_direction])
  {
    switch (mix.algorithm[MixingParameters::charge]) {
      case 1 :
        mixing = new SimpleChargeDensityMixing(mix.mixingParameter[MixingParameters::charge]);
        if(iprint >= 0)
          printf("Mixing method     : simple\n");
        break;
      case 2 :
        mixing = new BroydenChargeDensityMixing(mix.mixingParameter[MixingParameters::charge]);
        if(iprint >= 0)
          printf("Mixing method     : broyden\n");
        break;
      default :
        mixing = new NoMixing;
        if(iprint >= 0)
          printf("Mixing method     : no mixing\n");
    }
    if(iprint >= 0){
      printf("Mixing quantity   : charge\n");
      printf("Mixing parameters : %4.2f\n", mix.mixingParameter[MixingParameters::charge]);
    }
  }
  // potential mixing
  else if (!mix.quantity[MixingParameters::no_mixing] &&
           !mix.quantity[MixingParameters::charge] && 
            mix.quantity[MixingParameters::potential] && 
           !mix.quantity[MixingParameters::moment_magnitude] && 
           !mix.quantity[MixingParameters::moment_direction])
  {
    switch (mix.algorithm[MixingParameters::potential]) {
      case 1 :
        if (mix.mixingParameter[MixingParameters::potential] == 0.0) {
          mixing = new FrozenPotential;
          if(iprint >= 0)
            printf("Mixing method     : frozen potential\n");
        }
        else {
          mixing = new SimplePotentialMixing(mix.mixingParameter[MixingParameters::potential]);
          if(iprint >= 0)
            printf("Mixing method     : simple\n");
        }
        break;
      case 2 :
        mixing = new BroydenPotentialMixing(mix.mixingParameter[MixingParameters::potential]);
        if(iprint >= 0)
          printf("Mixing method     : broyden\n");
        break;
      default :
        mixing = new FrozenPotential;
        if(iprint >= 0)
          printf("Mixing method     : frozen potential\n");
    }
    if(iprint >= 0){
      printf("Mixing quantity   : potential\n");
      printf("Mixing parameters : %4.2f\n", mix.mixingParameter[MixingParameters::potential]);
    }
  }
  else
  {
    if(iprint >= 0) {
      printf("Type of mixing is not supported.\n");
      for (int i = 0; i < mix.numQuantities; i++) {
        printf("quantity = %5d, algorithm = %5d, mixing parameter = %6.3f\n", 
               mix.quantity[i], mix.algorithm[i], mix.mixingParameter[i]);
      }
      exit(1);
    }
  }

}

