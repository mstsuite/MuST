// -*- mode: c++ -*-
// Graph1d class for 1d Wang-Landau Histogramm and DOS

#ifndef WL_GRAPH1D_H
#define WL_GRAPH1D_H

#include <limits>
#include <cmath>
#include <stdlib.h>
#include<vector>

template <class ValueType, class KeyType, class Int=long>
class Graph1d
{
protected:
  KeyType delta, minX, maxX;
  Int N;
  std::vector<ValueType> y;
public:
  Graph1d() : delta(0.1), N(0), y(0), minX(std::numeric_limits<KeyType>::max()) {;}
  Graph1d(KeyType _delta) : delta(_delta), N(0), y(0), minX(std::numeric_limits<KeyType>::max()) {;}
  void setRangeAndClear(KeyType min, KeyType max, Int _N) {
    N=_N;
    y.resize(N);
    minX=min; maxX=max; // minY=maxY=ValueType(0);
    delta=(maxX-minX)/KeyType(N);
  }
  void setRange(KeyType min, KeyType max) {
    minX=min; maxX=max; // minY=maxY=ValueType(0);
    delta=(maxX-minX)/KeyType(N);
  }
  void setDeltaAndClear(KeyType _delta) {
    N=0; delta=_delta;
  }
  // inline Int idx(KeyType x) {return (x<minX) ? Int(std::floor((x-minX)/delta)) : Int((x-minX)/delta);}
  inline Int idx(KeyType x) {return (N==0) ? -1 : Int(N*(x-minX)/(maxX-minX)); }  
  inline KeyType keyFromIdx(Int i) { return minX+(KeyType(i)+KeyType(0.5))*delta;}
  inline ValueType &operator[](Int i) {
//  if(i<0 || i>=N)
//    {std::cerr<<"Graph1d index out of range:"<<i<<std::endl; exit(1);}
    return y[i];}
  inline void extendTo(KeyType x) {
    Int i;
    if(N<1)
    {
	N=1;
	minX=x-0.5*delta;
	maxX=minX+delta;
        y.resize(N);
    }
    else
    {
	i=idx(x);
	if(i<0)
	{
            y.resize(N-i);
	    for(Int j=N-1; j>=0; j--)
            {
              y[j-i]=y[j];
              y[j]=ValueType(0);
            }
	    N=N-i; minX+=KeyType(i)*delta;
	    i=0;
	}
        else if(i>=N)
	{
	    Int n=i-N+1;
            y.resize(i+1,ValueType(0));
	    N=i+1; maxX+=KeyType(n)*delta;
        }
    } 
  }
  inline ValueType &operator()(KeyType x) {
    if(N==0 || x<minX || x> maxX) extendTo(x);
    return y[idx(x)];}
  inline KeyType getDelta() {return delta;}
  inline KeyType getMinX() {return minX;}
  inline KeyType getMaxX() {return maxX;}
  inline Int getN() {return N;}
  inline ValueType getMinY() { ValueType h=std::numeric_limits<ValueType>::max();
    for(Int i=0; i<N; i++) if(y[i]<h) h=y[i]; return h;}
  inline ValueType getMaxY() { ValueType h=std::numeric_limits<ValueType>::min();
    for(Int i=0; i<N; i++) if(y[i]>h) h=y[i]; return h;}
  inline void getMinMaxY(ValueType &hMin, ValueType &hMax) {
    hMin=std::numeric_limits<ValueType>::max();
    hMax=std::numeric_limits<ValueType>::min();
    for(Int i=0; i<N; i++)
    {
      if(y[i]<hMin) hMin=y[i];
      if(y[i]>hMax) hMax=y[i];
    } }
  inline ValueType getMeanY() {
    ValueType h=ValueType(0);
    for(Int i=0; i<N; i++) h+=y[i]; return h/ValueType(N);}
  inline ValueType getMinYInInterval(KeyType b, KeyType t) {
    ValueType h=std::numeric_limits<ValueType>::max();
    for(Int i=idx(b); i<=idx(t); i++) if(y[i]<h) h=y[i]; return h;}
  inline ValueType getMaxYInInterval(KeyType b, KeyType t) {
    ValueType h=std::numeric_limits<ValueType>::min();
    for(Int i=idx(b); i<=idx(t); i++) if(y[i]>h) h=y[i]; return h;}
  inline void getMinMaxYInInterval(KeyType b, KeyType t, ValueType &hMin, ValueType &hMax) {
    hMin=std::numeric_limits<ValueType>::max();
    hMax=std::numeric_limits<ValueType>::min();
    for(Int i=idx(b); i<=idx(t); i++)
    {
      if(y[i]<hMin) hMin=y[i];
      if(y[i]>hMax) hMax=y[i];
    } }
  inline ValueType getMeanYInInterval(KeyType b, KeyType t) {
    ValueType h=ValueType(0);
    for(Int i=idx(b); i<=idx(t); i++) h+=y[i]; return h/ValueType(idx(t)-idx(b)+1);}
  void scale(ValueType s) {for(Int i=0; i<N; i++) y[i]*=s;}
  void clear() {for(Int i=0; i<N; i++) y[i]=ValueType(0);}
};

template<class ValueType, class KeyType, class Int>
inline bool syncronizeGraphs(Graph1d<ValueType, KeyType, Int> &a, Graph1d<ValueType, KeyType, Int> &b)
{
  if(a.getDelta()!=b.getDelta()) return false; // we cant sync if the intervals are different
  if(a.getMinX()<b.getMinX()) b.extendTo(a.getMinX()); else a.extendTo(b.getMinX());
  if(a.getMaxX()>b.getMaxX()) b.extendTo(a.getMaxX()); else a.extendTo(b.getMaxX());
  if(a.getN()!=b.getN()) return false; // should not happen!
  return true;
}

template<class ValueType, class KeyType, class Int=long>
class Kernel1d: public Graph1d<ValueType,KeyType,Int>
{
private:
  KeyType width;
  Int center;
public:
  Kernel1d() {;}
  Kernel1d(KeyType _delta, KeyType _width) : Graph1d<ValueType,KeyType,Int>(_delta) {
    width=_width;
    this->extendTo(KeyType(0));
    this->extendTo(width); this->extendTo(-width); center=this->idx(KeyType(0));}
  void setWidthAndClear(KeyType _delta, KeyType _width) {
    this->setDeltaAndClear(_delta);
    width=_width;
    this->extendTo(KeyType(0));
    this->extendTo(width); this->extendTo(-width); center=this->idx(KeyType(0));}
  KeyType getWidth() {return width;}
  Int getCenter() {return center;}
};

typedef enum {None, Epanechnikov, Quartic, TriWight, Triangle, Uniform, Gaussian, Cosine} KernelType;

template<class ValueType, class KeyType>
inline ValueType epanechnikov(KeyType x)
{
  ValueType h=ValueType(1)-ValueType(x*x);
  return h<0 ? ValueType(0) : ValueType(0.75)*h;
}

template<class ValueType, class KeyType>
inline ValueType quartic(KeyType x)
{
  ValueType h=ValueType(1)-ValueType(x*x);
  return ValueType(0.9375)*h*h; // 0.9375 = 15/16
}

template<class ValueType, class KeyType>
inline ValueType triwight(KeyType x)
{
  ValueType h=ValueType(1)-ValueType(x*x);
  return ValueType(1.09375)*h*h*h; // 1.09375 = 35/32
}

template<class ValueType, class KeyType>
inline ValueType triangle(KeyType x)
{
  return ValueType(1)-std::abs(ValueType(x));
}

template<class ValueType, class KeyType>
inline ValueType uniform(KeyType x)
{
  return std::abs(x)<1 ? ValueType(0.5) : ValueType(0);
}

template<class ValueType, class KeyType>
inline ValueType gaussian(KeyType x)
{
  return ValueType(std::sqrt(1.0/(2.0*M_PI)))*std::exp(ValueType(-0.5*x*x));
}

template<class ValueType, class KeyType>
inline ValueType cosine(KeyType x)
{
  return ValueType(0.25*M_PI)*std::cos(ValueType(0.5*M_PI*x));
}

template<class ValueType, class KeyType, class Int>
void initEpanechnikov(Kernel1d<ValueType, KeyType, Int> &k)
{
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*epanechnikov<ValueType,KeyType>(l*k.keyFromIdx(i));
}

template<class ValueType, class KeyType, class Int>
void initQuartic(Kernel1d<ValueType, KeyType, Int> &k)
{
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*quartic<ValueType,KeyType>(l*k.keyFromIdx(i));
}

template<class ValueType, class KeyType, class Int>
void initTriWight(Kernel1d<ValueType, KeyType, Int> &k)
{
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*triwight<ValueType,KeyType>(l*k.keyFromIdx(i));
}

template<class ValueType, class KeyType, class Int>
void initTriangle(Kernel1d<ValueType, KeyType, Int> &k)
{
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*triangle<ValueType,KeyType>(l*k.keyFromIdx(i));
}

template<class ValueType, class KeyType, class Int>
void initUniform(Kernel1d<ValueType, KeyType, Int> &k)
{
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*uniform<ValueType,KeyType>(l*k.keyFromIdx(i));
}

template<class ValueType, class KeyType, class Int>
void initGaussian(Kernel1d<ValueType, KeyType, Int> &k)
{
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*gaussian<ValueType,KeyType>(l*k.keyFromIdx(i));
}

template<class ValueType, class KeyType, class Int>
void initCosine(Kernel1d<ValueType, KeyType, Int> &k)
{
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*cosine<ValueType,KeyType>(l*k.keyFromIdx(i));
}

template<class ValueType, class KeyType, class Int>
void initNone(Kernel1d<ValueType, KeyType, Int> &k)
{
  k.setWidthAndClear(k.getDelta(),0.0);
  k[0]=1.0;
//  std::cout<<"Initialized 'None' kernel:\n  N="<<k.getN()
//           <<"\n  delta="<<k.getDelta()<<"\n  width="<<k.getWidth()
//           <<"\n  center="<<k.getCenter()<<std::endl;
}

template<class ValueType, class KeyType, class Int>
void initKernel(KernelType t, Kernel1d<ValueType, KeyType, Int> &k)
{
  switch(t)
  {
    case None: initNone(k); break;
    case Epanechnikov: initEpanechnikov(k); break;
    case Quartic: initQuartic(k); break;
    case TriWight: initTriWight(k); break;
    case Triangle: initTriangle(k); break;
    case Uniform: initUniform(k); break;
    case Gaussian: initGaussian(k); break;
    case Cosine: initCosine(k); break;
    default: std::cerr<<"Unknown kernel type in Graph1d!\n"; exit(1);
  }
}

void getKernelName(KernelType t, std::string &n)
{
  switch(t)
  {
    case None: n="None"; break;
    case Epanechnikov: n="Epanechnikov"; break;
    case Quartic: n="Quartic"; break;
    case TriWight: n="TriWight"; break;
    case Triangle: n="Triangle"; break;
    case Uniform: n="Uniform"; break;
    case Gaussian: n="Gaussian"; break;
    case Cosine: n="Cosine"; break;
    default: std::cerr<<"Unknown kernel type in Graph1d!\n"; exit(1);
  }
}

KernelType getKernelType(std::string &n)
{
  if(n=="None") return None;
  else if(n=="Epanechnikov") return Epanechnikov;
  else if(n=="Quartic") return Quartic;
  else if(n=="TriWight") return TriWight;
  else if(n=="Triangle") return Triangle;
  else if(n=="Uniform") return Uniform;
  else if(n=="Gaussian") return Gaussian;
  else if(n=="Cosine") return Cosine;
  return Epanechnikov;
}

template<class ValueType, class KeyType, class Int>
inline void addKernel(Graph1d<ValueType, KeyType, Int> &g, Kernel1d<ValueType, KeyType, Int> &k, KeyType x)
{
  //g.extendTo(x-k.getWidth()); g.extendTo(x+k.getWidth());
  Int i0 = g.idx(x-k.getWidth());
  for(Int i=0; i<k.getN(); i++) g[i+i0]+=k[i];
}

#endif
