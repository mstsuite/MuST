// -*- mode: c++ -*-
// Graph1dMoments class for 1d Wang-Landau Histogramm and DOS
// multiple linked 1d Graphs:
// [y, <m^1>, ... , <m^k>, M](x)

#ifndef WL_GRAPH1DMOMENTS_H
#define WL_GRAPH1DMOMENTS_H

#include "Graph1d.hpp"

#include <limits>
#include <cmath>
#include <stdlib.h>
#include <vector>

template <class ValueType, class KeyType, class Int=long>
class Graph1dMoments : public Graph1d<ValueType, KeyType, Int>
{
protected:
  int k;
  std::vector<Int> M;
  std::vector<std::vector<ValueType> > m;
public:
  Graph1dMoments() : k(0) {;}
  Graph1dMoments(KeyType _delta) : Graph1d<ValueType, KeyType, Int>(_delta), k(0) {;}

  bool test(void)
  {
    // check sizes
    bool notOK=false;
    size_t sN,sy,sM,sm;
    sN=Graph1d<ValueType, KeyType, Int>::N;
    sy=Graph1d<ValueType, KeyType, Int>::y.size();
    if(k!=m.size()) {printf("!!!!!!! Graph1dMoments::m size differs from k m:%d k:%d\n",m.size(),k); notOK=true;}
    sM=sm=sN;
    if(k>0)
    {
      sM=M.size();
      sm=m[0].size();
      for(int i=1;i<k;i++)
      {
        if(sm!=m[i].size()) {printf("!!!!!! Graph1dMoments::m of different sizes!!!!!\n"); notOK=true;}
      }
    }
    if(sy!=sN) {printf("!!!!!!! Graph1dMoments::y size differs y:%d N:%d\n",sy,sN); notOK=true;}
    if(sM!=sN) {printf("!!!!!!! Graph1dMoments::M size differs M:%d N:%d\n",sM,sN); notOK=true;}
    if(sm!=sN) {printf("!!!!!!! Graph1dMoments::m[*] sizes differ m[*]:%d N:%d\n",sm,sN); notOK=true;}
    return notOK;
  }

  void setRangeAndClear(KeyType min, KeyType max, Int _Nval) {
    Graph1d<ValueType, KeyType, Int>::setRangeAndClear(min,max,_Nval);
    if(k>0)
    {
      M.resize(_Nval);
      for(int i=0; i<k; i++) m[i].resize(_Nval);
    }
  }
  void setNumberOfMoments(int _k)
  {
    if(k!=0) {std::cerr<<"Graph1dMoments number of moments can be set only once!"<<std::endl; exit(1);}
    k=_k;
    m.resize(k);
    M.resize(Graph1d<ValueType, KeyType, Int>::N);
    for(int i=0; i<k; i++) m[i].resize(Graph1d<ValueType, KeyType, Int>::N);
  }
  int getNumberOfMoments() {return k;}
  inline void addMomentsAtIdx(Int i, ValueType v)
  {
    if(k>0)
    {
      ValueType vv=v;
      M[i]++;
      for(int ii=0; ii<k; ii++)
      {
        // if(i>=m[ii].size()) {printf("!!!! too large index in m[%d]:%d i:%d\n",ii,m[ii].size(),i); exit(1);}
        m[ii][i]+=(vv-m[ii][i])/((double)M[i]);
        vv*=v;
      }
    }
  }
  inline ValueType getMomentAtIdx(Int i, int _k)
  {
    return m[_k][i];
  }
  inline Int getNumberOfSamplesAtIdx(Int i)
  {
    return M[i];
  }
  inline void setMomentAtIdx(Int i, int _k, ValueType v)
  {
    // if(_k>=k) {printf("!!!!!! too large moment index!! _k:%d k:%d\n",_k,k); exit(1);}
    m[_k][i]=v;
  }
  inline void setNumberOfSamplesAtIdx(Int i, Int _M)
  {
    M[i]=_M;
  }

 inline void extendTo(KeyType x) {
    Int i;
    if(Graph1d<ValueType, KeyType, Int>::N<1)
      {
	Graph1d<ValueType, KeyType, Int>::N=1;
	Graph1d<ValueType, KeyType, Int>::minX=x-0.5*Graph1d<ValueType, KeyType, Int>::delta;
	Graph1d<ValueType, KeyType, Int>::maxX=Graph1d<ValueType, KeyType, Int>::minX+Graph1d<ValueType, KeyType, Int>::delta;
        Graph1d<ValueType, KeyType, Int>::y.resize(Graph1d<ValueType, KeyType, Int>::N);
        if(k>0)
        {
          M.resize(Graph1d<ValueType, KeyType, Int>::N);
          for(int i=0; i<k; i++) m[i].resize(Graph1d<ValueType, KeyType, Int>::N);
        }
      } else
      {
	i=this->idx(x);
	if(i<0)
	  {
            Graph1d<ValueType, KeyType, Int>::y.resize(Graph1d<ValueType, KeyType, Int>::N-i);
            if(k>0)
            {
              M.resize(Graph1d<ValueType, KeyType, Int>::N-i);
              for(int ii=0; ii<k; ii++) m[ii].resize(Graph1d<ValueType, KeyType, Int>::N-i);
            }
	    for(Int j=Graph1d<ValueType, KeyType, Int>::N-1; j>=0; j--)
            {
              Graph1d<ValueType, KeyType, Int>::y[j-i]=Graph1d<ValueType, KeyType, Int>::y[j];
              Graph1d<ValueType, KeyType, Int>::y[j]=ValueType(0);
              if(k>0)
              {
                M[j-i]=M[j];
                M[j]=0;
                for(int ii=0; ii<k; ii++)
                {
                  m[ii][j-i]=m[ii][j];
                  m[ii][j]=ValueType(0);
                }
              }
            }
	    Graph1d<ValueType, KeyType, Int>::N=Graph1d<ValueType, KeyType, Int>::N-i; Graph1d<ValueType, KeyType, Int>::minX+=KeyType(i)*Graph1d<ValueType, KeyType, Int>::delta;
	    i=0;
	  } else if(i>=Graph1d<ValueType, KeyType, Int>::N)
	  {
	    Int n=i-Graph1d<ValueType, KeyType, Int>::N+1;
            Graph1d<ValueType, KeyType, Int>::y.resize(i+1,ValueType(0));
            if(k>0)
            {
              M.resize(i+1,0);
              for(int ii=0; ii<k; ii++) m[ii].resize(i+1,ValueType(0));
            }
	    Graph1d<ValueType, KeyType, Int>::N=i+1;
            Graph1d<ValueType, KeyType, Int>::maxX+=KeyType(n)*Graph1d<ValueType, KeyType, Int>::delta;
	  }
      } }


  inline ValueType &operator()(KeyType x) {
    if(Graph1d<ValueType, KeyType, Int>::N==0 || x<Graph1d<ValueType, KeyType, Int>::minX || x> Graph1d<ValueType, KeyType, Int>::maxX) extendTo(x);
    return Graph1d<ValueType, KeyType, Int>::y[Graph1d<ValueType, KeyType, Int>::idx(x)];}
};

template<class ValueType, class KeyType, class Int>
inline bool syncronizeGraphs(Graph1dMoments<ValueType, KeyType, Int> &a, Graph1dMoments<ValueType, KeyType, Int> &b)
{
  if(a.getDelta()!=b.getDelta()) return false; // we cant sync if the intervals are different
  if(a.getMinX()<b.getMinX()) b.extendTo(a.getMinX()); else a.extendTo(b.getMinX());
  if(a.getMaxX()>b.getMaxX()) b.extendTo(a.getMaxX()); else a.extendTo(b.getMaxX());
  if(a.getN()!=b.getN()) return false; // should not happen!
  return true;
}

template<class ValueType, class KeyType, class Int>
inline void addKernel(Graph1dMoments<ValueType, KeyType, Int> &g, Kernel1d<ValueType, KeyType, Int> &k, KeyType x)
{
  //g.extendTo(x-k.getWidth()); g.extendTo(x+k.getWidth());
  Int i0 = g.idx(x-k.getWidth());
  for(Int i=0; i<k.getN(); i++) g[i+i0]+=k[i];
}

#endif
