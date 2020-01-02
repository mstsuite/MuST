// -*- mode: c++ -*-
// Graph2d class for 2d Wang-Landau Histogramm and DOS

#ifndef WL_GRAPH2D_H
#define WL_GRAPH2D_H

#include <limits>
#include <cmath>
#include <stdlib.h>
#include <vector>


template <class ValueType, class KeyType, class Int=long>
class Graph2d
{
private:
  KeyType deltaX, minX, maxX, deltaY;
  Int Nx;
  std::vector< std::vector<ValueType>* > val;
  std::vector<KeyType> minY, maxY;
  std::vector<Int> Ny;
public:
  Graph2d() : deltaX(0.1), deltaY(0.1), Nx(0) {;}
  Graph2d(KeyType _deltaX, KeyType _deltaY) : deltaX(_deltaX), deltaY(_deltaY), Nx(0) {;}
  ~Graph2d() {freeVal(); }
  void freeVal() {for(Int i=0; i<val.size(); i++) {delete val[i]; val[i]=NULL; Ny[i]=0;} }
  void setXRangeAndClear(KeyType min, KeyType max, Int _Nval) {
    Nx=_Nval;
    freeVal();
    val.resize(Nx); Ny.resize(Nx);
    minX=min; maxX=max; // minY=maxY=ValueType(0);
    deltaX=(maxX-minX)/KeyType(Nx);
  }
  void setXRange(KeyType min, KeyType max) {
    minX=min; maxX=max; // minY=maxY=ValueType(0);
    deltaX=(maxX-minX)/KeyType(Nx);
  }
  void setDeltaAndClear(KeyType _deltaX, KeyType _deltaY) {
    freeVal();
    Nx=0;
    deltaX=_deltaX;
    deltaY=_deltaY;
  }
  void setYminAndClear(Int ix, KeyType _minY, Int _Nyval)
  {
    val[ix] = new std::vector<ValueType>(_Nyval);
    minY[ix]= _minY;
    maxY[ix]= _minY + deltaY*KeyType(_Nyval);
  }
  inline Int idxX(KeyType x)
  {
    return (x<minX) ? Int(std::floor((x-minX)/deltaX)) :
      Int((x-minX)/deltaX);
  }

  inline Int idxY(Int ix, KeyType y)
  {
    return (y<minY[ix]) ? Int(std::floor((y-minY[ix])/deltaY)) :
      Int((y-minY[ix])/deltaY);
  }

  inline KeyType keyFromIdxX(Int ix)
  {
    return minX+(KeyType(ix)+KeyType(0.5))*deltaX;
  }
  inline KeyType keyFromIdxY(Int ix, Int iy)
  {
    return minY[ix]+(KeyType(iy)+KeyType(0.5))*deltaY;
  }

  inline std::vector<ValueType> &operator[](Int ix) {
//  if(i<0 || i>=N)
//    {std::cerr<<"Graph2d index out of range:"<<i<<std::endl; exit(1);}
    return *(val[ix]);}

  inline Int extendToX(KeyType x)
  {
    Int ix;
    if(Nx<1)
      {
	Nx=1;
        KeyType h = std::floor(x/deltaX+0.5);
	minX=(h-0.5)*deltaX;
	maxX=minX+deltaX;
        val.resize(Nx); Ny.resize(Nx);
	minY.resize(Nx); maxY.resize(Nx);
	ix=0;
      } else
      {
	ix=idxX(x);
	if(ix<0)
	  {
            val.resize(Nx-ix); Ny.resize(Nx-ix);
	    minY.resize(Nx-ix); maxY.resize(Nx-ix);
	    for(Int j=Nx-1; j>=0; j--)
	      {
		val[j-ix]=val[j]; val[j]=NULL;
		Ny[j-ix]=Ny[j]; Ny[j]=0;
		minY[j-ix]=minY[j]; maxY[j-ix]=maxY[j];
	      }
	    Nx=Nx-ix; minX+=KeyType(ix)*deltaX;
	    ix=0;
	  } else if(ix>=Nx)
	  {
	    Int n=ix-Nx+1;
            val.resize(ix+1,NULL);
	    Ny.resize(ix+1,0);
	    minY.resize(ix+1); maxY.resize(ix+1);
	    Nx=ix+1; maxX+=KeyType(n)*deltaX;
	  }
      }
    return ix;
  }

  inline void extendToY(Int ix, KeyType y) {
    Int iy;

    if(Ny[ix]<1)
      {
	Ny[ix]=1;
        KeyType h = std::floor(y/deltaY+0.5);
	minY[ix]=(h-0.5)*deltaY;
	maxY[ix]=minY[ix]+deltaY;
	if(val[ix]==NULL) val[ix]=new std::vector<ValueType>(Ny[ix]);
        else val[ix]->resize(Ny[ix]);
	iy=0;
      } else
      {
	iy=idxY(ix,y);
	if(iy<0)
	  {
            val[ix]->resize(Ny[ix]-iy);
	    for(Int j=Ny[ix]-1; j>=0; j--)
	      { (*val[ix])[j-iy]=(*val[ix])[j]; (*val[ix])[j]=ValueType(0); }
	    Ny[ix]=Ny[ix]-iy; minY[ix]+=KeyType(iy)*deltaY;
	    iy=0;
	  } else if(iy>=Ny[ix])
	  {
	    Int n=iy-Ny[ix]+1;
            val[ix]->resize(iy+1,ValueType(0));
	    Ny[ix]=iy+1; maxY[ix]+=KeyType(n)*deltaY;
	  }
      }
  }

  inline void extendTo(KeyType x, KeyType y) {
    Int ix=extendToX(x);
    extendToY(ix, y);
  }

  inline ValueType &operator()(KeyType x , KeyType y) {
    Int ix = idxX(x);
    return (*val[ix])[idxY(ix,y)];}
  inline KeyType getDeltaX() {return deltaX;}
  inline KeyType getDeltaY() {return deltaY;}
  inline KeyType getMinX() {return minX;}
  inline KeyType getMinY(Int ix) {return minY[ix];}
  inline KeyType getMaxX() {return maxX;}
  inline KeyType getMaxY(Int ix) {return maxY[ix];}
  inline Int getNx() {return Nx;}
  inline Int getNy(Int ix) {return Ny[ix];}
  inline ValueType getMinVal() { ValueType h=std::numeric_limits<ValueType>::max();
    for(Int ix=0; ix<Nx; ix++) for(Int iy=0; iy<Ny[ix]; iy++) if((*val[ix])[iy]<h) h=(*val[ix])[iy]; return h;}
  inline ValueType getMaxVal() { ValueType h=std::numeric_limits<ValueType>::min();
    for(Int ix=0; ix<Nx; ix++) for(Int iy=0; iy<Ny[ix]; iy++) if((*val[ix])[iy]>h) h=(*val[ix])[iy]; return h;}
  inline ValueType getMinValWithBorders(KeyType xBorder, KeyType yBorder)
  {
    Int xStart=idxX(minX+xBorder);
    Int xEnd=idxX(maxX-xBorder);
    ValueType h=std::numeric_limits<ValueType>::max();
    for(Int ix=xStart; ix<xEnd; ix++)
    {
      Int yStart=idxY(ix,minY[ix]+yBorder);
      Int yEnd=idxY(ix,maxY[ix]-yBorder);
      for(Int iy=yStart; iy<yEnd; iy++)
        if((*val[ix])[iy]<h) h=(*val[ix])[iy];
    }
    return h;
  }
  inline ValueType getMaxValWithBorders(KeyType xBorder, KeyType yBorder)
  {
    Int xStart=idxX(minX+xBorder);
    Int xEnd=idxX(maxX-xBorder);
    ValueType h=std::numeric_limits<ValueType>::max();
    for(Int ix=xStart; ix<xEnd; ix++)
    {
      Int yStart=idxY(ix,minY[ix]+yBorder);
      Int yEnd=idxY(ix,maxY[ix]-yBorder);
      for(Int iy=yStart; iy<yEnd; iy++)
        if((*val[ix])[iy]>h) h=(*val[ix])[iy];
    }
    return h;
  }
  void scale(ValueType s)
  {
    for(Int ix=0; ix<Nx; ix++)
      {for(Int iy=0; iy<Ny[ix]; iy++) (*val[ix])[iy]*=s;}
  }
  void clear()
  {
    for(Int ix=0; ix<Nx; ix++)
      {for(Int iy=0; iy<Ny[ix]; iy++) (*val[ix])[iy]=ValueType(0);}
  }
};

template<class ValueType, class KeyType, class Int>
inline bool syncronizeGraphs(Graph2d<ValueType, KeyType, Int> &a, Graph2d<ValueType, KeyType, Int> &b)
{
  if(a.getDelta()!=b.getDelta()) return false; // we cant sync if the intervals are different
  if(a.getMinX()<b.getMinX()) b.extendTo(a.getMinX()); else a.extendTo(b.getMinX());
  if(a.getMaxX()>b.getMaxX()) b.extendTo(a.getMaxX()); else a.extendTo(b.getMaxX());
  if(a.getN()!=b.getN()) return false; // should not happen!
  return true;
}

template<class ValueType, class KeyType, class Int=long>
class Kernel2d: public Graph2d<ValueType,KeyType,Int>
{
private:
  KeyType widthX, widthY;
  Int centerX, centerY;
public:
  Kernel2d() {;}
  Kernel2d(KeyType _deltaX, KeyType _deltaY, KeyType _widthX, KeyType _widthY)
    : Graph2d<ValueType,KeyType,Int>(_deltaX, _deltaY), widthX(_widthX), widthY(_widthY) {
    this->extendTo(KeyType(0), KeyType(0));
    /*
    this->extendTo(widthX,widthY); this->extendTo(widthX,-widthY);
    this->extendTo(-widthX, widthY); this->extendTo(-widthX,-widthY);
    */
    this->extendToX(widthX); this->extendToX(-widthX);
    for(Int ix=0; ix<this->getNx(); ix++)
      {
	this->extendToY(ix,widthY); this->extendToY(ix,-widthY);
      }
    centerX=this->idxX(KeyType(0)); centerY=this->idxY(centerX,KeyType(0));
  }
  void setWidthAndClear(KeyType _deltaX, KeyType _deltaY, KeyType _widthX, KeyType _widthY) {
    this->setDeltaAndClear(_deltaX, _deltaY);
    widthX=_widthX; widthY=_widthY;
    this->extendTo(KeyType(0), KeyType(0));
    /*
    this->extendTo(widthX,widthY); this->extendTo(widthX,-widthY);
    this->extendTo(-widthX, widthY); this->extendTo(-widthX,-widthY);
    */
    this->extendToX(widthX); this->extendToX(-widthX);
    for(Int ix=0; ix<this->getNx(); ix++)
      {
	this->extendToY(ix,widthY); this->extendToY(ix,-widthY);
      }
    centerX=-this->idxX(KeyType(0)); centerY=-this->idxY(centerX,KeyType(0)); }
  KeyType getWidthX() {return widthX;}
  KeyType getWidthY() {return widthY;}
  Int getCenterX() {return centerX;}
  Int getCenterY() {return centerY;}
};

// these are also defined in Graph1d
#if !defined(WL_GRAPH1D_H) && !defined(WL_GRAPH1DMOMENTS_H)
typedef enum {Epanechnikov, Quartic, TriWight, Triangle, Uniform, Gaussian, Cosine} KernelType;

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

void getKernelName(KernelType t, std::string &n)
{
  switch(t)
  {
    case Epanechnikov: n="Epanechnikov"; break;
    case Quartic: n="Quartic"; break;
    case TriWight: n="TriWight"; break;
    case Triangle: n="Triangle"; break;
    case Uniform: n="Uniform"; break;
    case Gaussian: n="Gaussian"; break;
    case Cosine: n="Cosine"; break;
    default: std::cerr<<"Unknown kernel type in Graph2d!\n"; exit(1);
  }
}

KernelType getKernelType(std::string &n)
{
  if(n=="Epanechnikov") return Epanechnikov;
  else if(n=="Quartic") return Quartic;
  else if(n=="TriWight") return TriWight;
  else if(n=="Triangle") return Triangle;
  else if(n=="Uniform") return Uniform;
  else if(n=="Gaussian") return Gaussian;
  else if(n=="Cosine") return Cosine;
  return Epanechnikov;
}

#endif

template<class ValueType, class KeyType, class Int>
void initKernelFromFunction(Kernel2d<ValueType, KeyType, Int> &k, ValueType (*f)(KeyType))
{
  KeyType lX=KeyType(1)/KeyType(k.getWidthX());
  KeyType lY=KeyType(1)/KeyType(k.getWidthY());
  for(Int ix=0; ix<k.getNx(); ix++)
    {
      for(Int iy=0; iy<k.getNy(ix); iy++)
	{
	  k[ix][iy]=(*f)(std::sqrt(lX*k.keyFromIdxX(ix)*lX*k.keyFromIdxX(ix)+
				   lY*k.keyFromIdxY(ix,iy)*lY*k.keyFromIdxY(ix,iy)));
	  // std::cout<<"k["<<ix<<"]["<<iy<<"]="<<k[ix][iy]<<std::endl;
	}
    }
}


template<class ValueType, class KeyType, class Int>
void initEpanechnikov(Kernel2d<ValueType, KeyType, Int> &k)
{
  initKernelFromFunction(k,&epanechnikov<ValueType,KeyType>);
  /*
  KeyType lX=KeyType(1)/KeyType(k.getWidthX());
  KeyType lY=KeyType(1)/KeyType(k.getWidthY());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*epanechnikov<ValueType,KeyType>(l*k.keyFromIdx(i));
  */
}

template<class ValueType, class KeyType, class Int>
void initQuartic(Kernel2d<ValueType, KeyType, Int> &k)
{
  initKernelFromFunction(k,&quartic<ValueType,KeyType>);
  /*
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*quartic<ValueType,KeyType>(l*k.keyFromIdx(i));
  */
}

template<class ValueType, class KeyType, class Int>
void initTriWight(Kernel2d<ValueType, KeyType, Int> &k)
{
  initKernelFromFunction(k,&triwight<ValueType,KeyType>);
  /*
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*triwight<ValueType,KeyType>(l*k.keyFromIdx(i));
  */
}

template<class ValueType, class KeyType, class Int>
void initTriangle(Kernel2d<ValueType, KeyType, Int> &k)
{
  initKernelFromFunction(k,&triangle<ValueType,KeyType>);
  /*
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*triangle<ValueType,KeyType>(l*k.keyFromIdx(i));
  */
}

template<class ValueType, class KeyType, class Int>
void initUniform(Kernel2d<ValueType, KeyType, Int> &k)
{
  initKernelFromFunction(k,&uniform<ValueType,KeyType>);
  /*
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*uniform<ValueType,KeyType>(l*k.keyFromIdx(i));
  */
}

template<class ValueType, class KeyType, class Int>
void initGaussian(Kernel2d<ValueType, KeyType, Int> &k)
{
  initKernelFromFunction(k,&gaussian<ValueType,KeyType>);
  /*
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*gaussian<ValueType,KeyType>(l*k.keyFromIdx(i));
  */
}

template<class ValueType, class KeyType, class Int>
void initCosine(Kernel2d<ValueType, KeyType, Int> &k)
{
  initKernelFromFunction(k,&cosine<ValueType,KeyType>);
  /*
  KeyType l=KeyType(1)/KeyType(k.getWidth());
  for(Int i=0; i<k.getN(); i++) k[i]=ValueType(l)*cosine<ValueType,KeyType>(l*k.keyFromIdx(i));
  */
}

template<class ValueType, class KeyType, class Int>
void initKernel(KernelType t, Kernel2d<ValueType, KeyType, Int> &k)
{
  switch(t)
  {
    case Epanechnikov: initEpanechnikov(k); break;
    case Quartic: initQuartic(k); break;
    case TriWight: initTriWight(k); break;
    case Triangle: initTriangle(k); break;
    case Uniform: initUniform(k); break;
    case Gaussian: initGaussian(k); break;
    case Cosine: initCosine(k); break;
    default: std::cerr<<"Unknown kernel type in Graph2d!\n"; exit(1);
  }
}

template<class ValueType, class KeyType, class Int>
inline void addKernel(Graph2d<ValueType, KeyType, Int> &g, Kernel2d<ValueType, KeyType, Int> &k, KeyType x, KeyType y)
{
  g.extendToX(x-k.getWidthX()); g.extendToX(x+k.getWidthX());
  Int ix0 = g.idxX(x-k.getWidthX());
  for(Int ix=0; ix<k.getNx(); ix++)
    {
      g.extendToY(ix+ix0,y-k.getWidthY()); g.extendToY(ix+ix0,y+k.getWidthY());
      Int iy0 = g.idxY(ix+ix0,y-k.getWidthY());
      for(Int iy=0; iy<k.getNy(ix); iy++)
	{
	  // std::cout<<"ix0="<<ix0<<" iy0="<<iy0<<" ix="<<ix<<" iy="<<iy<<std::endl;
	  // std::cout<<"k.Nx="<<k.getNx()<<" g.Nx="<<g.getNx()<<std::endl;
	  // std::cout<<"k.Ny[ix]="<<k.getNy(ix)<<" g.Ny[ix+ix0]="<<g.getNy(ix+ix0)<<std::endl;
	  // std::cout<<"k[ix][iy]="<<k[ix][iy]<<std::endl;
	  g[ix+ix0][iy+iy0]+=k[ix][iy];
	}
    }
}

#endif
