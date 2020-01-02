c...............................................
      function ylag(xi,x,y,ind1,n1,imax,iex)
c
c     lagrangian interpolation
c   modified from subroutine polint from Numerical Recipe, by Press et al.
c     xi is intepolated entry into x-array
c     n is the order of lagrangran interpolation
c     y is array from which ylag is obtained by interpolation
c     ind is the min-i for x(i).gt.xi
c     if ind=0,x-array will be searched
c     imax is max index of x-and y-arrays
c     extrapolation can occur,iex=-1 or +1
c
      implicit real*8(a-h,o-z)
      parameter(nmax=10)
      real*8 x(*),y(*)
      real*8 c(nmax),d(nmax),dx(nmax)
c
      ind=ind1
      n=n1
      iex=0
      if (n.le.imax) go to 10
      n=imax
      iex=n
   10 if (ind.gt.0) go to 40
      do 20 j = 1,imax
         if (xi-x(j) .lt. 0.0) goto 30
         if (xi-x(j) .eq. 0.0) goto 130
   20 continue
      iex=1
      go to 70
   30 ind=j
   40 if (ind.le.1) iex=-1
      inl=max(0,ind-ishft(n+1,-1))
      inu=inl+n
      if (inu.le.imax) go to 80
   70 inl=imax-n
      inu=imax
   80 dif=abs(xi-x(inl+1))
      ns=1
      do 105 j=1,inu-inl
	 j1=j+inl
	 dx(j)=xi-x(j1)
         dift=abs(dx(j))
         if(dift.lt.dif) then
           ns=j
           dif=dift
         endif
         c(j)=y(j1)
         d(j)=y(j1)
  105 continue
      ylag=y(inl+ns)
      ns=ns-1
      do 110 j=inl+1,inu-1
         do 100 i=1,inu-j
	    ij=i+j
            den=x(i+inl)-x(ij)
            if(den.eq.0.d0) stop'Two xs are the same in ylag.'
            den=(d(i)-c(i+1))/den
            d(i)=dx(ij-inl)*den
            c(i)=dx(i)*den
  100    continue
         if(ishft(ns,1).lt.inu-j) then
           ylag=ylag+c(2*ns+1)
         else
           ylag=ylag+d(ns)
           ns=ns-1
         endif
  110 continue
      return
  130 ylag=y(j)
      return
      end
