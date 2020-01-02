!     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function ylag(xi,x,y,ind1,n1,imax,iex) result(lag)
!     ================================================================
!
!     lagrangian interpolation
!   modified from subroutine polint from Numerical Recipe, by Press et al.
!     xi is intepolated entry into x-array
!     n is the order of lagrangran interpolation
!     y is array from which ylag is obtained by interpolation
!     ind is the min-i for x(i).gt.xi
!     if ind=0,x-array will be searched
!     imax is max index of x-and y-arrays
!     extrapolation can occur,iex=-1 or +1
!
!     ================================================================
      use KindParamModule, only : IntKind, RealKind
!
      implicit none
!
      integer (kind=IntKind), intent(in) :: ind1, n1, imax
      integer (kind=IntKind), intent(out) :: iex
!
      real (kind=RealKind), intent(in) :: xi, x(*), y(*)
!
      integer (kind=IntKind), parameter :: nmax=10
      integer (kind=IntKind) :: ind, n, j, inu, ns, inl, j1, i, ij
!
      real (kind=RealKind) :: lag, dift, dif, den
      real (kind=RealKind) :: c(nmax), d(nmax), dx(nmax)
!
      ind=ind1
      n=n1
      iex=0
      if (n.le.imax) go to 10
      n=imax
      iex=n
   10 if (ind.gt.0) go to 40
      do 20 j = 1,imax
!        if (xi-x(j)) 30,130,20
         if (xi-x(j) < 0.0d0) then
            goto 30
         else if (xi-x(j) > 0.0d0) then
            goto 20
         else
            goto 130
         endif
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
      do j=1,inu-inl
         j1=j+inl
         dx(j)=xi-x(j1)
         dift=abs(dx(j))
         if(dift.lt.dif) then
           ns=j
           dif=dift
         endif
         c(j)=y(j1)
         d(j)=y(j1)
      enddo
      lag=y(inl+ns)
      ns=ns-1
      do j=inl+1,inu-1
         do i=1,inu-j
            ij=i+j
            den=x(i+inl)-x(ij)
            if(den.eq.0.d0) stop 'Two xs are the same in ylag.'
            den=(d(i)-c(i+1))/den
            d(i)=dx(ij-inl)*den
            c(i)=dx(i)*den
         enddo
         if(ishft(ns,1).lt.inu-j) then
           lag=lag+c(2*ns+1)
         else
           lag=lag+d(ns)
           ns=ns-1
         endif
      enddo
      return
  130 lag=y(j)
      return
      end function ylag
