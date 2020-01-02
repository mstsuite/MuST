      subroutine plglmax(lmax,x,plm)
c     ****************************************************************
c     Associated Legendre function....................................
c     Calclates all the p(l,m)'s up to lmax...........................
c     based on the formulae given in "Numerical Recipes" pages 180-183
c     (Equations 6.6.7, 6.6.8 and 6.6.9...............................
c     W. H. Press, B. P. Flannery, S A Teukolsky and W. T. Vetterling.
c     Cambridge Univ Press 1986.......................................
c     ****************************************************************
c
c Inputs:   lmax    integer scalar, max l for Plm
c           x       real*8 scalar, argument of Plm(x)
c Returns:  plm     real*8 array of ((lmax+1)*(lmax+2)/2), Plm(x)
c
      implicit   none
c
      integer    lmax
      integer    l,ll
      integer    m,mm
      integer    i
c
      real*8     plm((lmax+1)*(lmax+2)/2)
      real*8     x
      real*8     pmm
      real*8     somx2
      real*8     fact
      real*8     zero
      real*8     one
      real*8     tol
      real*8     two
c
      parameter  (zero=0.0d0)
      parameter  (one=1.0d0)
      parameter  (two=2.0d0)
      parameter  (tol=1.d-12)
c
c     ================================================================
      if(lmax.lt.0 .or. abs(x).gt.(one+tol)) then
         write(6,'(''plglmax:: bad arguments: lmax='',i5,'' x='',
     >         d14.6)') lmax,x
         stop'plglmax'
      endif
!      write(6,*) "plglmax: x=",x," lmax=",lmax
c
      if((one-abs(x)).le.tol) then
c        -------------------------------------------------------------
         call zeroout(plm,(lmax+1)*(lmax+2)/2)
c        -------------------------------------------------------------
	 if(x.lt.zero) then
            do l=0,lmax
	       i=(l+1)*(l+2)/2-l
	       plm(i)=one-two*mod(l,2)
            enddo
         else
            do l=0,lmax
	       i=(l+1)*(l+2)/2-l
	       plm(i)=one
            enddo
         endif
         return
      endif
c
c     ================================================================
c     begin calculation of p(l,m)'s...................................
      if(lmax.eq.0) then
c        =============================================================
c        special case lmax=0..........................................
         plm(1)=one
      else
c        =============================================================
c        minus sign added to be consistant with Numerical Recipes
c        which has (-1)^m factor in plm :  July 97  by xgz.............
c        =============================================================
         somx2=-sqrt((one-x)*(one+x))
         if(lmax.eq.1) then
c           ==========================================================
c           special case lmax=1.......................................
            plm(1)=one
            plm(2)=x
            plm(3)=somx2
         else
            do m=0,lmax
c              =======================================================
c                                       m       m
c              calculate the first two P   and P
c                                       m       m+1
c              =======================================================
               if(m.eq.0) then
                  plm(1)=one
                  plm(2)=x
               else
                  pmm=somx2
                  fact=one
                  do i=2,m
                     fact=fact+two
                     pmm=pmm*fact*somx2
                  enddo
		  mm=(m+1)*(m+2)/2
                  plm(mm)=pmm
                  if( mm+m+1.le.(lmax+1)*(lmax+2)/2 ) then
                     plm(mm+m+1)=x*(2*m+1)*pmm
                  end if
               endif
c              =======================================================
c                                  m        m
c              calculate the rest P     to P
c                                  m+2      lmax
c              =======================================================
               ll=(m+2)*(m+1)/2
               fact=(two*m+one)*x
               do l=m+2,lmax
                  pmm=(l+m-1)*plm(ll)
                  fact=fact+two*x
                  ll=ll+l-1
                  plm(ll+l)=( fact*plm(ll) - pmm )/dble(l-m)
               enddo
            enddo
         endif
      endif
c
      return
      end
