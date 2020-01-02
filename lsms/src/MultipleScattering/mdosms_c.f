c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine mdosms(zj_flag,n,dos,zz,zj,w1,pi,
     >                  iprint,istop)
c     ================================================================
c
      implicit   none
c
      character  istop*32
      character  sname*20
      parameter (sname='mdosms')
c
      integer    zj_flag
      integer    n,kkrsz_loc
      integer    iprint
      integer    i
c
      real*8     pi
c
      complex*16 zz(n)
      complex*16 zj
      complex*16 w1(n)
      complex*16 dos
      complex*16 czero
      complex*16 ctmp
c
      parameter (czero=(0.0d0,0.0d0))
c
c     ****************************************************************
c     calculates the density of states................................
c     written by w.a.s.,jr sept. 1, 1992
c     ****************************************************************
c
      dos=czero
      kkrsz_loc=sqrt(dble(n))+.5d0
c     ----------------------------------------------------------------
c  Backward sum try to cancel singularities in high l's first.
      do i=1,n
	dos=dos+zz(i)*w1(i)
      enddo
c       l component dos but no ss terms
c       if(iprint.ge.100) then
c          write(6,'(i5,2f14.8,18f9.4,'' dosd'')') 
c    >     mynod,-dos/pi,
c    >   (-zz((i-1)*kkrsz_loc+i)*w1((i-1)*kkrsz_loc+i)/pi,i=1,min(9,n))
c       endif
	if(zj_flag.eq.1) then
	  dos=dos-zj
	endif
      dos=-dos/pi
c
c     ================================================================
      if(istop.eq.sname) then
         call fstop(sname)
      endif
c
      return
      end
