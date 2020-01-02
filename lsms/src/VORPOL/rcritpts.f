c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rcritpts(rcrit,ncrit,
     >                    rmt2,rcs2,
     >                    xp,nbnd,
     >                    dc2,ncorn,
     >                    edgp,nedge,
     >                    runion)
c     ================================================================
c
      implicit   none
c
      integer    ncrit
      integer    nbnd
      integer    ncorn
      integer    nedge
      integer    iunion
      integer    i
      integer    sigma
c
      real*8     dc2(ncorn)
      real*8     runion(nbnd+ncorn+nedge+1)
      real*8     rcrit(*)
      real*8     rtol
      real*8     rmt2,rcs2
      real*8     xp(3,nbnd)
      real*8     edgp(3,nedge)
c
      parameter (rtol=1.0d-10)
c
c     ================================================================
c     joint tables of possible critical points........................
      iunion=1
      runion(1)=rmt2
      do i=1,nbnd
        if(sigma(xp(1,i),xp(2,i),xp(3,i),xp,nbnd,1).eq.1) then
          iunion=iunion+1
	  runion(iunion)=xp(1,i)*xp(1,i)+xp(2,i)*xp(2,i)+xp(3,i)*xp(3,i)
	endif
      enddo
      do i=1,ncorn
         runion(iunion+i)=dc2(i)
      enddo
      iunion=iunion+ncorn
      do i=1,nedge
         if(sigma(edgp(1,i),edgp(2,i),edgp(3,i),xp,nbnd,1).eq.1)
     &     then
           iunion=iunion+1
	   runion(iunion)=
     &       edgp(1,i)*edgp(1,i)+edgp(2,i)*edgp(2,i)+edgp(3,i)*edgp(3,i)
	 endif
      enddo
c     ================================================================
c     sort the table of possible critical points......................
c     ----------------------------------------------------------------
      call sort(iunion,runion)
c     ----------------------------------------------------------------
c     ================================================================
c     obtain list of unique critical points...........................
      ncrit=1
      rcrit(1)=runion(1)
      do i=2,iunion
         if(abs(rcrit(ncrit)-runion(i)).gt.rtol .and.
     >      (runion(i)-rcs2).le.rtol) then
            ncrit=ncrit+1
            rcrit(ncrit)=runion(i)
         endif
      enddo
      do i=1,ncrit
         rcrit(i)=sqrt(rcrit(i))
      enddo
c
c     ================================================================
c
      return
      end
