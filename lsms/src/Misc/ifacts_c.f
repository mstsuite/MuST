c
c     ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ifacts(lmax,illp,ilp1,iprint,istop)
c     =================================================================
c
      implicit   none 
c
      character  sname*20
      character  istop*32
c
      integer    iprint
      integer    lmax
      integer    l
      integer    m
      integer    lp
      integer    mp
      integer    lm
      integer    lmp
c
      complex*16 illp((lmax+1)**2,(lmax+1)**2)
      complex*16 ilp1(0:2*lmax)
      complex*16 cone
      complex*16 sqrtm1
c
      parameter (sqrtm1=(0.0d0,1.0d0))
      parameter (cone=(1.0d0,0.0d0))
      parameter (sname='ifacts')
c
c     ================================================================
c     set up factors i**l.............................................
      ilp1(0)=cone
      do l=1,2*lmax
         ilp1(l)=ilp1(l-1)*sqrtm1
      enddo
c     ================================================================
c     load the factors i**(l-lp)into the matrix illp(lm,lmp).........
      do l=0,lmax
         do m= -l,l
            lm=l*(l+1)+1+m
            do lp=0,lmax
               do mp= -lp,lp
                  lmp=lp*(lp+1)+1+mp
                  illp(lm,lmp)=ilp1(l)/ilp1(lp)
c                 ====================================================
c                 printout if needed..................................
                  if(iprint.ge.1) then
                     write(6,'('' l,m,lp,mp,lm,lmp,illp(lm,lmp):'',
     >               6i4,2f8.3)') l,m,lp,mp,lm,lmp,illp(lm,lmp)
                  endif
               enddo
            enddo
         enddo
      enddo
c     ================================================================
c     set up factors i**(l+1)..........................................
      do l=0,2*lmax
         ilp1(l)=ilp1(l)*sqrtm1
      enddo
c     ================================================================
c     printout if needed..............................................
      if(iprint.ge.1) then
         write(6,'('' ifacts:: l,i**(l+1):'',i3,2d12.4)') 
     >   (l,ilp1(l),l=0,2*lmax)
      endif
c
c     ================================================================
      if(istop.eq.sname) then
         call fstop(sname)
      else
         return
      endif
      end
