c
c    ===========================================================
      subroutine madewdj(lastrs,lastkn,myatom,num_atoms,
     >                  rslat_x,rslat_y,rslat_z,rssq,
     &                  knlat_x,knlat_y,knlat_z,xknsq,
     >                  atom_posi_x,atom_posi_y,atom_posi_z,
     >                  lmax,clm,lofj,mofj,ndlm,ndlmv,rmau,rcutau,
     >                  eta,omegbra,madmat,iprint,istop)
c    ===========================================================
c
c    *******************************************************************
c    program calculates the madelung potential from charges in periodically
c    repeated cells excluding the central unit cell; the potential 
c    written by w.a.s.,jr. 2-21-88 modified to give potential from
c    only atoms outside the unit cell and for l> 0 by DMN 1993.
c
c    *******************************************************************
c     input :
c     =====
c         1.) ewald card
cccccc        lastrs = number real-space vectors
cccccc        rsn= real-space vectors a.u.
cccccc        atom_posi_? =basis vectors in cells a.u.
cccccc        lastkn = number reciprocal-space vectors
cccccc        knlat = reciprocal-space vectors a.u.
cccccc        eta = (ewald parameter) in a.u.
cccccc        num_atoms = number of sublattices
cccccc        myatom = the sublattice for which the column of the 
cccccc                madelung matrix is to be calculated
cccccc        omegbra = unit cell volume a.u.
cccccc        rmau = radius at which potential is calculated a.u.
cccccc        rcutau = radius within which contributions from lattice points
cccccc                  are excluded   a.u.
cccccc        pi = 3.141....
cccccc        sqrtm1 = (0.0,1.0)
c     output :
c     =====
cccccc        rssq= real-space vectors magnitude squared a.u.
cccccc        xknsq= reciprocal-space vectors magnitude squared a.u.
c     internal :
c     =====
cccccc        aijau = distance between sublattices 
c
c   calling sequence
c   ================
c    |
c  pqintg                       ! evaluates integral over eta
c    |
c    |
c  madsum                       ! real & recip space sums
c
c***********************************************************************
      implicit none
c
      include 'madelung.h'
c
      character istop*32
      character sname*32
      parameter (sname='madewdj')
c
      integer    ipndlm
      parameter  (ipndlm=(iplmax_mad+1)*(iplmax_mad+2)/2)
c
      integer    lastrs
      integer    lastkn
      integer    num_atoms
      integer    iprint
      integer    myatom
      integer    nsub2
      integer    i
      integer    icell
      integer    l
      integer    m
      integer    ibegin
      integer    j
      integer    ndlm,ndlmv
      integer    lmax
      integer    lofj(ndlm)
      integer    mofj(ndlm)
c
      real*8     r
      real*8     erfc
      real*8     rdif
      real*8     atom_posi_x(num_atoms)
      real*8     atom_posi_y(num_atoms)
      real*8     atom_posi_z(num_atoms)
      real*8     jl(iplmax_mad+1,ipknlat_mad),work(iplmax_mad+1)
      real*8     aijau(3)
      real*8     plmk(ipndlm,ipknlat_mad)
      real*8     dqint(iprslat_mad)
      real*8     dqintl(0:iplmax_mad)
      real*8     knlat_x(lastkn)
      real*8     knlat_y(lastkn)
      real*8     knlat_z(lastkn)
      real*8     rslat_x(lastrs)
      real*8     rslat_y(lastrs)
      real*8     rslat_z(lastrs)
      real*8     rssq(lastrs)
      real*8     rcutau
      real*8     rcut2
      real*8     xknsq(lastkn)
      real*8     xknabs(ipknlat_mad)
      real*8     xknabsi
      real*8     clm(ndlm)
      real*8     xgs(ngauss_mad)
      real*8     wgs(ngauss_mad)
      real*8     one
      real*8     four
      real*8     pi
      real*8     rtmp1
      real*8     rtmp2
      real*8     rtmp3
      real*8     xk,rk
      real*8     yk
      real*8     zk
      real*8     term
      real*8     rmau
      real*8     omegbra
      real*8     eta
      real*8     twoetasq
      real*8     phi
      real*8     plm(ipndlm)
      real*8     z
      real*8     fnpi
      real*8     zero
      real*8     half
      real*8     two
      real*8     eight
      parameter (zero=0.0d0)
      parameter (half=0.5d0)
      parameter (one=1.0d0)
      parameter (two=2.0d0)
      parameter (four=4.0d0)
      parameter (eight=8.0d0)
c
      real*8     tmps,tmpc
c
      complex*16 sum1(ipndlm)
      complex*16 czero
      complex*16 dot
      complex*16 expidot
      complex*16 facl(0:iplmax_mad)
      complex*16 madmat(ndlmv,num_atoms)
      complex*16 dqintloc(ipndlm)
      complex*16 sum2
      complex*16 tmp
      complex*16 trans
      complex*16 cone
      complex*16 expphi
      complex*16 expfac
      complex*16 cylm
c     complex*16 cylmij(ipknlat_mad,ipndlm)
      complex*16 fack(ipknlat_mad,ipndlm)
      complex*16 expphim(0:iplmax_mad)
      complex*16 sqrtm1
      parameter (cone=(1.0d0,0.0d0))
      parameter (czero=(0.0d0,0.0d0))
      parameter (sqrtm1=(0.0d0,1.0d0))
c
      pi=fnpi()
      twoetasq=(two*eta)**2
      rcut2=(rcutau)**2
c     get weights for integration over cos(theta). used for removing
c     contribution of gaussian charges inside rcut which are included in 
c     k space sum.
!     call gauleg(-one,one,xgs,wgs,ngauss_mad)
      call gauss_legendre_points(-one,one,xgs,wgs,ngauss_mad)
c
c     facl apears in expresion for plane wave in terms of bessel functions
      do l=0,lmax
        facl(l)=sqrtm1**l
      enddo
c     store all factors indepedent of contributing site
c     store bessel functions and legendre polynomials for all k except k=0
c     could move outside madsum
c     knlat is in a.u.; it is [vj X vk]/omega*2pi
c     ----------------------------------------------------------
      call zeroout(jl,(iplmax_mad+1)*ipknlat_mad)
c     ----------------------------------------------------------
      do i = 2,lastkn
c       length of k squared
        xknsq(i)=knlat_x(i)*knlat_x(i) +
     >           knlat_y(i)*knlat_y(i) +
     >           knlat_z(i)*knlat_z(i)
c       length of k
        xknabs(i)=sqrt(xknsq(i))
c       ----------------------------------------------------------
        call bessj(lmax+1,xknabs(i)*rmau,jl(1,i),work)
c       ----------------------------------------------------------
        xknabsi=one/xknabs(i)
        zk=knlat_z(i)*xknabsi
        xk=knlat_x(i)*xknabsi
        yk=knlat_y(i)*xknabsi
c       ----------------------------------------------------------
! meis: changed to normalized associated Legendre functions
        call plm_normalized(lmax,zk,plmk(1,i))
c       ----------------------------------------------------------
C       if(yk.eq.zero.and.xk.eq.zero)then
C         phi=zero
C       else
C         phi=atan2(yk,xk)
C       endif
C       expphi=exp(sqrtm1*phi)
        rk=sqrt(xk*xk+yk*yk)
	if(rk.eq.zero) then
	expphi=cone
	else
        expphi=dcmplx(xk/rk,yk/rk)
	endif
        tmp=cone
        expphim(0)=cone
        do m=1,lmax
          tmp=tmp*expphi
          expphim(m)=tmp
        enddo
        expfac=(exp( -xknsq(i)/twoetasq )/xknsq(i) )
        do j=1,ndlm
          l=lofj(j)
          m=mofj(j)
          cylm=conjg(clm(j)*expphim(m)*plmk(j,i))
          fack(i,j)=
     >    expfac*jl(l+1,i)*(four*pi)**2/omegbra*facl(l)*cylm
        enddo
      enddo
c
c    ------------------------------------------------------------------
c    start madmat calculation
c    ------------------------------------------------------------------
c     loop over contributing atom positions based on basis in distant cell
c     we are actually subtracting contributions from inside rcut which
c     are put in by k-space sum.
      do nsub2 =1,num_atoms
        call zeroout(dqintloc,2*ndlm)
        if(nsub2.eq.myatom)then
          ibegin=2
        else
          ibegin=1
        endif
c       vector from contributing site to my site.
        aijau(1) = atom_posi_x(myatom) - atom_posi_x(nsub2)
        aijau(2) = atom_posi_y(myatom) - atom_posi_y(nsub2)
        aijau(3) = atom_posi_z(myatom) - atom_posi_z(nsub2)
c    ------------------------------------------------------------------
c    subrtact rsn from aij to recalculate rssq which is used in 
c    calculating the real-space integral
c    ------------------------------------------------------------------
c           rdif is distance to center of gaussian to be removed
        do icell = 1,lastrs
c         vectors from mysite to contributing sites
          rtmp1 = (rslat_x(icell) - aijau(1) )
          rtmp2 = (rslat_y(icell) - aijau(2) )
          rtmp3 = (rslat_z(icell) - aijau(3) )
          rssq(icell) = rtmp1*rtmp1 + rtmp2*rtmp2 + rtmp3*rtmp3
c         if rtmp is inside rcut we need to subtract the potential of a
c         gaussian charge at that point.
          if(rssq(icell).lt.rcut2)then
            rdif=sqrt(rssq(icell))
c           phi is the angle to the center
C           if(rtmp2.eq.zero.and.rtmp1.eq.zero)then
C             phi=zero
C           else
C             phi=atan2(rtmp2,rtmp1)
C           endif
c           load exp(i*m*phi) for rotation matrix
C           expphi=exp(sqrtm1*phi)
        rk=sqrt(rtmp1*rtmp1+rtmp2*rtmp2)
	if(rk.eq.zero) then
	expphi=cone
	else
        expphi=dcmplx(rtmp1/rk,rtmp2/rk)
	endif
            tmp=cone
            expphim(0)=cone
            do m=1,lmax
              tmp=tmp*expphi
              expphim(m)=tmp
            enddo
c           find cos(theta) to center of gaussian
            if(rdif.ne.zero) then
              z=rtmp3/rdif
            elseif(rtmp3.eq.zero) then
              z=zero
	    else
	      stop'madewd_f'
            endif
c           find the integral of the potential times spherical harmonic
c           over a sphere of radius rm around mysite from a gaussian located
c           a distance rdif along the z-axiz.  If rdif is large compared
c           to 1/eta use interf. if rdif is small use interfsmr. only m=0
c           harmonics give non-zero values because of azmuthal symetry.
c           dqintl contains these integrals indexed by l=0,lmax
C           if(rdif*eta.gt.one)then
              call interf(plm,
     >        twoetasq,rmau,lmax,ndlm,mofj,lofj,clm,xgs,wgs,ngauss_mad,
     >        rdif,dqintl)
C           else
C             call interfsmr(plm,
C    >        twoetasq,rmau,lmax,ndlm,mofj,lofj,clm,
C    >        xgs,wgs,ngauss_mad,ngaussr_mad,rdif,dqintl)
C           endif
c           get legendre polynomial for cos(theta)
c       ----------------------------------------------------------
! meis: changed to normalized associated Legendre functions
            call plm_normalized(lmax,z,plm)
c       ----------------------------------------------------------
c           rotate dqintl by theta,phi using trans; which is particularly
c           simple because dqintl has only m=0 parts. trans is a spherical
c           harmonic times sqrt(fpi/(2l+1))
            do j=1,ndlm
              l=lofj(j)
              m=mofj(j)
              trans=
     >        clm(j)*plm(j)*expphim(m)*sqrt(four*pi/(2*l+1))
	      trans=conjg(trans)
              dqintloc(j)=dqintloc(j)+trans*dqintl(l)
            enddo
          endif
        enddo
c***************************************************************
c       calculate the difference in the potential of a gaussian
c       charge and a point charge at a distance squared of rssq
c       the potential at only the myatom lattice point is
c       calculated.  this is the l=0 part of the potential at rm.
c       the higher l's should be small and vanish as rcut tends to
c       infinity.  We are ignoring the higher l contributions.
c       call pqintg(twoetasq,rssq,lastrs,dqint,ibegin)
c       dqint is feed to madsum where the sum over distant atoms is made.
c       ------------------------------------------------------------------
c       perform the reciprocal-space sum and real-space sum
c       ----------------------------------------------------------------
c       madsum gives a matrix which generates the potential from charges outside
c       the central unit cell plus the potential from gaussian charges of width
c       determined by eta located in the central cell at basis positions of 
c       atoms in distant unit cells.
c      -----------------------------------------------------------------
c       calculate real space sum (sum2) and recip. space
c       sum (sum1).
c      -----------------------------------------------------------------
c       note sum starts at 2 since kn=0.0 of 1/(rs-aij) is
c       canceled by kn=0.0 of the 1/rs sum.
c      -----------------------------------------------------------------
c       do k space sum
        do j=1,ndlm
          sum1(j)=czero
          do i = 2,lastkn
c           k dot r(mysite)-r(i2)
            dot = 
     >      knlat_x(i)*aijau(1)+knlat_y(i)*aijau(2)+knlat_z(i)*aijau(3)
            expidot=exp( sqrtm1*dot )
            sum1(j) = sum1(j) + fack(i,j)*expidot
          enddo
        enddo
        sum2=czero
        do i=ibegin,lastrs
	   if(rssq(i).gt.rcut2)then
             rdif=sqrt(rssq(i))
             sum2 = sum2 
     >       +erfc(rdif*eta)/rdif
	   endif
	enddo
       if(iprint.ge.1) then
        write(6,'(''madewd_f::nsub2,sum2#'',i5,2d16.8)')nsub2,real(sum2)
       endif
c       sum1 is a potential expanded in spherical harmonics,  the l=0 term
c       must be multiplied by y0 inorder to yeild contribution to potential
c       the sum2 term is approximated by its l=0 term only, it must be
c       multiplied by 1/y0 to be consistent with the sum1 term.
        sum2=sqrt(four*pi)*sum2
c       add the k space contribution stored in sum1.
        do j=1,ndlm
          madmat(j,nsub2) = sum1(j)
        enddo
c       add the real space sum from sites outside rcut
c       recall that it is approximated by its l=0 part
        madmat(1,nsub2)=madmat(1,nsub2)+sum2
c       add the constant which sets the constant in the potential 
c       and makes it stationary in eta.
        term = four*pi/( omegbra*twoetasq )
        madmat(1,nsub2)=madmat(1,nsub2)-term*sqrt(four*pi)
c       subtract contributions from gaussian charges at perfect crystal
c       sites inside the celtral cell so that madmat willcorrespond to
c       contributions from atoms outside the central unit cell.
        do j=1,ndlm
          madmat(j,nsub2)=madmat(j,nsub2)+dqintloc(j)
        enddo
       if(iprint.ge.1) then
         write(6,'(''sum1,sum2,dqintloc='',6d16.8)')
     >   real(sum1(1)),real(sum2),real(dqintloc(1))
         write(6,'('' mad=='',5f12.4)')(madmat(j,nsub2),j=1,ndlm)
       endif
c    ----------------------------------------------------------------
      enddo                                  ! end do loop over nsub2
c    ----------------------------------------------------------------
      if( istop .eq. sname ) then
        call fstop(sname)
      endif
c    ----------------------------------------------------------------
      return
      end
