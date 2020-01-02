      subroutine rwave(srelflg,irregflg,rr,ri,drr,dri,matom,l,energy
     *     ,pot,nr,rad,potshift,jmt,kmax,emach,eunit,ws,wsdim,r_sph)
c
      implicit none
c
c---- dimensioning:
c
      integer kmax              ! largest grid point 
c
c---- mode flags:
c
      logical srelflg           ! scalar relativistic if true
      logical irregflg          ! find irregular solution if true
c
c---- output:
c
      complex*16 rr(kmax)      ! return regular solution
      complex*16 ri(kmax)      ! return irregular solution
      complex*16 drr(kmax)     ! return derivativ of regular solution
      complex*16 dri(kmax)     ! return derivativ of irregular solution
      complex*16 rrsph          !  rr at r_sph
      complex*16 drrsph          ! drr at r_sph
      complex*16 risph          !  rr at r_sph
      complex*16 drisph          ! drr at r_sph
      complex*16 matom          ! inverse of t(l) 
      complex*16 dummy          !  store derivatives
      real*8 work(2*kmax)          !  store derivatives
c
c---- input:
c
      integer l                 ! angular momentum number
      complex*16 energy         ! energy
      integer nr
      real*8 pot(nr,kmax)          ! potential
      real*8 rad(kmax)          ! radial grid
      complex*16 potshift       ! potential shift (may be complex)
      integer jmt             ! 
      real*8 emach              ! machine precision
      real*8 eunit             ! 1=Ryd, 2=Hartree
      real*8 r_sph             ! radius of spherical potential
c
! -- replace commonblock /dfbits/ in initwave, dfv:
      real*8 b,allp1

c
c---- work space:
c
      integer wsdim
      complex*16 ws(0:(wsdim/2-1),2)       ! workspace for bessel functions
c
c---- local variables:
c
      complex*16 einside,ka,kar,bj,dbj,cnorm,tatoml,dfi,fi,
     *     bh1,dbh1,eiz,cdum
      integer nv
      parameter (nv=4)
      real*8 y(nv),dy(nv),yp(nv),al,r,
     *     rfrom,rto,htry,hdid,hnext
      real*8 escale
      integer i,j,jm,jp,idid
      integer istart
c     ==================================================
c     dimensioned to stop run time array bounds problems
c     ==================================================
      integer isymm(1)
      real*8 vso(1,1)
c
      external dfv
c
c
      if(wsdim .lt. (max(l+1,2) * 2))then
         write(6,*) 'ERROR: wsdim=',wsdim,' < minimum =',
     *        (max(l+1,2) * 2)
         stop 'stoped in rwave (see unit 6)'
      endif
c     aplying a shift to the potential is like shifting the energy
c     inside the sphere 
      einside=energy+potshift
c     energy outisde sphere (used for bessel functions) is not shifted 
      ka=sqrt(eunit*energy+(0.0d0,1.0d0)*emach)
c
      escale=0.5d0*eunit
      al=dfloat(l)
      if(rad(1).eq.0.d0) then
      rr(1)=0.0d0
      drr(1)=0.0d0
      if(l.eq.0) drr(1)=1.0d0
      istart=2
      else
      istart=1
      endif
      call initwave(srelflg,rad,pot(1,1),nr,nv,istart,l,einside,
     &     y,escale,.false.,b,allp1)
c
      jp=istart
      call dfv(rad(istart),y,dy,nv,nr,rad(jp),einside,pot(1,jp),
     &     vso,1,isymm,eunit,b,allp1)
      rr(istart) =dcmplx(y(1),y(2))
      drr(istart)=dcmplx(dy(1),dy(2))
c
c----- integration proceeds out to grid point jmt, let routine decide
c      where to do intermediate calculations
c
      do j=istart+1,jmt
        yp(1)=y(1)
        yp(2)=y(2)
        yp(3)=y(3)
        yp(4)=y(4)
        jm=j-1
        jp=j-1
        rfrom=rad(jm)
        rto=rad(j)
        htry=rto-rfrom
        call bulirsch_stoer_integrator(rfrom, rto, y, dy, nv,
     &      rad(jp),einside,pot(1,jp),eunit,b,allp1,nv)
            rr(j) =dcmplx(y(1),y(2))
            drr(j)=dcmplx(dy(1),dy(2))
      enddo
c
c----- calculate the t-matrix and normalisation
c
c----- dfi is the logarithmic deriavtive at the sphere radius
c      fi is the wavefunction at the sphere radius
c
      dfi=drr(jmt)/rr(jmt)   !-1.0d0/rad(jmt)
      fi=rr(jmt)/rad(jmt)
      kar=ka*rad(jmt)
      call cinterp(rad,drr,jmt,r_sph,drrsph,dummy,work)
      call cinterp(rad,rr,jmt,r_sph,rrsph,dummy,work)
      dfi=drrsph/rrsph   !-1.0d0/rad(jmt)
      fi=rrsph/r_sph
      kar=ka*r_sph
c
c----- match phi to (jl + tl*hl)
c
c------ bj=ws(l,1), bh=ws(l,2)*eiz/kar, eiz=exp(i*kar)
c
      call zsphbesjh(max(l,1),kar,ws(0,1),ws(0,2),eiz,1)
c
c------ in subsequent code will define bj=jl, bh1=h1*r
c
      bj=ws(l,1)
      bh1=ws(l,2)*eiz/ka
c
c------ dbj=djl+jl/r, dbh1=r*dh1+h1, note that dfi above is really
c       drr/rr all is now consistent
c
      if(l.eq.0) then
        dbj=-ka*(ws(1,1)-ws(0,1)/kar)
        dbh1=-(ws(1,2)-ws(0,2)/kar)*eiz
      else
        cdum=al/kar
        dbj=ka*(ws(l-1,1)-ws(l,1)*cdum)
        dbh1=(ws(l-1,2)-ws(l,2)*cdum)*eiz
      endif
c
      cdum=-(bj*dfi-dbj)/(bh1*dfi-dbh1)
      cnorm=(bj+cdum*bh1)/fi
c     tatoml=rad(jmt)*cdum
      tatoml=r_sph*cdum
      matom=(0.d0,-1.0d0)*ka/tatoml
c
      do j=1,jmt
        rr(j) =cnorm*rr(j)
        drr(j)=cnorm*drr(j)
      enddo
c
c------ only find irregular solution if required to
c
      if(irregflg) then
c
        r=r_sph
c       r=rad(jmt)
        kar=ka*r
c bh1 is missing eiz/(ka*r)
      call zsphbesjh(max(l,1),kar,ws(0,1),ws(0,2),eiz,1)
      bh1=ws(l,2)*eiz/ka
      if(l.eq.0) then
        dbh1=-ws(1,2)*eiz
      else
        cdum=(l+1)/kar
        dbh1=(ws(l-1,2)-ws(l,2)*cdum)*eiz
      endif
      y(1)=dreal(bh1)
      y(2)=dimag(bh1)
      y(3)=dreal(dbh1)
      y(4)=dimag(dbh1)
      jp=jmt
c     call dfv(r,y,dy,nv,nr,rad(jp),einside,pot(1,jp),vso,1,isymm,
c    &     eunit)
      call dfv(r,y,dy,nv,nr,r_sph,einside,pot(1,jp),vso,1,isymm,
     &     eunit,b,allp1)
      ri(jmt) =dcmplx(y(1),y(2))
      dri(jmt)=dcmplx(dy(1),dy(2))
      if(rad(1).eq.0.d0) then
      ri(1) =0.0d0
      dri(1)=0.0d0
      endif
c
c----- start the inwards integration of irregular solution
c
      do j=jmt-1,istart,-1
        yp(1)=y(1)
        yp(2)=y(2)
        yp(3)=y(3)
        yp(4)=y(4)
        jm=j
        jp=j
        rfrom=rad(j+1)
        if(j+1.eq.jmt)rfrom=r_sph
        rto=rad(jm)
        htry=rto-rfrom
        call bulirsch_stoer_integrator(rfrom, rto, y, dy, nv,
     &      rad(jp),einside,pot(1,jp),eunit,b,allp1,nv)
            ri(j) =dcmplx(y(1),y(2))
            dri(j)=dcmplx(dy(1),dy(2))
      enddo
c
      endif !irregflg
c
      do i=jmt,kmax
c     do i=jmt+1,kmax
        kar=rad(i)*ka
        call zsphbesjh(max(l,1),kar,ws(0,1),ws(0,2),eiz,1)
        bj=ws(l,1)
        bh1=ws(l,2)*eiz/ka
        if(l.eq.0) then
          dbj=-ka*(ws(1,1)-ws(0,1)/kar)
          dbh1=-(ws(1,2)-ws(0,2)/kar)*eiz
        else
          cdum=al/kar
          dbj=ka*(ws(l-1,1)-ws(l,1)*cdum)
          dbh1=(ws(l-1,2)-ws(l,2)*cdum)*eiz
        endif
        rr(i)=bj*rad(i)+bh1*tatoml
        ri(i)=bh1
        drr(i)=dbj*rad(i)+dbh1*tatoml
        dri(i)=dbh1
      enddo
c
      return
      end
