c
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
c      single scatterer .......   self-consistent charge density
c
c         new radial r grid (doubling grid:  hermann-skillman style)
c
c                        written july,aug,dec,88 jan,feb,89   fjp
c =====================================================================
        implicit real*8 (a-h,o-z)
c
c     subroutine  .....    1  redin   .....    1
c     subroutine  .....    2  radgrd  .....    2
c     subroutine  .....    3  mtwsrd  .....    3
c     subroutine  .....    4  rdpot   .....    4
c     subroutine  .....    5  oldrdp  .....    5
c     subroutine  .....    6  atmrdp  .....    6
c     subroutine  .....    7  stdrdp  .....    7
c     subroutine  .....    8  stdwpt  .....    8
c     subroutine  .....    9  oldwpt  .....    9
c     subroutine  .....   10  genpte  .....   10
c     subroutine  .....   11  alpha2  .....   11
c     subroutine  .....   12  zeropt  .....   12
c     subroutine  .....   13  keepem  .....   13
c     subroutine  .....   14  ensum   .....   14
c     subroutine  .....   15  egrid   .....   15
c     subroutine  .....   16  solvit  .....   16
c     subroutine  .....   17  ricbes  .....   17
c     subroutine  .....   18  findns  .....   18
c     subroutine  .....   19  nbound  .....   19
c     subroutine  .....   20  bstate  .....   20
c     subroutine  .....   21  deepst  .....   21
c     subroutine  .....   22  shalst  .....   22
c     subroutine  .....   23  outsch  .....   23
c     subroutine  .....   24  inwsch  .....   24
c     subroutine  .....   25  qfit    .....   25
c     subroutine  .....   26  qrfit   .....   26
c     subroutine  .....   27  hank    .....   27
c     subroutine  .....   28  hcoef   .....   28
c     subroutine  .....   29  gamhnk  .....   29
c     subroutine  .....   30  bessj   .....   30
c     subroutine  .....   31  ptable  .....   31
c     subroutine  .....   32  intqud  .....   32
c     subroutine  .....   33  smpold  .....   33
c     subroutine  .....   34  intquc  .....   34
c     subroutine  .....   35  derv5   .....   35
c     subroutine  .....   36  derv4   .....   36
c     subroutine  .....   37  derv3   .....   37
c     subroutine  .....   38  ylag    .....   38
c     subroutine  .....   39  ylagc   .....   39
c     subroutine  .....   40  dgamix  .....   40
c     subroutine  .....   41  dpdif   .....   41
c     subroutine  .....   42  savit   .....   42
c     subroutine  .....   43  sumand  .....   43
c     subroutine  .....   44  mxlneq  .....   44
c     subroutine  .....   45  summag  .....   45
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
c
c /////////////////////////////////////////////////////////////////////
c      .......  parameters ........
c      ipcomp: number of components
c      ipspin: number of spins
c      iplmax: max. ang. mom.
c      ipmesh: number of points of r-grid
c      ipsize: size of block in radial grid
c      ipsts: number of core and semi-core states
c      ipits: number of iterations used in dganderson mixing method
c /////////////////////////////////////////////////////////////////////
c
        parameter (ipmesh=641,ipsize=40,iplmax=3)
        parameter(ipspin=2,ipcomp=4,ipsts=20,ipdim=ipmesh*ipcomp*ipspin)
        parameter (ipdim2=ipcomp*ipspin,ipits=4)
c
c /////////////////////////////////////////////////////////////////////
c      .......... names .........
c      ptflin:  potential file (input)
c      ptflot:  potential file (output)
c      keefil:  keep file (total energy info)
c      title:   title for run
c      comnam, comsym: component chemical name and symbol
c /////////////////////////////////////////////////////////////////////
c
        character*20 ptflin,ptflot,keefil
        character*78 title
        character*15 comnam(ipcomp)
        character*2  comsym(ipcomp)
c
c /////////////////////////////////////////////////////////////////////
c       ............   arrays   ..........................
c         r:  radial r grid (atomic units)
c         v:  potential actually r*v in rydbergs
c         vnew:  potential from new charge density (as above)
c         rho:  charge density   4*pi*r*r*rho
c         rhoc:  core charge density
c         rhosc:  semi-core (semi-valence) charge density
c         vchd:  valence charge density
c         rhomix: mixed charge density
c         rhoins: input charge densities
c         rhoots: output charge densites
c       ............   arrays   ..........................
c /////////////////////////////////////////////////////////////////////
c
        dimension     r(ipmesh)
        dimension     v(ipmesh,ipcomp,ipspin)
        dimension  vnew(ipmesh,ipcomp,ipspin)
        dimension   rho(ipmesh,ipcomp,ipspin)
        dimension  rhoc(ipmesh,ipcomp,ipspin)
        dimension rhosc(ipmesh,ipcomp,ipspin)
        dimension  vchd(ipmesh,ipcomp,ipspin)
        dimension rhomix(ipmesh+1,ipcomp,ipspin)
        dimension rhoins(ipits,ipmesh+1,ipcomp,ipspin)
        dimension rhoots(ipits,ipmesh+1,ipcomp,ipspin)
c
c /////////////////////////////////////////////////////////////////////
c       ............   arrays   ..........................
c       core   and  semi-core   information
c            numc and numscs:  number of bound states
c            npc  and npsc  :  princ. quantum numbers
c            lc   and lsc   :  ang. mom. numbers
c            ec   and esc   :  eigen-energies
c            occ  and ocsc  :  occupation (2*l+1)
c /////////////////////////////////////////////////////////////////////
c
        dimension numc(ipcomp,ipspin),numsc(ipcomp,ipspin)
        dimension    npc(ipsts,ipcomp,ipspin)
        dimension     lc(ipsts,ipcomp,ipspin)
        dimension     ec(ipsts,ipcomp,ipspin)
        dimension   npsc(ipsts,ipcomp,ipspin)
        dimension    lsc(ipsts,ipcomp,ipspin)
        dimension    esc(ipsts,ipcomp,ipspin)
        dimension    occ(ipsts,ipcomp,ipspin)
        dimension   ocsc(ipsts,ipcomp,ipspin)
        dimension qsum(ipcomp,ipspin),qsums(ipcomp,ipspin)
c
c /////////////////////////////////////////////////////////////////////
c       ............   arrays   ..........................
c        zed and z1:  atomic number  of each component
c        conc: concentration  of each component
c        zv1:   total valence per component (summed over spins)
c        zv2:    valence per componet per spin
c        xmom1:  moment for each component
c                 only average moment will be fixed
c        zcore:  core charge per component per spin
c        zsemi:  semi-core (semi-valence) charge per comp. per spin
c        ef1 and efest:  fermi energy for each spin
c                        efest is current estimate,
c                        ef1 is only the initial guess
c        ztot and zto1:  charge in cell of each component per spin
c                        ztot is used in fixed moment scheme
c                        zto1 is used with the currect charge density
c                        ztot  corresponds  to correct moment
c                        zto1  corresponds  to current moment
c        vws:    valence and semi-core charge in cell per spin
c        vws1:  same as above but corresponds to old charge density
c        zc1:  core charge (total)
c        zsc1: semi-core charge (total)
c        zval: valence charge (total)
c        qmt:  total charge within muffin-tin radius
c        xmmt: moment within the muffin-tin radius
c
c        esemi:    one-electrom sum (over semi-core states)
c        evsum: one-electron energies summed over semi-core and valence
c        rhsrho:  root-mean-squared difference in the charge density
c
c /////////////////////////////////////////////////////////////////////
c
c
        dimension zed(ipcomp),conc(ipcomp),xmom1(ipcomp),zval(ipcomp)
        dimension qmt(ipcomp),zcore(ipcomp,ipspin),zsemi(ipcomp,ipspin)
        dimension xmmt(ipcomp),zv2(ipcomp,ipspin),rmsrho(ipcomp,ipspin)
        dimension z1(ipcomp),zc1(ipcomp),zsc1(ipcomp),zv1(ipcomp)
        dimension vws(ipcomp,ipspin),vws1(ipcomp,ipspin)
        dimension ztot(ipcomp,ipspin),zto1(ipcomp,ipspin)
        dimension esemi(ipcomp,ipspin),evsum(ipcomp,ipspin),ef1(ipspin)
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
c
        data rho/ipdim*0.0/,rhoc/ipdim*0.0/,rhosc/ipdim*0.0/
        data v/ipdim*0.0/,vnew/ipdim*0.0/,vchd/ipdim*0.0/
        data evsum/ipdim2*0.0/,thresh/1.e-10/,etopc/-10.0/
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
c
c   mesh is the number of points on the r-grid
c   the r-grid is broken into blocks (nsize in length)
c
           mesh=ipmesh
           nsize=ipsize
c
c .....................................................................
c
c  read the pertinent information for this run
c    ipr is the print level
c   alat is the lattice constant
c   lmax is the largest ang. momemtum number
c  icrys indicates fcc (icry=1), bcc (icrys=2), or sc (icrys=3)
c  ncomp is the number of components ( max. set by ipcomp )
c nspins is the number of spins ( 1 or 2 )
c  z1 are the atomic numbers
c  conc store the concentrations
c  iexc determine which exchange-correlation potential
c  xmom1 stores the moments (for fixed moment calculations)
c  rmt and rws are the muffin-tin and wigner sietz radii
c  ipot and jpot control potential input and output
c  iunin, iunot  are units for potential input and output
c  iunkp is the unit for the keep file (record energy and press)
c  ndiv,ebott,eimag,ef1  are parameters for energy grid
c  nscfit is the number of scf iterations
c  alpa is the mixing parameter
c
c .....................................................................
         call redin(ipr,alat,lmax,icrys,ncomp,conc,nspins,iexc,z1,
     >      xmom1,rmt,rws,ipot,jpot,ndiv,ebott,eimag,ef1,nscfit,alpa,
     >             iunin,iunot,iunkp)
c .....................................................................
c
c  set up various parameters: rmt, rws, oma, cay
c     rmt and rws are the muffin-tin and wigner-sietz radii
c     oma*alat**3 is the interstial volume
c     cay is the constant for the madelung term in the energy
c
c .....................................................................
            call mtwsrd(icrys,alat,rmt,rws,oma,cay,ipr)
c .....................................................................
c  set up the r-grid:  doubling grid
c    jmt is the index corresponding to the muffin-tin radius
c    jws is the index corresponding to the grid point is
c          closet to the w-s radius (but larger)
c
c .....................................................................
            call radgrd(r,mesh,nsize,mblock,rmt,jmt,rws,jws,ipr)
            nblock=mesh/nsize
c
c .....................................................................
c   read in    the potentials: v
c        the charge densities: rho
c     the core charge density: rhoc
c     the core and semi-core bound state info:
c         numc,npc,lc,ec,occ    and    numsc,mpsc,lsc,esc,ocsc
c
c .....................................................................
             do 6 ispin=1,nspins
               is=ispin
               do 5 icomp=1,ncomp
               ic=icomp
c .....................................................................
               call rdpot(ipot,iunin,r,mesh,nsize,mblock,jmt,jws,ic,is,
     >                    v(1,ic,is),v0,rho(1,ic,is),rhoc(1,ic,is),
     >                    alat,zed(ic),zcore(ic,is),vws(ic,is),efx,
     >                     numc,npc,lc,ec,occ,
     >                     numsc,npsc,lsc,esc,ocsc,ipr)
c .....................................................................
                zto1(ic,is)=zcore(ic,is)+vws(ic,is)
                ztot(ic,is)=.5*z1(ic)+(nspins-1)*(1.5-is)*xmom1(ic)
                if(ipr.gt.1) then
                     write(6,*)' ----------- total charge: site=',is,
     >              '  spin=',is,'    z=',ztot(ic,is),'  -------------'
                endif
c .....................................................................
c  get chemical symbol and name for periodic table
c
c .....................................................................
               iz=zed(ic)
               call ptable(iz,comsym(ic),comnam(ic))
c .....................................................................
c
c     determine difference in potential zeroes:  vdif
c
c +++++++++++++++++++++++++++++  vdif= v0(down)-v0(up)  +++++++++++++++
               v0up=(2-is)*v0
               v0dn=(is-1)*v0
  5        continue
           if(ef1(is).lt.0.001) ef1(is)=efx
  6        continue
c
c
c ---------------------------------------------------------------------
             close(iunin)
             vdif=(nspins-1)*(v0dn-v0up)
             if(nspins.gt.1)
     >         write(6,*) ' muffin-tin-zero difference =  ',vdif
c .....................................................................
       if(ipot.lt.0) then
c                                 starting with atomic ch. density
c                                        regenerate the potentials
c _____________________________________________________________________
c
            do 2 ic=1,ncomp
            zsemi(ic,nspins)=0
            zsemi(ic,  1) =0
            zv2(ic,nspins)=0
            zv2(ic,1)=0
            totalz=zto1(ic,1)+zto1(ic,nspins)
            ztot(ic,1)=.5*totalz+.5*(nspins-1)*xmom1(ic)
            ztot(ic,nspins)=totalz-ztot(ic,1)
    2       continue
c _____________________________________________________________________
       itoten=2
      call genpte(vnew,vdfnw,v,vdif,ncomp,nspins,mesh,nsize,nblock,
     >  jmt,jws,rmt,rws,r,rho,rhoc,rmsrho,
     >  conc,zed,vws,zcore,qmt,xmmt,
     >  cay,oma,icrys,alat,iexc ,qsum,evsum,itoten,ipr)
c _____________________________________________________________________
             vdif=vdfnw
c
            do is=1,nspins
            do ic=1,ncomp
            do i=1,mesh
            v(i,ic,is)=vnew(i,ic,is)
            enddo
            enddo
            enddo
c
        endif
c _____________________________________________________________________
c
c                 ------  beginning of  scf loop  --------
c                 ------  beginning of  scf loop  --------
c                 ------  beginning of  scf loop  --------
c                 ------  beginning of  scf loop  --------
c                 ------  beginning of  scf loop  --------
c                 ------  beginning of  scf loop  --------
c                 ------  beginning of  scf loop  --------
c                 ------  beginning of  scf loop  --------
c
c _____________________________________________________________________
c
      do 500 it=1,nscfit
      iter=it
      itermx=it
c .....................................................................
c
c   store the  input charge density in rhoins
c .....................................................................
       call savit(itermx,rho,jws,zto1,ncomp,nspins,rhoins)
c .....................................................................
c
       do is=1,nspins
       do ic=1,ncomp
          qsum(ic,is)  = 0.0
          qsums(ic,is) = 0.0
       enddo
       enddo
c .....................................................................
c
c  find all bound states:  both core and semi-core
c  all states found below etopc (-10 ryd) are treated as core
c  all other states below ebott are combined with the valence
c      ebott is bottom of contour used in integrating
c      over the valence states
c
      etopsc = ebott
      do ispin=1,nspins
      is=ispin
      do icomp=1,ncomp
      ic=icomp
c /////////////////////////////////////////////////////////////////////
c .....................................................................
c             determine the number of core states
c .....................................................................
         call findns(r,mesh,jmt,jws,rws,ic,is,v(1,ic,is),zed(ic),
     >        numc,npc,lc,ec,etopc,etopc,thresh,ipr)
c .....................................................................
         if(numc(ic,is).gt.0) then
           if(ipr.gt.0) then
              write(6,'(/,1x,78(''_''))')
              write(6,
     &    '(/,18x,''core states:  component #'',i2,6x,''spin #'',i2,/)')
     &            ic,is
               endif
c .....................................................................
c     find the core states
c .....................................................................
               call bstate(r,mesh,jmt,jws,rws,ic,is,v(1,ic,is),
     >                    rhoc(1,ic,is),zed(ic),numc,npc,lc,ec,occ,
     >                    etopc,thresh,qsum(ic,is),ipr)
c .....................................................................
               if(ipr.gt.0) then
                  write(6,'(/,4x,
     &   ''n'',4x,''l'',2x,''energy(ryd)'',3x,''occ''/,3x,27(''-'') )')
                  nmcs=numc(ic,is)
                  write(6,'(2i5,f13.4,f7.2 )')
     >              (npc(i,ic,is),lc(i,ic,is),ec(i,ic,is),
     >                 occ(i,ic,is),i=1,nmcs)
               endif
c
          else
c
              if(ipr.gt.0) then
                  write(6,'(1x,78(''/''))')
                  write(6,*)' core charge density is set to zero'
                  write(6,'(1x,78(''/''))')
              endif
              zcore(ic,is)=0.
              do i=1,mesh
              rhoc(i,ic,is)=0.0
              enddo
c
          endif
c /////////////////////////////////////////////////////////////////////
c .....................................................................
c    determine the number of semi-core states
c .....................................................................
          call findns(r,mesh,jmt,jws,rws,ic,is,v(1,ic,is),zed(ic),
     >        numsc,npsc,lsc,esc,etopsc,etopc,thresh,ipr)
c .....................................................................
          if(numsc(ic,is).gt.0) then
             if(ipr.gt.0) then
                write(6,'(/,1x,78(''_''))')
                write(6,'(/,13x,a,i2,6x,a,i2,/)')
     >            'semi-core states:  component #',ic,'spin #',is
             endif
c .....................................................................
             call bstate(r,mesh,jmt,jws,rws,ic,is,v(1,ic,is),
     >                   rhosc(1,ic,is),zed(ic),numsc,npsc,lsc,esc,ocsc,
     >                   etopsc,thresh,qsums(ic,is),ipr)
c .....................................................................
             if(ipr.gt.0) then
                  write(6,'(/,4x,''n'',4x,''l'',2x,''energy(ryd)'',3x,
     >            ''occ''/,3x,27(''-'') )')
                  nmcs=numsc(ic,is)
                  write(6,'(2i5,f13.4,f7.2)')
     >              (npsc(i,ic,is),lsc(i,ic,is),esc(i,ic,is),
     >              ocsc(i,ic,is),i=1,nmcs)
             endif
c .....................................................................
         else
c
             if(ipr.gt.0) then
                   write(6,'(1x,78(''/''))')
                   write(6,*)' semi-core charge density is set to zero'
                   write(6,'(1x,78(''/''))')
             endif
             zsemi(ic,is)=0.
             do i=1,mesh
             rhosc(i,ic,is)=0
             enddo
c
         endif
c /////////////////////////////////////////////////////////////////////
c
c .....................................................................
c    found all bound states now determine number of valence
c       electrons for charge neutrality
c .....................................................................
c
           sum=0.0
           do i=1,numc(ic,is)
             sum=sum+occ(i,ic,is)
           enddo
           zcore(ic,is)=sum
           sum=0.0
           esum=0.0
           do i=1,numsc(ic,is)
ccc           ebott = min(ebott,esc(i,ic,is))
              esum=esum+ocsc(i,ic,is)*esc(i,ic,is)
              sum=sum+ocsc(i,ic,is)
           enddo
           esemi(ic,is)=esum
           zsemi(ic,is)=sum
           zv2(ic,is)=ztot(ic,is)-zcore(ic,is)-zsemi(ic,is)
c
       enddo ! ispin
       enddo ! icomp
c .....................................................................
                if(ndiv.eq.0) stop
c .....................................................................
c
c      calculate the valence charge density:  vchd
c      determine the fermi energy and one-electron energy sum
c
c .....................................................................
cc       ebott = max(-0.5d0,ebott)
         call ensum(vchd,evsum,ef1,zv2,vws,conc,
     >     ncomp,nspins,lmax,r,v,mesh,nsize,jmt,jws,rws,vdif,
     >     ndiv,ebott,eimag,ipr)
c .....................................................................
       write(6,'(1x,78('':''))')
c .....................................................................
       do 35 ic=1,ncomp
         do 30 is=1,nspins
              evsum(ic,is)=evsum(ic,is)+esemi(ic,is)
              vws(ic,is)=vws(ic,is)+zsemi(ic,is)
              ztot(ic,is)=zcore(ic,is)+vws(ic,is)
              sumdrh=0
              do 25 i=1,jws
                 rhtmp=rho(i,ic,is)
                 rho(i,ic,is)=rhoc(i,ic,is)+rhosc(i,ic,is)+vchd(i,ic,is)
                 sumdrh=sumdrh+(rho(i,ic,is)-rhtmp)**2
   25         continue
              rmsrho(ic,is)=sqrt( sumdrh/jws )
   30     continue
c _____________________________________________________________________
          zval(ic)=zed(ic)-(zcore(ic,1)+zcore(ic,nspins))
     >                   -(zsemi(ic,1)+zsemi(ic,nspins))
c _____________________________________________________________________
          write(6,'(3x,''component #'',i2,2x,''z='',
     >       f4.0,3x,a2,4x,''concentration='',f8.4,3x,a15 /
     >       3x,''atomic number='',f4.0,5x,''core: '',f4.0,
     >       3x,''semi-core: '',f4.0,3x,''valence: '',f4.0)')
     >       ic,zed(ic),comsym(ic),conc(ic),comnam(ic),
     >       zed(ic),zcore(ic,1)+zcore(ic,nspins),
     >       zsemi(ic,1)+zsemi(ic,nspins),zval(ic)
          write(6,'(1x,78('':''))')
 35   continue
c
c .....................................................................
c
c        store the output charge density in rhoots
c .....................................................................
       call savit(itermx,rho,jws,ztot,ncomp,nspins,rhoots)
c .....................................................................
c                                        regenerate the potentials
c                                       and calculate total energy
c _____________________________________________________________________
c
      itoten=1
      call genpte(vnew,vdfnw,v,vdif,ncomp,nspins,mesh,nsize,nblock,
     >  jmt,jws,rmt,rws,r,rho,rhoc,rmsrho,
     >  conc,zed,vws,zcore,qmt,xmmt,
     >  cay,oma,icrys,alat,iexc ,qsum,evsum,itoten,ipr)
c _____________________________________________________________________
c
c .....................................................................
c
c    use the d.g.anderson scheme to find optimize mixing
c     currently uses 4 previous iterations
c
c .....................................................................
       call dgamix(itermx,jws,ncomp,nspins,rhoins,rhoots,rhomix,alpa)
c .....................................................................
c
        do is=1,nspins
        do ic=1,ncomp
           zto1(ic,is)=rhomix(jws+1,ic,is)
           vws1(ic,is)=zto1(ic,is)-zcore(ic,is)
           if(ipr.gt.0) write(6,*)' ic=',ic,'  is=',is,'   after mix:'
     >            ,'  zto1=',zto1(ic,is),'   vws1=',vws1(ic,is)
           do i=1,jws
               rho(i,ic,is)=rhomix(i,ic,is)
           enddo
        enddo
        enddo
c
c      calculate new potential from mixed charge density
c
c _____________________________________________________________________
c
      itoten=2
      call genpte(vnew,vdfnw,v,vdif,ncomp,nspins,mesh,nsize,nblock,
     >  jmt,jws,rmt,rws,r,rho,rhoc,rmsrho,
     >  conc,zed,vws1,zcore,qmt,xmmt,
     >  cay,oma,icrys,alat,iexc ,qsum,evsum,itoten,ipr)
c _____________________________________________________________________
c
      vdif=vdfnw
      do ispin=1,nspins
         is=ispin
         vd=(is-1)*vdfnw
         do icomp=1,ncomp
            ic=icomp
            chmtcs=.5*qmt(ic)+(nspins-1)*(1.5-is)*xmmt(ic)
            do 88 i=1,jmt
               v(i,ic,is)=vnew(i,ic,is)
              if ( it .eq. nscfit ) then
                rhoc(i,ic,is)=rhoc(i,ic,is)+rhosc(i,ic,is)
              endif
 88         continue
c
            if ( it .eq. nscfit ) then
              do 89 i=jmt+1,jws
                rhoc(i,ic,is)=rhoc(i,ic,is)+rhosc(i,ic,is)
 89           continue
            endif
c
c .....................................................................
            if(jpot.le.1) then
               if(ipr.gt.1) then
                 call stdwpt(   6 ,ic,ncomp,    conc,comsym,comnam,
     >           zed(ic),vws(ic,is),   alat,ef1(is),mesh,nsize,mblock,
     >           jmt,jws,r,vnew(1,ic,is),vd,rho(1,ic,is),rhoc(1,ic,is)
     >           ,numc(ic,is), npc(1,ic,is), lc(1,ic,is), ec(1,ic,is),
     >           numsc(ic,is),npsc(1,ic,is),lsc(1,ic,is),esc(1,ic,is)
     >            ,chmtcs,is,nspins)
               endif
c
c .....................................................................
               call stdwpt(iunot,ic,ncomp,    conc,comsym,comnam,
     >         zed(ic),vws(ic,is),   alat,ef1(is),mesh,nsize,mblock,
     >         jmt,jws,r,vnew(1,ic,is),vd,rho(1,ic,is),rhoc(1,ic,is),
     >         numc(ic,is), npc(1,ic,is), lc(1,ic,is), ec(1,ic,is),
     >         numsc(ic,is),npsc(1,ic,is),lsc(1,ic,is),esc(1,ic,is)
     >           ,chmtcs,is,nspins)
c .....................................................................
            else
c ---------------------------------------------------------------------
               if(ipr.gt.1) then
                  call oldwpt(   6 ,ic,ncomp,    conc,comsym,comnam,
     >            zed(ic),vws(ic,is),  alat,ef1(is),mesh,nsize,mblock,
     >            jmt,jws,rws,r,vnew(1,ic,is),vd,
     >            rho(1,ic,is),rhoc(1,ic,is),2*zcore(ic,is),
     >            numc(ic,is),npc(1,ic,is),lc(1,ic,is),ec(1,ic,is),
     >            numsc(ic,is),npsc(1,ic,is),lsc(1,ic,is),esc(1,ic,is)
     >             ,chmtcs,is,nspins)
                endif
c
c .....................................................................
                call oldwpt(iunot,ic,ncomp,    conc,comsym,comnam,
     >          zed(ic),vws(ic,is),  alat,ef1(is),mesh,nsize,mblock,
     >          jmt,jws,rws,r,vnew(1,ic,is),vd,
     >          rho(1,ic,is),rhoc(1,ic,is),2*zcore(ic,is),
     >          numc(ic,is),npc(1,ic,is),lc(1,ic,is),ec(1,ic,is),
     >          numsc(ic,is),npsc(1,ic,is),lsc(1,ic,is),esc(1,ic,is)
     >            ,chmtcs,is,nspins)
c .....................................................................
c
            endif
         enddo
      enddo
c
        rewind(iunot)
c
c .............................   end of scf loop  ..................
c .............................   end of scf loop  ..................
c .............................   end of scf loop  ..................
c .............................   end of scf loop  ..................
c .............................   end of scf loop  ..................
c .............................   end of scf loop  ..................
c .............................   end of scf loop  ..................
c .............................   end of scf loop  ..................
c
c _____________________________________________________________________
  500    continue
c
      stop
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine redin(ipr,alat,lmax,icrys,ncomp,conc,nspins,iexc,z1,
     >       xmom1,rmt,rws,ipot,jpot,ndiv,ebott,eimag,ef1,nscfit,alpa,
     >             iunin,iunot,iunkp)
c =====================================================================
        implicit real*8 (a-h,o-z)
        character*20 infile,otfile,blank,ptflin,ptflot,keefil
        character*78 title
c
        parameter (ipcomp=4,ipspin=2)
c
        dimension conc(ipcomp),ef1(ipspin)
        dimension z1(ipcomp),xmom1(ipcomp)
c
        data blank/'                    '/
c
        isave=6
c
c ---------------------------------------------------------------------
c                     ntract=1   interactive mode for input
c                           =0   batch mode
c ---------------------------------------------------------------------
c
c
        write(6,*)' single scatterer program,  scf-total energy'
ccc     write(6,*)' do you wish to run this program in an interactive'
ccc  >       ,' mode? (1=yes,0=no)'
ccc     read(5,*) ntract
        read(5,'(a)') title
        ntract = 0
c
c ---------------------------------------------------------------------
        if(ntract.eq.0) then
c
ccc         write(6,*)' input for this program will be read from a file'
ccc         write(6,*)' please type the name of this file. (format:a20)'
ccc         read(5,'(a20)') infile
c
       else
c
            write(6,*)' program will read input interactively'
            write(6,*)' for future use the input information will be'
     >               ,' recorded in a file'
            write(6,*)' please type the name of this file. (format:a20)'
            read(5,'(a20)') infile
            if(infile.ne.blank) then
                 isave=10
                 open(unit=10,file=infile,form='formatted',status='old')
            endif
         endif
c ---------------------------------------------------------------------
c
ccc      write(6,*) 'please type name of output file. (format:a20)'
ccc      write(6,*) 'if blank: the program will continue typing to the'
ccc  >             ,' screen'
         read(5,'(a20)') otfile
c
c ---------------------------------------------------------------------
        if(ntract.eq.0) then
ccc        open(unit=5,file=infile,form='formatted',status='old')
           if(otfile.ne.blank) then
              open(unit=6,file=otfile,form='formatted',status='unknown')
           endif
            read(5,*) ipr
            read(5,*) icrys
            read(5,*) lmax
            read(5,*) alat
            read(5,*) nspins
            read(5,*) iexc
            read(5,*) ncomp
            do 1 ic=1,ncomp
            read(5,*) z1(ic),xmom1(ic)
  1         continue
            read(5,*) (conc(i),i=1,ncomp)
            read(5,*) rmt,rws
ccc         read(5,*) ipot,jpot
            ipot = -nspins; jpot = 2
            read(5,'(a20)') ptflin
            read(5,'(a20)') ptflot
            read(5,'(a20)') keefil
            read(5,*) ndiv,ebott,eimag
            read(5,*) (ef1(is),is=1,nspins)
            read(5,*) nscfit,alpa
            close(5)
c ---------------------------------------------------------------------
        else
c ---------------------------------------------------------------------
            write(6,*) ' input: title of run.  (format:a78)'
            read(5,'(a78)') title
            write(6,*)' input: print level (0=min,1=intermediate,2=max)'
            read(5,*) ipr
            write(6,*)' input: crystal type (1=fcc,2=bcc)'
            read(5,*) icrys
            write(6,*)' input: lattice constant (a.u.)'
            read(5,*) alat
            write(6,*)' input: max. ang.mom. quant. numb (lmax=0,1,2,3)'
            read(5,*) lmax
            write(6,*)' input: number of spins (nspins=1,2)'
            read(5,*) nspins
            write(6,*)' input: ex-corr. potential (1=vb-h,2=vosko)'
            read(5,*) iexc
            write(6,*)' input: number of components (ncomp=1,2,3,4)'
            read(5,*) ncomp
c
            csum=0.
            do 2 ic=1,ncomp
               write(6,*) ' ===========> component #',ic
               if(ncomp.eq.1) then
                   conc(1)=1.
                   write(6,*)' single component: concentration = 1.0'
               else
                   write(6,*) ' input: concentration'
                   read(5,*) conc(ic)
               endif
               write(6,*)' input: atomic number (zed)'
               read(5,*) z1(ic)
               write(6,*) ' input: moment desired'
               read(5,*) xmom1(ic)
               csum=csum+conc(ic)
  2         continue
c
  6         continue
            if(abs(1.-csum).gt.0.00001) then
                 write(6,*)' discrepancy in concentrations'
                 do 8 ic=1,ncomp
                   write(6,*) ' ===========> component #',ic
                   write(6,*) ' z=',z1(ic),'   moment=',xmom1(ic)
                   write(6,*) ' input: concentration'
                   read(5,*) conc(ic)
   8             continue
                 go to 6
            endif
c
            write(6,*)' input: mt radius and ws radius'
            write(6,*)'    (if zero, they will take'
     >          ,' on values corresponding to touching m.tins'
            read(5,*) rmt,rws
ccc         write(6,*)' input: type of input potential '
ccc         write(6,*)'        (0=new std,1=old para,2=old ferro)'
ccc         write(6,*)'         (-1=atomic para,-2=atomic ferro)'
ccc         read(5,*) ipot
ccc         write(6,*)' input: type of output potential',
ccc  >                       ' (1=new std,2=old)'
ccc         read(5,*) jpot
            ipot = -nspins; jpot = 2
            write(6,*)' input: beginning potential'
            write(6,*)' please type the name of this file. (format:a20)'
            read(5,'(a20)') ptflin
            write(6,*)' input:     output potential'
            write(6,*)' please type the name of this file. (format:a20)'
            read(5,'(a20)') ptflot
            write(6,*)' input: keep file (summary file containing '
     >                ,' total energy and rms info'
            write(6,*)' please type the name of this file. (format:a20)'
            read(5,'(a20)') keefil
            write(6,*)' input: energy grid information'
            write(6,*)' input: number of divisions per .1 rydbergs'
            write(6,*)'        positive for single site'
            write(6,*)'        negative for phase shifts'
            write(6,*)'        zero:  stops after finding bound states'
            read(5,*) ndiv
            write(6,*)' input: bottom of contour'
            read(5,*) ebott
            write(6,*)' input: shift off real axis'
            read(5,*) eimag
            do 9 is=1,nspins
                write(6,*) ' input: estimate of fermi energy'
     >                   ,' for spin ',is
                write(6,*) '        (if zero: will be adjusted to'
     >                   ,' correspond to value in potential file)'
                write(6,*) '         (if reading atomic potentials'
     >                  ,' should not be zero, use .5)'
                read(5,*) ef1(is)
   9        continue
            write(6,*) ' input: number of scf iterations'
            read(5,*) nscfit
            write(6,*) ' input: mixing parameter (typically .1)'
            read(5,*) alpa
        endif
        if(ndiv.le.0) ipr=max(1,ipr)
c .....................................................................
         write(6,'(a78)') title
         write(6,*) 'ipr=',ipr
         write(6,*) 'icrys=', icrys
         write(6,*) 'lmax=', lmax
         write(6,*) 'alat=',alat
         write(6,*) 'nspins=',nspins
         write(6,*) 'iexc=',iexc
         write(6,*) 'ncomp=',ncomp
         write(6,*) (' i=',i,' conc(i)=',conc(i),i=1,ncomp)
         write(6,*) 'ipot=',ipot
         write(6,*) 'jpot=',jpot
         write(6,*)' contour info:  bottom=',ebott,'    shift=',eimag,
     >              '   top=',(ef1(is),is=1,nspins),'    ndiv=',ndiv
         write(6,*) ' number of scf iterations is ',nscfit,
     >            '     with mixing parameter alpha=',alpa
         write(6,'(1x,'' input potential file is '',a20)') ptflin
         write(6,'(1x,''output potential file is '',a20)') ptflot
         write(6,'(1x,''            keep file is '',a20)') keefil
c .....................................................................
c
          write(isave,'(a78)') title
          write(isave,'(i5,t26,''print level'')') ipr
          write(isave,'(i5,t26,''crystal type (1=fcc,2=bcc)'')') icrys
          write(isave,'(i5,t26,''lmax'')') lmax
          write(isave,'(f10.5,15x,''lattice constant'')') alat
          write(isave,'(i5,t26,''number of spins'')') nspins
          write(isave,'(i5,t26,''ex-cor type (1=vb-hedin,2=vosko)'')')
     >                             iexc
          write(isave,'(i5,t26,''number of components'')') ncomp
          do 7 ic=1,ncomp
              write(isave,'(2(1x,f9.5),t26,''z(at.num),  mom'')')
     >         z1(ic),xmom1(ic)
  7       continue
          if(ncomp.eq.1) then
               write(isave,'(1x,f9.5,t26,''concentrations'')')
     >                  (conc(ic),ic=1,ncomp)
               go to 12
          endif
          if(ncomp.eq.2) then
               write(isave,'(2(1x,f9.5),t26,''concentrations'')')
     >                  (conc(ic),ic=1,ncomp)
               go to 12
          endif
          if(ncomp.eq.3) then
               write(isave,'(3(1x,f9.5),5x,''concentrations'')')
     >                  (conc(ic),ic=1,ncomp)
               go to 12
          endif
          if(ncomp.eq.4) then
               write(isave,'(4(1x,f9.5),5x,''concentrations'')')
     >                  (conc(ic),ic=1,ncomp)
               go to 12
          endif
          write(6,*) 'stopping, ncomp.ne.1,2,3 or 4'
          stop
c ----------------------------------------------------------------------
  12      continue
c ----------------------------------------------------------------------
            write(isave,'(1x,f9.5,1x,f9.5,5x,''mt radius, ws radius'')')
     >                  rmt,rws
            write(isave,'(2i5,t26,''pot type: in, out'')') ipot,jpot
            write(isave,'(a20,5x,''pot file  (input)'')') ptflin
            write(isave,'(a20,5x,''pot file (output)'')') ptflot
            write(isave,'(a20,5x,''keep file'')') keefil
            write(isave,'(i5,1x,f7.4,1x,f7.4,t26,
     >               ''e-grid: ndiv, bott, shift'')')
     >                 ndiv,ebott,eimag
            if(nspins.eq.1) then
                      write(isave,'(1x,f9.5,t26,''ef (est)'')')
     >                  (ef1(is),is=1,nspins)
            else
                      write(isave,'(2(1x,f9.5),t26,''ef (est)'')')
     >                  (ef1(is),is=1,nspins)
            endif
            write(isave,'(i5,1x,f9.5,t26,
     >             ''num of scf its, mix parameter'')')  nscfit,alpa
           if(isave.ne.6) close(isave)
c .....................................................................
           if(ntract.eq.1) then
                write(6,*)' do you wish to continue?  (1=yes,0=no)'
                read(5,*) icon1
                if(icon1.eq.0) stop
           endif
c .....................................................................
           iunin=33
           iunot=34
           iunkp=20
c
         open(unit=iunin,file=ptflin,form='formatted',status='old')
         open(unit=iunot,file=ptflot,form='formatted',status='unknown')
         open(unit=iunkp,file=keefil,form='formatted',status='unknown')
c        write(6,*) iunin,iunot,iunkp
c .....................................................................
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        subroutine radgrd(r,mesh,nsize,mblock,rmt,jmt,rws,jws,ipr)
c =====================================================================
        implicit real*8 (a-h,o-z)
        parameter (ipmesh=641)
        dimension r(ipmesh)
c _____________________________________________________________________
c
c     routine set up the r-grid, which is a
c     "doubling grid" similar to hermann-skillman type
c     each block contains nsize (typically 40) points
c     grid contains mesh (641) points, ipmesh/nsize=nblock blocks (16)
c     first grid point is zero, the second is the first increment
c     the increment doubles each block but never exceeds .01 au
c     beyond the mblock block the increment is fixed
c     the muffin-tin radius is set to be a the end of a block
c     the grid extends beyond the wigner-seitz radius by a few points
c _____________________________________________________________________
c
        nblock=mesh/nsize
        drmax=.01
        dr=rmt/nsize
c ---------------------------------------------------------------------
        mblock=nblock
c
        if(dr.gt.drmax) then
            xidr=dr/drmax
            idr = xidr + 2
            mblock=nblock - idr
c
            if(ipr.gt.2) then
                write(6,'(1x,78(''=''))')
                write(6,'('' adjusting the doubling'',
     >                     '' in the radial grid'')')
                write(6,'('' will double increment only for the first'',
     >                     i3,'' blocks'')' ) mblock
                write(6,'('' keep increment constant for all'',
     >                      '' additional blocks'')')
            endif
c
            if(mblock.le.3) then
                write(6,*)' trouble: too few blocks for specified dr'
                write(6,*) ' stopping in subroutine radgrd'
                stop
             endif
        endif
c ---------------------------------------------------------------------
        i=1
        id=1
c
        ir=0
        r(1)=0.0
c
        do 50 j=1,nblock
           do 40 k=1,nsize
             i=i+1
             ir=ir+id
             r(i)=ir
   40      continue
           if(j.lt.mblock) id=id+id
   50   continue
c .....................................................................
c                                      muffin tin  radius
         jmt=mesh-3*nsize
         rfactr=rmt/r(jmt)
c .....................................................................
c                                     wigner-seitz radius
         do 60 i=1,mesh
            r(i)=r(i)*rfactr
            if(r(i).lt.rws) jws=i
  60     continue
         jws=jws+1
c .....................................................................
c                                      decrease mesh size if possible
         del1=r(2)
         delm=id*rfactr
  61     continue
c
         i=(jws-1)/nsize + 2
         mesh=min(ipmesh,i*nsize+1)
         nblock=mesh/nsize
c
c
c ---------------------------------------------------------------------
c           now print out information about radial grid
c
      write(6,'(1x,78(''='')/'' radial grid parameters: mesh size is'',
     >      i4,'' broken up into'',i3,'' blocks of'',i3,'' points''/
     >      1x,78(''-'')/''         initial increment is'',e16.9,
     >      '', which doubles until block'',i3/
     >      '' ending with the increment at'',e16.9/1x,78(''='')  )')
     >       mesh,nblock,nsize,del1,mblock,delm
c
c ---------------------------------------------------------------------
          mpr=20
          write(6,'(1x,i3,''.r='',f8.5,2x,i3,''.r='',f8.5,
     >    2x,i3,''.r='',f8.5,2x,i3,''.r='',f8.5,2x,i3,''.r='',f8.5)')
     >             (k,r(k),k=1,mesh,mpr)
          write(6,'(1x,78(''=''))')
          write(6,'(4x,''  muffin tin  radius r(mt)='',f8.5,
     >         7x,''  at i='',i3,'' r('',i3,'')='',f8.5)')
     >                     rmt,jmt,jmt,r(jmt)
          write(6,'(4x,'' wigner-sietz radius r(ws)='',f8.5,
     >      8x,'' at i='',i3,'' r('',i3,'')='',f8.5  /
     >      ''                               '',12x,
     >      ''and  at i='',i3,'' r('',i3,'')='',f8.5)')
     >         rws,jws-1,jws-1,r(jws-1),jws,jws,r(jws)
          write(6,'(1x,78(''=''))')
c ---------------------------------------------------------------------
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
         subroutine mtwsrd(icrys,alat,rmt,rws,oma,cay,ipr)
c =====================================================================
         implicit real*8 (a-h,o-z)
         parameter(three=3.,one=1.,third=one/three)
         parameter(pi=3.14159265358979)
c _____________________________________________________________________
c
c  routine sets muffin tin and wigner seitz radius, oma and cay
c   rmt is the muffin-tin radius
c   vol is the wigner-seitz volume
c   oma is the volume of the wigner-seitz cell outside the muffin tin
c   cay is the dimensionless constant which enters the madelung
c          term in the potential and total energy
c _____________________________________________________________________
         if(rmt.lt.0.01) rmt=alat*sqrt(icrys+one)/4.
         vol=2**icrys*(one/8.)
         if(rws.lt.0.01) rws= alat *  ( (3.*vol)/(4.*pi) ) ** third
         if(icrys.eq.1) then
             cay=2.41603350
         endif
         if(icrys.eq.2) then
             cay=2.0427605
         endif
         if(icrys.eq.3) then
             cay=1.5583429
         endif
         oma=vol-(4.*pi/3.)*(rmt/alat)**3
         if(ipr.lt.1) return
         write(6,'('' exiting mtwsrd: oma='',f12.8,3x,''cay='',f12.8)')
     >              oma,cay
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        subroutine rdpot(ipot,iunit,r,mesh,nsize,mblock,jmt,jws,
     >          ic,is,vr,v0,rho,rhoc,
     >          alat,zed,zcore,vws,ef,
     >          numbs,npbs,lbs,ebs,occbs,
     >          numsc,npsc,lsc,esc,occsc,ipr)
c =====================================================================
       implicit real*8 (a-h,o-z)
c
       parameter (ipmesh=641)
       parameter (ipspin=2,ipcomp=4,ipbst=20)
c
       dimension r(ipmesh),vr(ipmesh),rho(ipmesh),rhoc(ipmesh)
       dimension numbs(ipcomp,ipspin),numsc(ipcomp,ipspin)
       dimension npbs(ipbst,ipcomp,ipspin),ebs(ipbst,ipcomp,ipspin),
     >           lbs(ipbst,ipcomp,ipspin),occbs(ipbst,ipcomp,ipspin)
       dimension npsc(ipbst,ipcomp,ipspin),esc(ipbst,ipcomp,ipspin),
     >           lsc(ipbst,ipcomp,ipspin),occsc(ipbst,ipcomp,ipspin)
c _____________________________________________________________________
c
c       read the potentials
c       ===================
c
c      five (5)  options     ipot= -2,-1,0,1,2
c _____________________________________________________________________
c
c       ipot lt zero   atomic potentials,
c            -1   for spin restricted potentials
c            -2   for spin polarized potentials
c _____________________________________________________________________
c
c       ipot gt zero   log grid type --- old standard
c             1   for spin-restricted potentials
c             2   for spin polarized potentials
c _____________________________________________________________________
c
c .......    ipot equal zero   new standard      ......................
c                              ============
c
c  first line (i1) of pot deck tells if potential is spin-polarized
c _____________________________________________________________________
c
       if(abs(ipot).gt.2) then
             write(6,*) ' the value of ipot corresponds to an option'
             write(6,*) ' that has not been implemented: stopping'
             stop
       endif
c _____________________________________________________________________
c
c      if(ipot) 60,50,40
       if(ipot < 0) then
          goto 60
       else if (ipot == 0) then
          goto 50
       else
          goto 40
       endif
c _____________________________________________________________________
c
 40   call oldrdp(ipot,iunit,mesh,r,jmt,vr,v0,jws,zcore,vws,rho,rhoc,
     >   ef,numbs(ic,is),npbs(1,ic,is),lbs(1,ic,is),ebs(1,ic,is),
     >   numsc(ic,is),npsc(1,ic,is),lsc(1,ic,is),esc(1,ic,is) ,ipr)
      go to 999
c _____________________________________________________________________
c
 60   call atmrdp(ipot,iunit,mesh,r,vr,rho,numc,zcore,vws,
     >       npbs(1,ic,is),lbs(1,ic,is),ebs(1,ic,is),is,xmom,ipr)
c _____________________________________________________________________
c
      vws=0
      ef=0.0
      v0 = 0.0d0
      numsc(ic,is) = 0
      go to 999
c _____________________________________________________________________
c
 50   call stdrpt(iunit,zed,zcore,vws,alat,ef,
     >      mesh,nsize,mblock,jmt,jws,    r,vr,v0,rho,rhoc,
     >   numbs(ic,is),npbs(1,ic,is),lbs(1,ic,is),ebs(1,ic,is) ,
     >   numsc(ic,is),npsc(1,ic,is),lsc(1,ic,is),esc(1,ic,is) ,ipr)
c _____________________________________________________________________
c
  999 continue
c _____________________________________________________________________
c
        mpr=80
        if(ipr.ge.1) mpr=40
        if(ipr.ge.2) mpr=20
        if(ipr.ge.1) then
                   write(6,'(1x,78(''=''))')
                   write(6,'(1x,i3,''.r='',f8.5,2x,i3,''.r='',f8.5,
     >                2x,i3,''.r='',f8.5,2x,i3,''.r='',f8.5,2x,i3,
     >                ''.r='',f8.5)') (k,r(k),k=1,mesh,mpr)
        endif
c _____________________________________________________________________
c
        write(6,'(1x,78(''=''))')
        write(6,'(1x,i3,''.v='',f8.2,2x,i3,''.v='',f8.2,
     >               2x,i3,''.v='',f8.2,2x,i3,''.v='',f8.2,2x,i3,
     >               ''.v='',f8.2)') (k,vr(k),k=1,mesh,mpr)
        write(6,'(1x,78(''=''))')
c _____________________________________________________________________
c
        if(ipr.ge.1) then
                   write(6,'(1x,i3,''rh='',f8.5,2x,i3,''rh='',f8.5,
     >               2x,i3,''rh='',f8.5,2x,i3,''rh='',f8.5,2x,i3,
     >               ''rh='',f8.5)') (k,rho(k),k=1,mesh,mpr)
                   write(6,'(1x,78(''=''))')
c
                   write(6,'(1x,i3,''rc='',f8.5,2x,i3,''rc='',f8.5,
     >               2x,i3,''rc='',f8.5,2x,i3,''rc='',f8.5,2x,i3,
     >               ''rc='',f8.5)') (k,rhoc(k),k=1,mesh,mpr)
                    write(6,'(1x,78(''=''))')
         endif
c _____________________________________________________________________
c                                                atomic number
        ized=-.5*vr(1) + .0000001
        zed=ized
c _____________________________________________________________________
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine oldrdp(ipot,iunit,mesh,r,jmt,vr,v0,jws,zcore,vws,rho,
     >        rhoc,ef,numbs,npbs,lbs,ebs,numsc,npsc,lsc,esc,ipr )
c =====================================================================
       implicit real*8 (a-h,o-z)
       parameter (ipmesh=641,ipold=551,ipbst=20)
       dimension r(ipmesh),vr(ipmesh),rho(ipmesh),rhoc(ipmesh)
       dimension rold(ipold),vrold(ipold),
     >          rhold(ipold),corold(ipold)
       dimension ebs(ipbst),npbs(ipbst),lbs(ipbst)
       dimension esc(ipbst),npsc(ipbst),lsc(ipbst)
c _____________________________________________________________________
c
c <<<<<<<<<<<<<<<<<<<<<  old spin-polarized standard  >>>>>>>>>>>>>>>>>
c
c       ipot gt zero   log grid type --- old standard
c             1   for spin-restricted potentials
c             2   for spin polarized potentials
c _____________________________________________________________________
c
      nin=iunit
      nshrt=6
c
      read(nin,'(a80)') jttl2
      if(ipr.ge.2) write(nshrt,'(a80)') jttl2
      read(nin,'(f5.0,17x,f12.5,f5.0,e20.13)') azed,aa,zcd,ef
      if(ipr.ge.2)
     >   write(nshrt,'(f5.0,17x,f12.5,f5.0,e20.13)') azed,aa,zcd,ef
      read(nin,'(17x,2e20.13,i5)') xst,xmt,jmt0
      if(ipr.ge.2)  write(nshrt,'(17x,2e20.13,i5)') xst,xmt,jmt0
      read(nin,'(4e20.13)') (vrold(j),j=1,jmt0)
      read(nin,'(35x,e20.13)') v0
      if(ipr.ge.2) write(nshrt,'(35x,e20.13)') v0
c
      hh=(xmt-xst)/(jmt0-1)
      do 3 j=1,jmt0
        rold(j)=exp(xst+(j-1)*hh)
        vrold(j)=vrold(j)-v0*rold(j)
  3   continue
      do 33 j=jmt0+1,ipold
        rold(j)=exp(xst+(j-1)*hh)
        vrold(j)=0.0
 33   continue
c
      read(nin,'(i5,e20.13)') ntrho,vws
      vws=ipot*.5*vws
      if(ipr.ge.0) write(nshrt,'(i5,e20.13)') ntrho,vws
      if(ntrho.gt.0) read(nin,'(4e20.13)') ( rhold(j),j=1,ntrho)
      read(nin,'(3i5)') numbs,ntcor,numsc
      if(ipr.ge.1) write(nshrt,'('' numc='',i5,''  ntcor='',i5,
     >  ''  numcsm='',i5)') numbs,ntcor,numsc
c
      zcore=0
      if(numbs.gt.0) then
            do 20 ist=1,numbs
            read(nin,'(2i5,f12.5)') npbs(ist),lbs(ist),ebs(ist)
            zcore=zcore+2*lbs(ist)+1
   20       continue
      endif
      if(numsc.gt.0) then
            do 21 i=1,numsc
            read(nin,'(2i5,f12.5)')  npsc(i),lsc(i),esc(i)
   21       continue
      endif
      if(ntcor.gt.0) then
            read(nin,'(4e20.13)')   (corold(j),j=1,ntcor)
      endif
c
c .............................................................
c
      do 50 i=1,jmt
      vr(i)=ylag(r(i),rold, vrold, 0,3,jmt0,iex)
   50 continue
      if(ntrho.gt.0) then
          do 51 i=1,jws
          rho(i)=ylag(r(i),rold, rhold,0,3,ntrho,iey) /(3-ipot)
   51     continue
      endif
      if(ntcor.gt.0) then
          do 52 i=1,mesh
          rhoc(i)=ylag(r(i),rold,corold,0,3,ntcor,iey)/(3-ipot)
   52     continue
      endif
c
c .....................................................................
c
       return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine atmrdp(ipt,iunit,mesh,r,vr,rho,numc,zcore,vws,
     >       npbs,lbs,ebs,is,xmom,ipr )
c =====================================================================
       implicit real*8 (a-h,o-z)
       parameter (ipmesh=641,ipold=1001,ipbst=20)
c
       save
c
       dimension r(ipmesh),vr(ipmesh),rho(ipmesh)
       dimension rold(ipold),vrolup(ipold), rholup(ipold)
       dimension vroldn(ipold), rholdn(ipold)
       dimension ebs(ipbst),npbs(ipbst),lbs(ipbst)
       dimension ebsu(ipbst),npbsu(ipbst),lbsu(ipbst),occbsu(ipbst)
       dimension ebsd(ipbst),npbsd(ipbst),lbsd(ipbst),occbsd(ipbst)
c _____________________________________________________________________
c
c <<<<<<<<<<<<<<<<<<<<<  standard atomic potential  >>>>>>>>>>>>>>>>>>>
c
c       ipot lt zero   atomic potentials,
c            -1   for spin restricted potentials
c            -2   for spin polarized potentials
c _____________________________________________________________________
c
        ipot=abs(ipt)
        if(is.eq.1) then
c
             read(iunit,*) numbsu,(occbsu(i),i=1,ipbst)
             if(ipr.ge.2) write(6,*) numbsu,occbsu
             read(iunit,*)
     x    (npbsu(i),i=1,ipbst),(lbsu(i),i=1,ipbst),(ebsu(i),i=1,ipbst)
             if(ipr.ge.2) write(6,*) npbsu,lbsu,ebsu
            if(ipot.ne.1) then
             read(iunit,*) numbsd,(occbsd(i),i=1,ipbst)
             read(iunit,*)
     x    (npbsd(i),i=1,ipbst),(lbsd(i),i=1,ipbst),(ebsd(i),i=1,ipbst)
            endif
      read(iunit,'(4e20.13)')(vrolup(i),i=1,ipold),(vroldn(i),i=1,ipold)
             read(iunit,'(4e20.13)')(rold(i),i=1,ipold)
      read(iunit,'(4e20.13)')(rholup(i),i=1,ipold),(rholdn(i),i=1,ipold)
c
             limr=ipold
             do 40 i=1,mesh
                vr(i)=ylag(r(i),rold, vrolup, 0,3,limr,iex)
                rho(i)=ylag(r(i),rold, rholup,0,3,limr,iey)
   40        continue
c
             xmom=0
             xnumup=0
             xnumdn=0
             numbs=numbsu
             do 45 ist=1,numbsu
                npbs(ist)=npbsu(ist)
                lbs(ist)=lbsu(ist)
                ebs(ist)=ebsu(ist)
                xnumup=xnumup+occbsu(ist)
   45        continue
             zcore=xnumup
             if(abs(ipot).ne.2) return
             do 49 ist=1,numbsd
                xnumdn=xnumdn+occbsd(ist)
   49        continue
             xmom=xnumup-xnumdn
c
        else
c
             do 50 i=1,mesh
                vr(i)=ylag(r(i),rold, vroldn, 0,3,limr,iex)
                rho(i)=ylag(r(i),rold, rholdn,0,3,limr,iey)
   50        continue
c
             numbs=numbsd
             do 55 ist=1,numbsd
                npbs(ist)=npbsd(ist)
                lbs(ist)=lbsd(ist)
                ebs(ist)=ebsd(ist)
   55        continue
             zcore=xnumdn
c
      endif
c
c .....................................................................
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine stdrpt(iunit,       z1,zcore,vws,alat1,ef,
     >   mesh1,nsize1,mblck1,jmt1,jws1,    r,  v,vdif,chden,corden,
     >   ncore,npc,lc,ec,      nsemi,npsc,lsc,esc   ,ipr    )
c =====================================================================
      implicit real*8 (a-h,o-z)
      character*78 title
      character*1 lsymc,lsymsc
c
      parameter (iprold=641)
c
      real*8 z,vws,alat,ef,rmt
      real*8 z1,alat1,ef1,divr,vdif
c
      dimension lsymc(20),lsymsc(20)
      dimension npc(20),npsc(20),   lc(20),lsc(20)
      dimension ec(20),esc(20)
      dimension rold(iprold),vold(iprold),chold(iprold),corold(iprold)
      dimension r(mesh1),v(mesh1),chden(mesh1),corden(mesh1)
c _____________________________________________________________________
c
      read(iunit,'(a78)') title
      read(iunit,'(1x,f4.0,f15.12,1x,f9.5,1x,f14.11)')
     >          z,vws,alat,ef
      read(iunit,'(7i5,f15.10)')
     >          nv,nr,nc,mesh,nsize,mblock,jmt,rmt
      read(iunit,'(4e20.13)')
     >          (vold(i),i=1,nv)
      read(iunit,'(e20.13)')
     >          vdif
      if(nr.gt.0) read(iunit,'(4e20.13)')
     >          (chold(i),i=1,nr)
      read(iunit,'(2i5)')
     >          ncore,nsemi
      zcore=0
      if(ncore.gt.0) then
              do 20 i=1,ncore
              read(iunit,'(2i5,f12.5,4x,a1)')
     >                  npc(i),lc(i),ec(i),lsymc(i)
              zcore=zcore+2*lc(i)+1
   20         continue
      endif
      if(nsemi.gt.0) then
             do 21 i=1,nsemi
             read(iunit,'(2i5,f12.5,4x,a1)')
     >                  npsc(i),lsc(i),esc(i),lsymsc(i)
   21        continue
      endif
      read(iunit,'(4e20.13)')
     >          (corold(i),i=1,nc)
c _____________________________________________________________________
c
      ef1=ef
c
      i=1
      id=1
      ir=0
      rold(1)=0.0
      nblock=mesh/nsize
c
      do 50 j=1,nblock
          do 40 k=1,nsize
          i=i+1
          ir=ir+id
          rold(i)=ir
  40      continue
          if(j.lt.mblock) id=id+id
  50  continue
c
      divr=rmt/rold(jmt)
      do 60 i=1,mesh
          rold(i)=divr*rold(i)
          vold(i)=vold(i)-rold(i)*vdif
  60  continue
c
c ------------------------------------- check grid parameters and
c                                       interpolate if different
                                   interp=0
      if(mesh.ne.mesh1)            interp=1
      if(mblock.ne.mblck1)         interp=1
      if(nsize.ne.nsize1)          interp=1
      if(jmt.ne.jmt1)              interp=1
      if(abs(alat-alat1).gt.1.e-8) interp=1
c
      if(interp.eq.0) then
           if(ipr.ge.1)     write(6,
     >     '('' potential and charge density:  no interpolation'')')
           do 80 i=1,mesh
              v(i)=vold(i)
              chden(i)=chold(i)
              corden(i)=corold(i)
  80       continue
c
      else
           if(ipr.ge.1)     write(6,
     >      '('' potential and charge density:     interpolation'')')
           do 90 i=1,mesh1
              v(i)=ylag(r(i),rold,  vold,0,3,nv,iex1)
              chden(i)=ylag(r(i),rold, chold,0,3,nr,iex2)
              corden(i)=ylag(r(i),rold,corold,0,3,nc,iex3)
  90       continue
      endif
c _____________________________________________________________________
c
c ------------------------------------- now write out info
      mpr=80
      if(ipr.ge.1) mpr=40
      if(ipr.ge.2) mpr=20
      write(6,'('' potential and charge density file:''/a78)') title
      write(6,'(1x,f4.0,1x,f14.10,1x,f9.5,1x,f14.10)')
     >          z,vws,alat,ef
      write(6,'(6(1x,i4),1x,f14.10)')
     >          jmt1,jws1,mesh,nsize,mblock,jmt,rmt
      if(ipr.ge.2) write(6,'(4e20.13)')
     >          (v(i),i=1,mesh,mpr)
      if(ipr.ge.2) write(6,'(e20.13)')
     >          vdif
      if(ipr.ge.3)write(6,'(4e20.13)')
     >          (chden(i),i=1,mesh,mpr)
      if(ipr.ge.2) write(6,'(2i5)')
     >          ncore,nsemi
      if(ipr.ge.3) then
             do 100 i=1,ncore
             write(6,'(2i5,f12.5,1x,i2,1x,a1,6x,''core state'')')
     >           npc(i),lc(i),ec(i),npc(i),lsymc(i)
  100        continue
      endif
      if(ipr.ge.3) then
             do 101 i=1,nsemi
             write(6,'(2i5,f12.5,i3,1x,a1,'' semi-core state'')')
     >           npsc(i),lsc(i),esc(i),npsc(i),lsymsc(i)
  101        continue
      endif
      if(ipr.ge.3) write(6,'(4e20.13)')
     >           (corden(i),i=1,mesh,mpr)
c _____________________________________________________________________
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine stdwpt(iunit,ic,ncomp,conc,esym,ename,z,vws,alat,ef,
     >   mesh,nsize,mblock,jmt,jws,    r,  v,vdif,chden,corden,
     >   ncore,npc,lc,ec,      nsemi,npsc,lsc,esc,chmt,is,nspins )
c =====================================================================
      implicit real*8 (a-h,o-z)
      character*2  esym(ncomp)
      character*15 ename(ncomp)
      character*7  spn
      character*1  lsym
c
      real*8 z,vws,alat,ef,vdif
      dimension spn(3),lsym(0:3)
      dimension npc(20),npsc(20),   lc(20),lsc(20)
      dimension ec(20),esc(20)
      dimension conc(ncomp),r(mesh),v(mesh),chden(mesh),corden(mesh)
c
      data spn/'no spin','spin up','spin dn'/
      data lsym /'s','p','d','f'/
c _____________________________________________________________________
c
      if(iunit.eq.6) then
          mpr=20
          write(6,'(1x,78(''-''))')
          write(6,'(/5x,''potential and charge density''
     >              ,3x,''(format: std)'')')
      else
          mpr=1
      endif
c
      write(iunit,'(1x,a10,'':'',a7,'':'',4x,6( a2,''('',f6.2,'')''))')
     >          ename(ic),spn(is+nspins-1),(esym(i),conc(i),i=1,ncomp)
      write(iunit,'(1x,f4.0,f15.11,1x,f9.5,1x,f14.10,f15.2)')
     >          z,vws,alat,ef,chmt
      write(iunit,'(7(1x,i4),1x,f14.10)')
     >          jmt,jws,jws,mesh,nsize,mblock,jmt,r(jmt)
      write(iunit,'(4e20.13)')
     >          (v(i)+vdif*r(i),i=1,jmt,mpr)
      write(iunit,'(e20.13,5x,''muffin-tin zero'')')
     >          vdif
      write(iunit,'(4e20.13)')
     >          (chden(i),i=1,jws,mpr)
      write(iunit,'(2i5)')
     >          ncore,nsemi
c
      if(ncore.gt.0) then
           do 680 i=1,ncore
           write(iunit,'(2i5,f12.5,1x,i2,1x,a1,''      core state'')')
     >           npc(i),lc(i),ec(i),npc(i),lsym(lc(i))
  680      continue
      endif
c
      if(nsemi.gt.0) then
           do 690 i=1,nsemi
           write(iunit,'(2i5,f12.5,1x,i2,1x,a1,'' semi-core state'')')
     >           npsc(i),lsc(i),esc(i),npsc(i),lsym(lsc(i))
  690      continue
      endif
c
      write(iunit,'(4e20.13)')
     >          (corden(i),i=1,jws,mpr)
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine oldwpt(iunit,ic,ncomp,conc,esym,ename,z,vws,alat,ef,
     >   mesh,nsize,mblock,jmt,jws,rws,r,  v,vdif,chden,corden,
     >   zcor,ncore,npc,lc,ec,      nsemi,npsc,lsc,esc,chmt,is,nspins )
c =====================================================================
      use SortModule, only : HeapSort
      use ChemElementModule, only : MaxLenOfAtomName, getName
      use ChemElementModule, only : getCoreStateN, getCoreStateL
      implicit real*8 (a-h,o-z)
      character*2  esym(ncomp)
      character*15 ename(ncomp)
      character*7  spn
      character*1  lsym
      character (len=MaxLenOfAtomName) :: atname
      logical :: found_i
      integer :: zi
c
      parameter (irold=551,incrs=501,xstart=-11.1309674)
c
      real*8 z,vws,alat,ef,vdif
      real*8 rmt,xlast,h,xtmp,xst
      dimension spn(3),lsym(0:3), idx(20), et(20)
      dimension npc(20),npsc(20),   lc(20),lsc(20)
      dimension ec(20),esc(20)
      dimension rold(irold),tmp(irold)
      dimension conc(ncomp),r(mesh),v(mesh),chden(mesh),corden(mesh)
c
      data spn/'no spin','spin up','spin dn'/
      data lsym /'s','p','d','f'/
c _____________________________________________________________________
c
      if(iunit.eq.6) then
          mpr=20
          write(6,'(1x,78(''-''))')
          write(6,'(/5x,''potential and charge density''
     >              ,3x,''(format: old)'')')
      else
          mpr=1
      endif
c
      if (is == 1) then
         write(iunit,'(a,6(a,f5.2,2x))')
     >               ' nspin   vdif :::::: Potentials for ',
     >               (trim(ename(i)),conc(i),i=1,ncomp)
         write(iunit,'(i5,2x,d20.13)')nspins,vdif
      endif
      write(iunit,'(3a,f3.0,4x,a)')
     >      ' NEWSS: ',ename(ic),' z = ',z,spn(is+nspins-1)
      write(iunit,'(f5.0,17x,f12.5,f5.0,e20.13)')z,alat,zcor,ef
c
c now setup old-style grid .............
c
      zi = z
      atname = getName(zi)
        rmt=r(jmt)
        xlast=log(rmt)
        write(iunit,'(17x,2e20.13,2i5)') xstart,xlast,incrs,incrs
        h=(xlast-xstart)/(incrs-1)
        xtmp=(log(rws)-xstart)/h
        num=xtmp + 1
        xst=xstart - h
        vmuf=(is-1)*vdif
        do 400 i=1,num
        rold(i)=exp(xst+i*h)
  400   continue
c
        do 500 i=1,num
        tmp(i)=ylag(rold(i),r,v,0,3,jws,iex) + vmuf*rold(i)
  500   continue
c
        write(iunit,'(4e20.13)') (tmp(i),i=1,incrs,mpr)
        write(iunit,'(35x,e20.13,2x,a)') vdif,'muffin-tin zero'
        write(iunit,'(i5,2e20.13)') num,(3-nspins)*vws,(3-nspins)*chmt
        do 600 i=1,num
        tmp(i)=(3-nspins)*ylag(rold(i),r,chden,0,3,jws,iex)
  600   continue
        write(iunit,'(4e20.13)') (tmp(i),i=1,num,mpr)
c
        n = 0
        do i = 1,ncore
           if (lc(i) == 0) then
              n = n + 1
           else
              n = n + 2
           endif
        enddo
        do i = 1,nsemi
           if (lsc(i) == 0) then
              n = n + 1
           else
              n = n + 2
           endif
        enddo
        write(iunit,'(2i5)') n,num
c       write(iunit,'(3i5)') ncore,num,nsemi
ccc     write(iunit,'(2i5)') ncore+nsemi,num
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  The following lines are modified so that the core states
c  are in an order compatible with ChemElementModule database
c  by Yang Wang, Nov 14, 2015
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        n = 1
        if(ncore.gt.0) then
           do i=1,ncore
              et(i) = ec(i)
           enddo
cccc       call HeapSort(ncore,et,idx)
           do j=1,ncore
cccc          i = idx(j)
              i = 0; found_i = .false.
              do while (i < ncore)
                 i = i + 1
                 if (npc(i)== getCoreStateN(atname,n) .and.
     >               lc(i) == getCoreStateL(atname,n)) then
                    found_i = .true.
                    exit
                 endif
              enddo
              if (found_i) then
                 if (lc(i) == 0) then
                    write(iunit,'(3i5,f12.5,1x,i2,1x,a1,6x,a)')
     >                 npc(i),lc(i),-1,ec(i),npc(i),lsym(lc(i)),
     >                 'core state'
                    n = n + 1
                 else
                    write(iunit,'(3i5,f12.5,1x,i2,1x,a1,6x,a)')
     >                 npc(i),lc(i),lc(i),ec(i),npc(i),lsym(lc(i)),
     >                 'core state'
                    write(iunit,'(3i5,f12.5,1x,i2,1x,a1,6x,a)')
     >                 npc(i),lc(i),-lc(i)-1,ec(i),npc(i),lsym(lc(i)),
     >                 'core state'
c                   write(iunit,'(3i5,f12.5,1x,i2,1x,a1,6x,a)')
c    >                 npc(i),lc(i),ec(i),npc(i),lsym(lc(i)),'core state'
                    n = n + 2
                 endif
              else
                 write(6,'(a)')
     >                 'A deep-core state is not found in database'
                 stop 'Error!'
              endif
           enddo
        endif
c
c       I took out the nsc and lsc comparison against the ChemElementModule 
c       database comparison due to the problem that the rule of the
c       semi-core configuration order is difficult to determine. For
c       exampl, for La, no 4f state under 5s in the semi-core states.
c       -Yang Wang, 05/19/2018............
        if(nsemi.gt.0) then
           do i=1,nsemi
              et(i) = esc(i)
           enddo
           call HeapSort(nsemi,et,idx)
           n = 0
           do j=1,nsemi
              i = idx(j)
              if (lsc(i) == 0) then
                 write(iunit,'(3i5,f12.5,1x,i2,1x,a1,1x,a)')
     >                 npsc(i),lsc(i),-1,esc(i),npsc(i),
     >                 lsym(lsc(i)),'semi-core state'
              else
                 write(iunit,'(3i5,f12.5,1x,i2,1x,a1,1x,a)')
     >                 npsc(i),lsc(i),lsc(i),esc(i),npsc(i),
     >                 lsym(lsc(i)),'semi-core state'
                 write(iunit,'(3i5,f12.5,1x,i2,1x,a1,1x,a)')
     >                 npsc(i),lsc(i),-lsc(i)-1,esc(i),npsc(i),
     >                 lsym(lsc(i)),'semi-core state'
              endif
           enddo
        endif
c
        do 700 i=1,num
        tmp(i)=(3-nspins)*ylag(rold(i),r,corden,0,3,jws,iex)
  700   continue
        write(iunit,'(4e20.13)') (tmp(i),i=1,num,mpr)
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine genpte(vnew,vdif,vold,vdold,
     >    ncomp,nspins,mesh,nsize,nblock,
     >    jmt,jws,rmt,rws,r,rho,rhoc,rmsrho,
     >    conc,zed,valtot,zcore,qmt,xmmt,
     >    cay,oma,icrys,alat,iexc,  qsum,esumv,itoten,ipr)
c =====================================================================
c                                             written  aug 88     fjp
c
c    generate new potential and calculates terms for the total energy
c    includes corrections to make the energy variational
c _____________________________________________________________________
c
      implicit real*8 (a-h,o-z)
c
      parameter (ipmesh=641,ipspin=2,ipcomp=4)
      parameter (one=1.,three=3.,third=one/three)
      parameter (pi=3.14159265358979,hr=1.,onemd=.999999999,rsmax=100.)
      parameter (onepd=1.0001)
c
      dimension conc(ipcomp),zed(ipcomp),r(ipmesh)
      dimension valtot(ipcomp,ipspin),zcore(ipcomp,ipspin)
      dimension qsum(ipcomp,ipspin),
     >           qs1(ipcomp,ipspin),qs2(ipcomp,ipspin)
      dimension rmspot(ipcomp,ipspin),dpotmx(ipcomp,ipspin)
      dimension term1a(ipcomp,ipspin),term1b(ipcomp,ipspin)
      dimension term2(ipcomp,ipspin),term3(ipcomp,ipspin)
      dimension correc(ipcomp,ipspin),cor3pv(ipcomp,ipspin)
      dimension ptrm2a(ipcomp,ipspin),ptrm2b(ipcomp,ipspin)
      dimension bt(ipcomp,ipspin),ptrm3(ipcomp,ipspin)
      dimension enxc(ipmesh,ipcomp),vmt(ipspin),vxcmt(ipspin)
      dimension   vx(ipmesh,ipcomp,ipspin)
      dimension vnew(ipmesh,ipcomp,ipspin)
      dimension vold(ipmesh,ipcomp,ipspin)
      dimension  rho(ipmesh,ipcomp,ipspin)
      dimension rhoc(ipmesh,ipcomp,ipspin)
      dimension rhov(ipmesh,ipcomp,ipspin)
      dimension   q1(ipmesh,ipcomp,ipspin)
      dimension   zr(ipmesh,ipcomp),q2(ipmesh,ipcomp)
      dimension  tmp(ipmesh),tmpi(ipmesh)
      dimension vh(ipmesh,ipcomp),tmpr(ipcomp*ipspin)
      dimension ztot(ipcomp,ipspin),qmt(ipcomp),xmmt(ipcomp)
      dimension zint(ipcomp,ipspin),qint(ipcomp),xmint(ipcomp)
      dimension ezpt(ipcomp,ipspin),tpzpt(ipcomp,ipspin)
      dimension thrpv(ipcomp,ipspin),esumv(ipcomp,ipspin)
      dimension ten(ipcomp,ipspin),rmsrho(ipcomp,ipspin)
c ____________________________________________________________________
c
c     generates r*v in rydbergs  change hr to .5 for hartrees
c ____________________________________________________________________
c
c  definitions of variables     in      <mt>,   <ws>,  <interstitial>
c         charge per spin per atom   q1(jmt),    ztot,    zint
c         charge per atom                qmt,             qint
c         moment per atom               xmmt,             xmint
c         average charge                tqmt, tz=tzed,    tqint
c         average moment                tmmt,    tmom,    tmint
c ____________________________________________________________________
c
c    valtot is the valence charge in the wigner-seitz  cell
c    q1(i): integrated charge per spin per atom within radius r(i)
c    q2(i): same as above except weighted by 1/r in the kernal
c ____________________________________________________________________
      tqmt=0.0
      tzed=0.0
      tz=0.0
      tqint=0.0
      tmmt=0.0
      tmom=0.0
      tmint=0.0
c ____________________________________________________________________
       rmsrh=0.0
       numrms=0
c
      do 50 ic=1,ncomp
        do 40 is=1,nspins
           numrms=numrms+1
           rmsrh=rmsrh+rmsrho(ic,is)
           call intqud(rho(1,ic,is),q1(1,ic,is),r,nblock,nsize)
           ztot(ic,is)=valtot(ic,is)+zcore(ic,is)
           zint(ic,is)=ztot(ic,is)-q1(jmt,ic,is)
           write(6,'('' ic,is,valtot,zcore '',2i5,2d12.5)')
     >     ic,is,valtot(ic,is),zcore(ic,is)
  40    continue
c
c
        qmt(ic)=q1(jmt,ic,1)+q1(jmt,ic,nspins)
        write(6,'('' qmtup and qmtdn '',2d12.5)') 
     >  q1(jmt,ic,1),q1(jmt,ic,nspins)
        xmmt(ic)=q1(jmt,ic,1)-q1(jmt,ic,nspins)
        qint(ic)=zint(ic,1)+zint(ic,nspins)
        write(6,'('' ztotup and ztotdn '',2d12.5)') 
     >  ztot(ic,1),ztot(ic,nspins)
        write(6,'('' zintup and zintdn '',2d12.5)') 
     >  zint(ic,1),zint(ic,nspins)
        write(6,'('' qintup and qintdn '',2d12.5)') 
     >  zint(ic,1),zint(ic,nspins)
        xmint(ic)=zint(ic,1)-zint(ic,nspins)
        tz=tz+conc(ic)*(ztot(ic,1)+ztot(ic,nspins))
        tzed=tzed+conc(ic)*zed(ic)
        tmom=tmom+conc(ic)*(ztot(ic,1)-ztot(ic,nspins))
        tmmt=tmmt+conc(ic)*xmmt(ic)
        tqmt=tqmt+conc(ic)*qmt(ic)
        tqint=tqint+conc(ic)*qint(ic)
        tmint=tmint+conc(ic)*xmint(ic)
  50  continue
      rmsrh=rmsrh/numrms
c ____________________________________________________________________
c
      write(6,'(1x,78(''_''))')
      write(6,'(4x,''conc'',t13,''at. no.'',t23,''wsc chg'',
     >        t33,''wsc mom'',t43,''mt chg'',t53,''mt mom'',
     >        t63,''out chg'',t73,''out mom'')')
       do 60 ic=1,ncomp
          write(6,'(1x,f8.5,7f10.5)')
     >    conc(ic),zed(ic),ztot(ic,1)+ztot(ic,nspins),
     >                  ztot(ic,1)-ztot(ic,nspins),
     >                  qmt(ic),xmmt(ic),qint(ic),xmint(ic)
   60  continue
c
       if(ncomp.gt.0) then
               write(6,'(1x,''averages:'',f9.5,6f10.5)')
     >                 tzed,tz,tmom,tqmt,tmmt,tqint,tmint
       endif
       write(6,'(1x,78(''_''))')
       if(abs(tzed-tz).gt.1.e-8) then
               write(6,*)' stopping: charge discrepency'
               do is=1,nspins
               do ic=1,ncomp
               do i=1,jmt,40
               write(6,*)' ic=',ic,'  is=',is,'  rh=',rho(i,ic,is)
               enddo
               enddo
               enddo
               stop
         endif
c ____________________________________________________________________
c
c     now calculate the interstitial potential vmt =  vxmt + vhout
c
      vhout=-2*tqint*cay/alat
      rsout=alat*(3*oma/(4*pi*abs(tqint)))**third
      rsout=min(rsout,rsmax)
      dzout=min(onemd,max(-onemd,tmint/tqint))
      do 100 is=1,nspins
        sp=3-2*is
        vxcmt(is)=alpha2(rsout,dzout,sp,iexc,excout)
        vmt(is)=vxcmt(is)+vhout
 100  continue
      vdif=vmt(nspins)-vmt(1)
c _____________________________________________________________________
      write(6,'(1x,''interstitial potential: v(mad.)='',f8.5,
     >          3x,''ex-corr: '',f8.5,''(up)'',2x,f8.5,''(dn)'')')
     >         vhout,(vxcmt(is),is=1,nspins)
      write(6,'(1x,''interstitial potential:         '',
     >         11x,''  total: '',f8.5,''(up)'',2x,f8.5,''(dn)'')')
     >         (vmt(is),is=1,nspins)
      if(nspins.gt.1) write(6,'(12x,''difference in muffin-tin'',
     >         '' zeroes:  vdif='',f12.8)') vdif
c ____________________________________________________________________
c     generate the potential
c
c
      do 150 ic=1,ncomp
         do 125 i=2,mesh
         zr(i,ic)=q1(i,ic,1)+q1(i,ic,nspins)
         tmp(i)=(rho(i,ic,1)+rho(i,ic,nspins))/r(i)
 125     continue
      tmp(1)=0.0
      zr(1,ic)=0.0
      call intqud(tmp,q2(1,ic),r,nblock,nsize)
 150  continue
c ____________________________________________________________________
c
      do 400 ic=1,ncomp
          do 225 is=1,nspins
          dpotmx(ic,is)=0.0
          rmspot(ic,is)=0.0
  225     continue
          vnew(1,ic,1)=-2*hr*zed(ic)
          vnew(1,ic,nspins)=vnew(1,ic,1)
          do 300 i=2,jmt
             totrh= rho(i,ic,1)+rho(i,ic,nspins)
             dz=(rho(i,ic,nspins)-rho(i,ic,1))/totrh
             dz=min(onemd,max(-onepd,dz))
             rs=((3*r(i)**2)/abs(totrh))**third
             rs=min(rs,rsmax)
             vh(i,ic)=   2 *   (  (zr(i,ic)-zed(ic))
     >                   + r(i)*(q2(jmt,ic)-q2(i,ic)) )
             do 250 is=1,nspins
                sp=3-2*is
                vx(i,ic,is)=alpha2(rs,dz,sp,iexc,enxc(i,ic))
                vnew(i,ic,is)=hr*(vh(i,ic)+r(i)*(vx(i,ic,is)-vmt(is)))
                dif2=(vnew(i,ic,is)-vold(i,ic,is))**2
                dpotmx(ic,is)=max(dpotmx(ic,is),dif2)
                rmspot(ic,is)=rmspot(ic,is)+dif2
  250        continue
  300     continue
          do 380 is=1,nspins
             dpotmx(ic,is)=sqrt(dpotmx(ic,is))
             rmspot(ic,is)=sqrt(rmspot(ic,is)/jmt)
 380      continue
 400  continue
c ____________________________________________________________________
c
      if(itoten.ne.1)  then
c
c ____________________________________________________________________
         mpr=80
         if(ipr.ge.1) mpr=40
         if(ipr.ge.2) mpr=20
         write(6,'(1x,78(''=''))')
         write(6,'(''   new potentials:'')')
         do i=1,jmt,mpr
            npr=0
            do ic=1,ncomp
            do is=1,nspins
               npr=npr+1
               tmpr(npr)=vnew(i,ic,is)
            enddo
            enddo
           jpr=(npr-1)/4
           do kpr=0,jpr
             if(kpr.eq.jpr) then
                jend=mod(npr-1,4)+1
             else
                jend=4
             endif
             write(6,'(1x,''i='',i4,3x,''r='',f10.7,2x,4f13.7)')
     >          i,r(i),(tmpr(kpr*4+j),j=1,jend)
           enddo
         enddo
c
         write(6,'(1x,78(''-''))')
         npr=0
         do ic=1,ncomp
         do is=1,nspins
           npr=npr+1
           tmpr(npr)=rmspot(ic,is)
         enddo
         enddo
         jpr=(npr-1)/4
         do kpr=0,jpr
         if(kpr.eq.jpr) then
              jend=mod(npr-1,4)+1
         else
              jend=4
         endif
         write(6,'(2x,''root-mean-squared diff'',4f13.7)')
     >         (tmpr(kpr*4+j),j=1,jend)
         enddo
         npr=0
         do ic=1,ncomp
         do is=1,nspins
           npr=npr+1
           tmpr(npr)=dpotmx(ic,is)
         enddo
         enddo
         jpr=(npr-1)/4
         do kpr=0,jpr
         if(kpr.eq.jpr) then
              jend=mod(npr-1,4)+1
         else
              jend=4
         endif
         write(6,'(6x,''maximum difference'',4f13.7)')
     >         (tmpr(kpr*4+j),j=1,jend)
         enddo
c
c
c
c
c
c
       else
c ________________________________            =========================
c                                             total energy and pressure
c                                             total energy and pressure
c                                             total energy and pressure
c                                             total energy and pressure
c ________________________________            =========================
c
         jblock=jmt/nsize
         do 750 is=1,nspins
            sp=3-2*is
            do ic=1,ncomp
               tmp(1)=0.0
               do i=2,jmt
               rhov(i,ic,is)=rho(i,ic,is)-rhoc(i,ic,is)
               enddo
c _____________________________________________________________________
c     - int 2pi r3dr rhocore dv/dr  =  int 2pi r2dr (v-d(rv)/dr)
c _____________________________________________________________________
c      where v is vold   and   vold = (vnew-vx) + (vold-vnew+vx)
c
            do 512 i=2,jmt
            tmp(i)=(vold(i,ic,is)-vnew(i,ic,is))+r(i)*vx(i,ic,is)
  512       continue
            call derv5(tmp,tmpi,r,jmt)
            do 515 i=2,jmt
            tmp(i)=rhoc(i,ic,is)*(tmpi(i)-vold(i,ic,is)/r(i)
     >                 +2*(q2(jmt,ic)-q2(i,ic)) - vmt(is) )
  515       continue
            call intqud(tmp,tmpi,r,jblock,nsize)
            qs1(ic,is)=-.5*tmpi(jmt)
c
            call derv5(vold(1,ic,is),tmpi,r,jmt)
            do 522 i=2,jmt
            tmp(i)=rhoc(i,ic,is)*(vold(i,ic,is)/r(i)-tmpi(i))
  522       continue
            call intqud(tmp,tmpi,r,jblock,nsize)
            qs2(ic,is)=.5*tmpi(jmt)
c
c            do 504 i=2,jmt
c 504        tmp(i)=vold(i,ic,is)/r(i)
c            call derv5(tmp,tmpi,r,jmt)
c            do 506 i=2,jmt
c 506        tmp(i)=rhoc(i,ic,is)*r(i)*tmpi(i)
c            call intqud(tmp,tmpi,r,jblock,nsize)
c            qs3(ic,is)=-.5*tmpi(jmt)
c _____________________________________________________________________
c
            bt(ic,is)=2*q2(jmt,ic)-vmt(is)
            term2(ic,is)= .5*(qint(ic)+sp*xmint(ic))*
     >                       (excout-tqint*(cay/alat))
     >        -r(jmt)*rho(jmt,ic,is)*(enxc(jmt,ic)-vx(jmt,ic,is))
            ptrm3(ic,is)=-.5*(qint(ic)+sp*xmint(ic))*
     >            (3*(excout-vxcmt(is))+tqint*(cay/alat))
     >        -r(jmt)*rho(jmt,ic,is)*(enxc(jmt,ic)-vx(jmt,ic,is))
c _____________________________________________________________________
c
c                                 - int dr 4pi r2 rhoval d(rv)/dr
            do 610 i=2,jmt
            tmp(i)=rhov(i,ic,is)*(2*q2(i,ic)-bt(ic,is))
  610       continue
c
            call intqud(tmp,tmpi,r,jblock,nsize)
            term1a(ic,is)=tmpi(jmt)
            call derv5(vx(1,ic,is),tmpi,r,jmt)
            do 620 i=2,jmt
            tmp(i)=-rhov(i,ic,is)*(vx(i,ic,is)+r(i)*tmpi(i))
  620       continue
c
            call intqud(tmp,tmpi,r,jblock,nsize)
            term1a(ic,is)=term1a(ic,is)+tmpi(jmt)
            call derv5(vnew(1,ic,is),tmpi,r,jmt)
            do 625 i=2,jmt
            tmp(i)=-rhov(i,ic,is)*tmpi(i)
  625       continue
            call intqud(tmp,tmpi,r,jblock,nsize)
            term1b(ic,is)=tmpi(jmt)
c _____________________________________________________________________
c
c                                   int dr 4pi r2 rho (4exc-3vxc)
            do 630 i=2,jmt
            tmp(i)=rho(i,ic,is)*(4*enxc(i,ic)-3*vx(i,ic,is))
  630       continue
            call intqud(tmp,tmpi,r,jblock,nsize)
            term3(ic,is)=tmpi(jmt)
c _____________________________________________________________________
c
c                                 - int dr 4pi r rhoval d(r2v)/dr
            do 640 i=2,jmt
            tmp(i)=r(i)**2*vx(i,ic,is)
  640       continue
            call derv5(tmp,tmpi,r,jmt)
            do 650 i=2,jmt
            tmp(i)=-rhov(i,ic,is)*(tmpi(i)+2*((zr(i,ic)-zed(ic))
     >            -r(i)*(2*q2(i,ic)-bt(ic,is)) ) )  /r(i)
  650       continue
            call intqud(tmp,tmpi,r,jblock,nsize)
            ptrm2a(ic,is)=tmpi(jmt)
c
            do 670 i=2,jmt
            tmp(i)=r(i)*vnew(i,ic,is)
  670       continue
            call derv5(tmp,tmpi,r,jmt)
            do 675 i=2,jmt
            tmp(i)=-rhov(i,ic,is)*tmpi(i)/r(i)
  675       continue
            call intqud(tmp,tmpi,r,jblock,nsize)
            ptrm2b(ic,is)=tmpi(jmt)
c _____________________________________________________________________
c
c            terms that make energy variational:  correc
c               and corrections to the pressure:  cor3pv
c
            do 710 i=2,jmt
            tmp(i)=(vold(i,ic,is)-vnew(i,ic,is))/r(i)
  710       continue
c
            call derv5(tmp,tmpi,r,jmt)
c
            do 715 i=2,jmt
            tmp(i)=r(i)*rhoc(i,ic,is)*tmpi(i)
  715       continue
c
            call intqud(tmp,tmpi,r,jblock,nsize)
            correc(ic,is)=tmpi(jmt)
c
            do 720 i=2,jmt
            tmp(i)=-rhov(i,ic,is)*(vold(i,ic,is)-vnew(i,ic,is))/r(i)
  720       continue
            call intqud(tmp,tmpi,r,jblock,nsize)
c
            cor3pv(ic,is)=correc(ic,is)+2.*tmpi(jmt)
            correc(ic,is)=correc(ic,is)+tmpi(jmt)
c
         enddo ! ic
  750    continue
c ____________________________________________________________________
c
         call zeropt(ezpt,tpzpt,alat,icrys,zed,ncomp,nspins)
c _____________________________________________________________________
c
         if(ipr.gt.1)  then
           do 850 ic=1,ncomp
           do is=1,nspins
             xqs=qsum(ic,is)
             xq1=qs1(ic,is)
             qd=xqs-xq1
             write(6,'(1x,78(''_''))')
             write(6,'(3x,''cmp '',i1,'' spn '',i1,''   qsum;'',
     >       ''  from core:'',f16.7,'' (ryd)''/19x, ''from integral1:''
     >       ,f16.7/14x,''difference (core-1):'',f16.7 )')
     >       ic,is,xqs,xq1,qd
             write(6,'(12x,''total energy:  term1a:'',f16.7/27x,
     >       ''term1b:'',f16.7/  15x,''difference (1a-1b):'',f16.7)')
     >       term1a(ic,is),term1b(ic,is),term1a(ic,is)-term1b(ic,is)
             write(6,'(12x,''total energy:   term2:'',f16.7)')
     >             term2(ic,is)
             write(6,'(12x,''total energy:   term3:'',f16.7)')
     >             term3(ic,is)
             write(6,'(12x,''total energy: zero-pt:'',f16.7)')
     >             ezpt(ic,is)
             write(6,'(12x,''total energy:  correc:'',f16.7)')
     >             correc(ic,is)
             write(6,'(16x,''pressure:  ptrm2a:'',f16.7/26x,
     >       '' ptrm2b:'',f16.7/  15x,''difference (2a-2b):'',f16.7)')
     >       ptrm2a(ic,is),ptrm2b(ic,is),ptrm2a(ic,is)-ptrm2b(ic,is)
             write(6,'(16x,''pressure:   ptrm3:'',f16.7)') ptrm3(ic,is)
             write(6,'(16x,''pressure: zero-pt:'',f16.7)') tpzpt(ic,is)
             write(6,'(16x,''pressure:  cor3pv:'',f16.7)')cor3pv(ic,is)
          enddo
  850     continue
         endif
c ____________________________________________________________________
         tpv=0.0
         toten=0.0
         do 875 ic=1,ncomp
             do 870 is=1,nspins
                ten(ic,is)=esumv(ic,is)+term1a(ic,is)+term2(ic,is)+
     >            term3(ic,is)+qsum(ic,is)+correc(ic,is)+ezpt(ic,is)
c
                  thrpv(ic,is)= 2.*esumv(ic,is)+ptrm2a(ic,is)+
     >                    ptrm3(ic,is)+tpzpt(ic,is)+cor3pv(ic,is)
  870        continue
c
             tpv=tpv + conc(ic)*(thrpv(ic,1)+thrpv(ic,nspins))
             toten=toten + conc(ic)*(ten(ic,1)+ten(ic,nspins))
c
 875     continue
c ____________________________________________________________________
c
         write(6,828)
         write(6,829)
         do 231 ic=1,ncomp
         do is=1,nspins
            write(6,830) ic,is,esumv(ic,is),qsum(ic,is),
     >             term1a(ic,is),term3(ic,is),term2(ic,is)
         enddo
  231    continue
c
         write(6,831)
         do 232 ic=1,ncomp
         do is=1,nspins
            write(6,833) ic,is,ezpt(ic,is),ptrm2a(ic,is),ptrm3(ic,is),
     >        thrpv(ic,is),ten(ic,is)
         enddo
  232    continue
         write(6,828)
         write(6,832) tpv,toten
         write(6,828)
c ____________________________________________________________________
c
         call keepem(icrys,alat,conc,zed,toten,tpv,rmsrh,ncomp,nspins)
c ____________________________________________________________________
      endif
c
  828  format(1x,78('*'))
  829  format(1x,'*cmp spn',t15,'esum',t28,'qsum',
     > t43,'term1',t56,'ex-corr',t71,'misc',t79,'*')
  830  format(1x,'*',2i3,3f14.5,f12.5,f14.5,t79,'*')
  833  format(1x,'*',2i3,3f14.5,f13.6,f14.6,t79,'*')
  831  format(1x,'*cmp spn',t16,'zero pt',t30,' ptrm2',
     >  t45,'misc',t56,'3pv (ry)',t67,'energy (ry)',t79,'*')
  832  format(1x,'*',3x,'averaged over components:',
     > 3x,'3pv=',f14.6,4x,'e(ry)=',f16.6,t79,'*')
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      function alpha2 (rs,dz,sp,iexch,exchg)
c =====================================================================
c
      implicit real*8 (a-h,o-z)
      parameter(one=1.,two=2.,three=3.,thrd=1./three,for3=4./three)
      dimension vxx(2),vxcs(2),g(3),dg(3),tbq(3),tbxq(3)
      dimension bxx(3),a(3),b(3),c(3),q(3),x0(3),bb(3),cx(3)
      save
c
c     data for von barth-hedin
c
      data ccp,rp,ccf,rf/0.045,21.,0.0225,52.916682/
c
c     data for vosko-wilk-nusair
c
      data incof/0/
      data  a/-0.033774,0.0621814,0.0310907/
      data  b/1.13107,3.72744,7.06042/
      data  c/13.0045,12.9352,18.0578/
      data x0/-0.0047584,-0.10498,-0.325/
      data cst,aip/1.92366105,0.91633059/
      data fnot,bip/1.70992095,0.25992105/
c
c  the following are constants needed to obtain potential
c  which are in g.s. painter's paper
c  =====given here for check=====(generated below)
c
c     data q/.7123108918e+01,.6151990820e+01,.4730926910e+01/
c     data bxx/-.4140337943e-03,-.3116760868e-01,-.1446006102/
c     data tbq/.3175776232e+00,.1211783343e+01,.2984793524e+01/
c     data tbxq/.3149055315e+00,.1143525764e+01,.2710005934e+01/
c     data bb/.4526137444e+01,.1534828576e+02,.3194948257e+02/
c
      go to (10,20) iexch
c
c     von barth-hedin  exch-corr potential
c     j. phys. c5,1629(1972)
c
  10  continue
      fm=two**for3-two
c     print *,' dz = ',dz
      fdz=( (one+dz)**for3 + (one-dz)**for3 - two)/fm
      ex=-0.91633/rs
      exf=ex*two**thrd
      xp=rs/rp
      xf=rs/rf
      gp=(one+xp**3)*log(one+one/xp) - xp*xp +xp/two - thrd
      gf=(one+xf**3)*log(one+one/xf) - xf*xf +xf/two - thrd
      exc=ex-ccp*gp
      excf=exf-ccf*gf
      dedz=for3*(excf-exc)*((one+dz)**thrd - (one-dz)**thrd)/fm
      gpp=three*xp*xp*log(one+one/xp)-one/xp + 1.5 - three*xp
      gfp=three*xf*xf*log(one+one/xf)-one/xf + 1.5 - three*xf
      depd=-ex/rs-ccp/rp*gpp
      defd=-exf/rs-ccf/rf*gfp
      decd=depd+(defd-depd)*fdz
c exchange-correlation energy
      exchg= exc + (excf-exc)*fdz
c exchange-correlation potential
      alpha2=exc+(excf-exc)*fdz - rs*decd*thrd + sp*(one-sp*dz)*dedz
      return
c
c
  20  continue
      if(incof.ne.0) go to 30
      incof=2
c
c   vosk0-wilk-nusair exch-corr potential
c    taken from g.s. painter
c   phys. rev. b24 4264,1981
c
c   generate constant coefficients for the parameterization (v-w-n)
c
      do 1 i=1,3
      cx(i)= x0(i)**2 + b(i)*x0(i) + c(i)
      bfc= 4*c(i) - b(i)**2
      q(i)= sqrt(bfc)
      bxx(i)= b(i)*x0(i)/cx(i)
      tbq(i)= 2*b(i)/q(i)
      tbxq(i)= tbq(i) + 4*x0(i)/q(i)
      bb(i)= 4*b(i)*( 1 - x0(i)*(b(i) + 2*x0(i))/cx(i) )
  1   continue
c
  30  continue
      zp1= 1 + dz
      zm1= 1 - dz
      xr=sqrt(rs)
      pex= -aip/rs
      xrsq= rs
c
c   generate g(i)=alpha,epsilon fct.s
c     and their derivatives dg(i)
c    1=alpha(spin stiffness)  2=ecp  3=ecf
c
      do 2 i=1,3
         qi=q(i)
         txb= 2*xr + b(i)
         fx= xrsq + xr*b(i) + c(i)
         arct= atan2(qi,txb)
         dxs= (xr-x0(i))**2/fx
c
         g(i)=a(i)*( log(xrsq/fx) + tbq(i)*arct-bxx(i)*
     >                 (log(dxs) + tbxq(i)*arct) )
c
         dg(i)=a(i)*( 2./xr - txb/fx -
     >      bxx(i)*(2./(xr-x0(i))-txb/fx) - bb(i)/(qi**2 + txb**2) )
c
   2  continue
c
      ecp=g(2)
      zp3=zp1**thrd
      zm3=zm1**thrd
      zp3m3=zp3-zm3
c part of last term in vx   eq(13)
      fx1=.5*for3*pex*zp3m3
      z4= dz**4
      fz= cst*(zp1**for3 + zm1**for3 - 2.)
      beta= fnot*( g(3)-g(2) )/g(1) - 1.
      ec= ecp + fz*g(1)*( 1 + z4*beta )/fnot
      ex= pex*( 1 + fz*bip )
      f3ex= for3*ex
c echange-correlation energy
      exchg= ec + ex
c exchange potential
      vxx(1)= f3ex + fx1*zm1
      vxx(2)= f3ex - fx1*zp1
c correlation potential
      vcc= ec - xr*( (1 - z4*fz)*dg(2) + z4*fz*dg(3)
     >                  +  (1 - z4)*fz*dg(1)/fnot  )/6
c
      facc= 4*g(1)*( dz**3*fz*beta + (1 + beta*z4 )*zp3m3/(6*bip) )/fnot
c
c exch-corr. potential for each spin as called in newpot
c
      vxcs(1)= vcc + zm1*facc + vxx(1)
      vxcs(2)= vcc - zp1*facc + vxx(2)
c
      if( sp.ge.0 ) alpha2= vxcs(1)
      if( sp.lt.0 ) alpha2= vxcs(2)
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine zeropt(ezpt,tpzpt,alat,icrys,z,nsites,mspn)
c =====================================================================
       implicit real*8 (a-h,o-z)
c
c ---------------------------------------------------------------------
       parameter (ipcomp=4,ipspin=2)
c ---------------------------------------------------------------------
c
       dimension z(ipcomp),ezpt(ipcomp,ipspin),tpzpt(ipcomp,ipspin)
       dimension idebye(49),grune(49),expvol(49)
c
        data bolts/6.33870e-06/
        data idebye/0,0,344,1440,0,0,0,0,0
     >  ,0,158,400,428,0,0,0,0,0,91
     >  ,230,360,420,380,630,410,467,445,450,343
     >  ,327,320,0,0,0,0,0,56,147,280
     >  ,291,275,450,350,600,480,274,225,209,108/
        data grune/0.,0.,1.18,1.18,0.,0.,0.,0.,0.
     >  ,0.,1.31,1.48,2.19,0.,0.,0.,0.,0.,1.37
     >  ,1.16,1.17,1.18,1.05,1.30,2.07,1.66,1.93,1.88,2.
     >  ,2.01,2.,0.,0.,0.,0.,0.,1.67,1.,.89
     >  ,.83,1.58,1.6,2.6,3.2,2.23,2.28,2.36,2.23,2.37/
        data expvol/0.,143.7,0.,54.54,0.,0.,0.,0.,0.
     >  ,0.,254.5,151.4,109.9,0.,0.,0.,0.,0.,481.3
     >  ,291.1,168.7,120.3,93.48,80.63,82.84,78.95,74.72,73.42,78.92
     >  ,99.35,132.4,0.,0.,0.,0.,0.,598.9,373.6,194.7
     >  ,139.9,119.2,102.7,97.25,92.54,93.70,99.67,111.9,142.9,179.2/
c
        factor=2**icrys
        vol=.125*factor*alat**3
c
        do 100 is=1,nsites
        iz=z(is)+.1
        do ns=1,mspn
c
           if (iz < 1 .or. iz > 49) then
              ezero=0.0d0
              tpvzer=0.0d0
           else
              if(expvol(iz).ne.0.0) then
                 ezero=1.125*bolts*idebye(iz)*
     >                 (expvol(iz)/vol)**grune(iz)
              else
                 ezero=0.0d0
              endif
              tpvzer =3.*grune(iz)*ezero
           endif
c
           ezpt(is,ns)=ezero
           tpzpt(is,ns)=tpvzer
c
        enddo
  100   continue
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine keepem(icrys,alat,conc,z,toten,tpv,rmsrh,nsites,mspn)
c =====================================================================
      implicit real*8 (a-h,o-z)
c
c ---------------------------------------------------------------------
      parameter (ipcomp=4)
c ---------------------------------------------------------------------
c
      character*8  hd1,hd2
      character heding*51,lines(100)*80
c
      dimension z(ipcomp),conc(ipcomp),iz(ipcomp)
c
      data hd1,hd2/' at. # =',' conc  ='/
      data heding/
     >'   #              energy(ryd)  press(mbar)   rmsrho'/
c _____________________________________________________________________
c
       factor=2**icrys
       vol=.125*factor*alat**3
       p=147.1*tpv/(3.*vol)
       eold=0.
       pold=0.
       do 1 i=1,ipcomp
         iz(i)=0
         if(i.le.nsites)iz(i)=z(i)+.01
  1    continue
c
       number=0
c
        rewind 20
        read(20,'(////16x,f13.5,f12.5)',err=200,end=200) eold,pold
c
        rewind 20
  101   i=number+1
        read(20,'(a80)',err=200,end=200) lines(i)
        number=number+1
        if(number.lt.100) go to 101
  200   continue
        num=max(number-3,1)
c
        write(6,'(20x,''output from subroutine keepem:'')')
        write(6,'(1x,78(''_''))')
        write(6,'('' icrys='',i3,'' alat='',f8.5,
     >   '' ncmp='',i3,'' nspn='',i3)') icrys,alat,nsites,mspn
        write(6,'(a8,i5,9i7)') hd1,(iz(i),i=1,nsites)
        write(6,'(a8,10f7.4)') hd2,(conc(i),i=1,nsites)
        write(6,'(a50)') heding
        write(6,'(1x,i3,12x,f13.5,f12.5,2x,e9.2)') num,toten,p,rmsrh
        write(6,'(1x,78(''_''))')
c
        rewind 20
        write(20,'('' icrys='',i3,'' alat='',f8.5,
     >   '' ncmp='',i3,'' nspn='',i3)') icrys,alat,nsites,mspn
        write(20,'(a8,10i7)') hd1,(iz(i),i=1,nsites)
        write(20,'(a8,10f7.4)') hd2,(conc(i),i=1,nsites)
        write(20,'(a50)') heding
c
        write(20,'(1x,i3,12x,f13.5,f12.5,2x,e9.2)') num,toten,p,rmsrh
        if (number.ge.5) write(20,'(a80)') (lines(i),i=5,number)
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
         difmax=0.000001
c
         dife=abs(toten-eold)
         difp=abs(p-pold)*3*vol/147.1
c
         if(dife.ge.difmax) return
         if(difp.ge.difmax) return
c
         write(6,*) ' converged total energy and pressure'
         stop
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine ensum(vchd,evsum,efest,zv,elec,conc,
     >    ncomp,nspins,lmax,r,v,mesh,nsize,jmt,jws,rws,vdif,
     >    ndiv,ebott,eimag,ipr)
c     = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
c
      implicit real*8(a-h,o-z)
      parameter (ipmesh=641,ipcomp=4,ipspin=2,ipepts=1949,iplmax=3)
c
      dimension r(ipmesh),v(ipmesh,ipcomp,ipspin)
      dimension evsum(ipcomp,ipspin),zsemi(ipcomp,ipspin)
      dimension vchd(ipmesh,ipcomp,ipspin),zv(ipcomp,ipspin)
      dimension elec(ipcomp,ipspin),efest(ipspin),conc(ipcomp)
c
      complex*16 chden(ipmesh,ipcomp)
      complex*16 engrid(ipepts),chgmt(ipcomp),chgws(ipcomp)
      complex*16 esolv,de,sqrtm1,delt(iplmax+1),tmat(iplmax+1)
      complex*16 deltl,cta(iplmax+1)
c
      parameter (sqrtm1=(0.,1.))
c
c .....................................................................
c
c       subroutine to sum over the energies (complex plane)
c                  forms valence charge density:  vchd
c                  calculated one-electron sum:  evsum
c                  updates fermi energy (stored in efest)
c       works in a fixed moment mode, fermi energy for
c             each spin may be different
c
c ......................................................................
              ndelt=lmax+1
              isws=0
              iphsft=0
              if(ndiv.lt.0) then
                   ndiv=abs(ndiv)
                   iphsft=1
                   isws=1
              endif
c             write(6,*) ' ndelt=',ndelt
c             write(6,*) ' jws=',jws
c
c
       do 700 is=1,nspins
c .................................  do 700   loop over spins
c
c
         itries=0
   66    continue
         itries=itries+1
c
         zvtot=0.0
         do 100 ic=1,ncomp
             zvtot=zvtot+conc(ic)*zv(ic,is)
             evsum(ic,is)=0.0
             elec(ic,is)=0.0
             do 75 i=1,mesh
             vchd(i,ic,is)=0.0
   75        continue
  100    continue
c -----------------------------------------------------------------------
c
             if(zvtot.le.1.e-12) then
                   if(abs(zvtot).lt.1.e-8) go to 700
                   write(6,*)' stopping in ensum: number of valence',
     >                       ' electrons is negative'
                   stop
              endif
c
c -----------------------------------------------------------------------
          call egrid(ndiv,ebott,efest(is),eimag,engrid,nepts,ipr)
c -----------------------------------------------------------------------
c
          do 400 iept=1,nepts
c ..................................  do 400  loop over energies
c
             esolv=engrid(iept)
             de=.5*(engrid(min(iept+1,nepts))-engrid(max(iept-1,1)))
c
             do 300 ic=1,ncomp
c               write(6,*)' entering solvit: e=',esolv
c -----------------------------------------------------------------------
                call solvit(esolv,r,v(1,ic,is),mesh,nsize,jmt,jws,rws,
     >          lmax,chden(1,ic),cta,delt,tmat,chgmt(ic),chgws(ic),isws)
c -----------------------------------------------------------------------
                if(iphsft.eq.1) then
                  do 105 ll=1,ndelt
                  deltl=delt(ll)
                  l=ll-1
                  write(6,'('' cmp'',i2,'' spn'',i2,''  e(ry)='',
     >               2f7.4,3x,''l='',i1,2x,''ph.shift='',2e13.6)')
     >                ic,is,esolv,l,deltl
  105             continue
                  go to 300
                endif
                do 200 i=1,jws
                 vchd(i,ic,is)=vchd(i,ic,is)+de*chden(i,ic)
  200           continue
                evsum(ic,is)=evsum(ic,is)+de*esolv*chgws(ic)
                elec(ic,is)=elec(ic,is)+de*chgws(ic)
  300        continue
  400     continue
c
           if(iphsft.eq.1) go to 700
c
           sum=0.0
           dosef=0.0
           do 500 ic=1,ncomp
             sum=sum+conc(ic)*elec(ic,is)
             dosef=dosef+conc(ic)*chgws(ic)
  500      continue
c
           difn=zvtot-sum
           dife=difn/dosef
           if(ipr.ge.1) then
              write(6,'(1x,78(''-''))')
              write(6,*)'   total val.charge=',zvtot
              write(6,*)'        integrals  =',sum
              write(6,*)'   d.o.s.  at  ef  =',dosef
              write(6,*)' adjustment in ef  =',dife
           endif
c
           sumn=0.0
           do 600 ic=1,ncomp
             do 550 i=1,jws
               vchd(i,ic,is)=vchd(i,ic,is)+dife*chden(i,ic)
  550        continue
             evsum(ic,is)=evsum(ic,is)+esolv*dife*chgws(ic)
             elec(ic,is)=elec(ic,is)+dife*chgws(ic)
             sumn=sumn+conc(ic)*elec(ic,is)
c            write(6,*) ' spin=',is,'   comp=',ic,'   zv=',elec(ic,is)
  600     continue
c
c          write(6,*) is,sumn,zvtot
c
           write(6,'(1x,78(''-''))')
           dfemx=.1*abs(efest(is))
           if(abs(dife).gt.dfemx)  then
                  imulp=1
                  if(abs(dife).gt.2*dfemx)imulp=2
                  dife=imulp*sign(dfemx/itries,dife)
                  efold=efest(is)
                  efest(is)=efold+dife
                  write(6,*)' estimate of ef is too crude, old ef='
     >                      ,efold,'   new try=',efest(is)
                 go to 66
            endif
            efest(is)=efest(is)+dife
            write(6,*)' new fermi energy is ',efest(is)
c
  700  continue
c                     if phase shifts are the only things required: stop
       if(iphsft.eq.1) stop
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine egrid(ndiv,ebott,etop,eimag,el,nepts,ipr)
c =====================================================================
c
      implicit real*8(a-h,o-z)
c
      complex*16 el,sqrtm1
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      parameter (ipepts=1949)
      parameter (tenth=0.10,sqrtm1=(0.,1.))
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      dimension el(ipepts)
c
c     ***************************************************************
c     constructs rectangular energy contour
c
c          elow is the bottom of the contour (on the real axis)
c          efest is the top of the contour (on the real axis)
c          eimag is the shift of the real axis
c          dele is the energy increment  =   .1 / ngrid
c          el is the grid
c          nepts is the total number of points on this grid
c          ndiv is the number of divisions in .1 rydberg
c
c     ****************************************************************
c
c
      elow=ebott
      efest=etop
      dele=tenth/ndiv
      divs=(efest-elow)/dele
      numr=divs + 1
      elow=efest-dele*numr
      divi=eimag/dele
      numi=divi
      isum= 2 + numr + 2*numi
c
      if(isum.gt.ipepts) go to 500
      el(1)= elow
      it=1
c
      if(numi.lt.1) go to 11
      do 10 n=1,numi
      it=it+1
      el(it)=el(it-1)+dele*sqrtm1
  10  continue
  11  continue
c
      do 20 n=1,numr
      it=it+1
      el(it)=el(it-1)+dele
  20  continue
c
      if(numi.lt.1) go to 31
      do 30 n=1,numi
      it=it+1
      el(it)=el(it-1)-sqrtm1*dele
  30  continue
  31  continue
c
      nepts=it
      if(ipr.gt.1) write(6,1111) (n,el(n),n=1,it)
 1111 format(1x,i5,2e20.5)
      return
c
 500  continue
      write(6,1122)
 1122 format('  too many points on energy grid: stop')
       stop
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine solvit(e,r,v,mesh1,nsize1,
     >        jmt,jws,rws,lmax,chden,cta,delt,tmat,chgmt,chgws,isws)
c =====================================================================
c
c solves the schrodinger equation at arbitrary energy (complex)
c solves for both the regular and irregular solutions
c sums over angular momentum numbers: zero through ndelt-1
c  and returns the charge density contribution at this energy.
c
c     needs:   e    the energy (complex) in rydbergs
c              v    r*v in rydbergs
c            mesh1   the dimension of the r-grid
c            ndelt  the number of angular momentum numbers
c
c  returns:   chden  the charge density, summed over angular momentum
c              cta   the cotangent of the phase shifts
c              delt  the phase shifts
c              tmat  the t matrix = -1./q  sin(delta)*exp(i*delta)
c
c ---------------------------------------------------------------------
           implicit complex*16 (a-h,o-z)
           complex*16 ii
           real*8 r,rws,v,pi,fac,b0,del,delold,ylag,one
c
           parameter (ipmesh=641,ipspin=2,iplmax=3)
           parameter (one=1.,sqrtm1=(0.,1.))
           parameter (pi=3.14159265358979)
c
           dimension bjl(iplmax+1,ipmesh),bnl(iplmax+1,ipmesh)
           dimension c1l(ipmesh),s1l(ipmesh)
           dimension c2l(ipmesh),s2l(ipmesh)
           dimension cd(ipmesh),sd(ipmesh)
           dimension y(ipmesh),ii(ipmesh)
           dimension r(ipmesh),v(ipmesh)
           dimension chden(ipmesh)
           dimension zzr(ipmesh),zjr(ipmesh),tmp(ipmesh)
           dimension zlzl(iplmax+1),zljl(iplmax+1)
           dimension zlzlws(iplmax+1),zljlws(iplmax+1)
           dimension cta(iplmax+1),delt(iplmax+1),tmat(iplmax+1),fac(11)
           dimension bj(iplmax+1),bh(iplmax+1),bn(iplmax+1)
           dimension dj(iplmax+1),dh(iplmax+1),dn(iplmax+1)
c
           data fac/ 1.0, 3.0, 15.0, 105.0, 945.0, 10395.0, 135135.0,
     >              2027025.0, 34459425.0, 654729075.0, 13749310575.0/
c _____________________________________________________________________
c
c                            statement function
c
       catang(gam) = (.0,.5)*log(((0.,1.)+gam)/((0.,1.)-gam))
c _____________________________________________________________________
c
c
        q=sqrt(e)
c
        ndelt=lmax+1
        mesh=mesh1
        nsize=nsize1
        nblock=mesh/nsize
        do 76 ir=1,mesh
            chden(ir)=(0.0,0.0)
   76   continue
c
        do 78 ir=2,mesh
          zz=q*r(ir)
          call ricbes(ndelt,zz,bj,bh,bn,dj,dh,dn)
          do lp=1,ndelt
            bjl(lp,ir)=bj(lp)
            bnl(lp,ir)=bn(lp)
          enddo
   78   continue
c
        b0=v(1)
c _____________________________________________________________________
c                                                           main loop
c                                                         =============
      do 500 loop=0,lmax
        l=loop
        m=l+1
c ---------------------------------------------------------------------
        b1=v(1)
        f1=(v(2)-v(1))/r(2)
        f2=(v(3)-v(1))/r(3)
        f3=(v(4)-v(1))/r(4)
        bx1=f1
        bx2=(f2-f1)/(r(3)-r(2))
        bx3=(f3-f1)/(r(4)-r(2))
        bx3=(bx3-bx2)/(r(4)-r(3))
        b4=bx3
        b3=bx2-bx3*(r(2)+r(3))
        b2=bx1-bx2*r(2)+bx3*r(2)*r(3)
        a=q**(l+1) / fac(l+1)
        b=b1*a/(2*l+2)
        c=(a*(b2-e)+b*b1)/(4*l+6)
        d=(a*b3+b*(b2-e)+c*b1)/(6*l+12)
        f=(a*b4+b*b3+c*(b2-e)+d*b1)/(8*l+20)
c ---------------------------------------------------------------------
c
        do 10 i=2,5
           zz=q*r(i)
           x=r(i)**l
           p1=a+r(i)*( b+ r(i)*( c+ r(i)*( d+ f*r(i) )))
           p2=b+r(i)*(2.*c+r(i)*(3.*d+4.*f*r(i) ))
           y(i)=x*r(i)*p1
           dyi=(l+1)*x*p1+x*r(i)*p2
c
           if(l.eq.0) then
              bjp=cos(zz)
              bnp=sin(zz)
           else
              bjp=(l*bjl(m,i)-zz*bjl(l,i))/zz
              bnp=(l*bnl(m,i)-zz*bnl(l,i))/zz
           endif
           c1l(i)=y(i)*bnp-bnl(m,i)*dyi/q
           s1l(i)=y(i)*bjp-bjl(m,i)*dyi/q
           vy=y(i)*v(i)/zz
           cd(i)=-bnl(m,i)*vy
           sd(i)=-bjl(m,i)*vy
   10   continue
c ---------------------------------------------------------------------
        delold=r(6)-r(5)
        ityp=0
        do 188 n=6,jmt
           zz=q*r(n)
           vl=v(n)/zz
           del=r(n)-r(n-1)
c ---------------------------------------------------------------------
c
           if(abs(del-delold).gt.1.e-10)  then
               ityp=1
               delold=del
               m4=n-7
               m3=n-5
               m2=n-3
               m1=n-1
c
            else
c
               if(ityp.eq.0) then
                  m4=n-4
                  m3=n-3
                  m2=n-2
                  m1=n-1
                endif
                if(ityp.eq.2) then
                  ityp=0
                  m4=n-5
                  m3=n-3
                  m2=n-2
                  m1=n-1
                endif
                if(ityp.eq.1) then
                  ityp=2
                  m4=n-6
                  m3=n-4
                  m2=n-2
                  m1=n-1
                endif
            endif
c ---------------------------------------------------------------------
c
            h24=del/24.0
            c1l(n)=c1l(m1)+h24*(55.0*cd(m1)-59.0*cd(m2)
     >                         +37.0*cd(m3)-9.0*cd(m4))
            s1l(n)=s1l(m1)+h24*(55.0*sd(m1)-59.0*sd(m2)
     >                      +37.0*sd(m3)-9.0*sd(m4))
            do 187 ijk=1,5
              y(n)=c1l(n)*bjl(m,n)-s1l(n)*bnl(m,n)
              cd(n)=-bnl(m,n)*y(n)*vl
              sd(n)=-bjl(m,n)*y(n)*vl
              c1l(n)=c1l(m1)+h24*(9.0*cd(n)+19.0*cd(m1)
     >                              -5.0*cd(m2)+cd(m3))
              s1l(n)=s1l(m1)+h24*(9.0*sd(n)+19.0*sd(m1)
     >                              -5.0*sd(m2)+sd(m3))
  187       continue
          y(n)=c1l(n)*bjl(m,n)-s1l(n)*bnl(m,n)
c ---------------------------------------------------------------------
  188     continue
c
c ---------------------------------------------------------------------
          do 245 i=jmt+1,mesh
            c1l(i)=c1l(jmt)
            s1l(i)=s1l(jmt)
            c2l(i)=(0.0,0.0)
            s2l(i)=(1.0,0.0)
            y(i)=c1l(i)*bjl(m,i)-s1l(i)*bnl(m,i)
                if(isws.ne.1) ii(i)=bjl(m,i)
  245     continue
c ---------------------------------------------------------------------
c
c                   calculate phase shifts
c ---------------------------------------------------------------------
          csl=c1l(jmt)
          snl=s1l(jmt)
          cta(m)=csl/snl
          tinv=q*( (0.0,1.0)-cta(m) )
          tmat(m)=1./tinv
          tandl=snl/csl
          delt(m)=catang(tandl)
c _____________________________________________________________________
          if(isws.eq.1) go to 500
c
c  find the irregular solution    <===================================
c  find the irregular solution    <===================================
c  find the irregular solution    <===================================
c  find the irregular solution    <===================================
c  find the irregular solution    <===================================
c
         nlim=3
         do 255 nn=0,nlim
            n=jmt-nn
            zz=q*r(n)
            ii(n)=bjl(m,n)
            vi=ii(n)*v(n)/zz
            cd(n)=-bnl(m,n)*vi
            sd(n)=-bjl(m,n)*vi
            c2l(n)=(1.0,0.0)
            s2l(n)=(0.0,0.0)
  255   continue
c ---------------------------------------------------------------------
        do 258 nn=1,nlim
            n=jmt+nn
            zz=q*r(n)
            vl=ylag(r(n),r,v,0,3,jmt,iex)/zz
            m4=n-4
            m3=n-3
            m2=n-2
            m1=n-1
            del=r(n)-r(n-1)
            h24=del/24.0
            c2l(n)=c2l(m1)+h24*(55.0*cd(m1)-59.0*cd(m2)
     >                      +37.0*cd(m3)-9.0*cd(m4))
            s2l(n)=s2l(m1)+h24*(55.0*sd(m1)-59.0*sd(m2)
     >                      +37.0*sd(m3)-9.0*sd(m4))
            do 257 ijk=1,5
               ii(n)=c2l(n)*bjl(m,n)-s2l(n)*bnl(m,n)
               cd(n)=-bnl(m,n)*ii(n)*vl
               sd(n)=-bjl(m,n)*ii(n)*vl
               c2l(n)=c2l(m1)+h24*(9.0*cd(n)+19.0*cd(m1)
     >                             -5.0*cd(m2)+cd(m3))
               s2l(n)=s2l(m1)+h24*(9.0*sd(n)+19.0*sd(m1)
     >                             -5.0*sd(m2)+sd(m3))
  257        continue
         ii(n)=bjl(m,n)
c ---------------------------------------------------------------------
  258    continue
c _____________________________________________________________________
c
         delold=r(n)-r(n-1)
         do 265 nn=1,jmt-2
            n=jmt-nn
            del=r(n+1)-r(n)
c
            m1=n+1
            m2=n+2
            m3=n+3
            m4=n+4
c ---------------------------------------------------------------------
           if(abs(del-delold).gt.1.e-10) then
                      a1=cd(m1)
                      a2=cd(m2)
                      a3=cd(m3)
                      call qfit(delold,a1,a2,a3,f1,f2)
                      cd(m2)=f1
                      cd(m3)=a2
                      cd(m4)=f2
c
                      a1=sd(m1)
                      a2=sd(m2)
                      a3=sd(m3)
                      call qfit(delold,a1,a2,a3,f1,f2)
                      sd(m2)=f1
                      sd(m3)=a2
                      sd(m4)=f2
                      delold=del
           endif
c ---------------------------------------------------------------------
           h24=-del/24.0
           zz=q*r(n)
           vl=v(n)/zz
           c2l(n)=c2l(m1)+h24*(55.0*cd(m1)-59.0*cd(m2)
     >                        +37.0*cd(m3)-9.0*cd(m4))
           s2l(n)=s2l(m1)+h24*(55.0*sd(m1)-59.0*sd(m2)
     >                        +37.0*sd(m3)-9.0*sd(m4))
           do 263 ijk=1,5
              ii(n)=c2l(n)*bjl(m,n)-s2l(n)*bnl(m,n)
              cd(n)=-bnl(m,n)*ii(n)*vl
              sd(n)=-bjl(m,n)*ii(n)*vl
              c2l(n)=c2l(m1)+h24*(9.0*cd(n)+19.0*cd(m1)
     >                               -5.0*cd(m2)+cd(m3))
              s2l(n)=s2l(m1)+h24*(9.0*sd(n)+19.0*sd(m1)
     >                               -5.0*sd(m2)+sd(m3))
  263     continue
c
          ii(n)=c2l(n)*bjl(m,n)-s2l(n)*bnl(m,n)
  265   continue
c --------------------------------------------------------------------
c            <==== found the irregular solution ====>
c            <==== found the irregular solution ====>
c            <==== found the irregular solution ====>
c            <==== found the irregular solution ====>
c            <==== found the irregular solution ====>
c --------------------------------------------------------------------
c
c  now check wronskian
c
c                          \        /
c                           \  /\  /
c                            \/  \/
c
c          do 267 n=2,jmt,20
c          w=(c2l(n)*s1l(n)-c1l(n)*s2l(n))/s1l(jmt)
c          write(6,*) ' n,r,w=',n,r(n),w
c 267      continue
c --------------------------------------------------------------------
c
c     chden = -1./pi * 2l+1  * 4 pi r**2  ( zzt-zj )
c
c         the 4*pi is cancelled by angular integration
c         the faulkner and stocks definition of z and j
c         include the spherical harmonic ylm
c
c --------------------------------------------------------------------
        do 275 n=2,mesh
          y(n)=-y(n)/snl
          zzr(n)=y(n)**2
          ii(n)=ii(n)/q
          zjr(n)=y(n)*ii(n)
          xxy = ((2*l+1)/pi)*(zzr(n)*tmat(m)-zjr(n))
          chden(n)=chden(n) +  sqrtm1*xxy
c         if(mod(n, 40).eq.1) write(6,'(3x,''at r('',i3,'')='',f9.5,
c    >               5x,''contrib to ch.den='',2e16.9)')
c    >                 n,r(n),-xxy
  275   continue
c _____________________________________________________________________
        call intquc(zzr,tmp,r,nblock,nsize)
        zlzl(m)=-tmp(jmt)/pi
        zlzlws(m)=-ylagc(rws,r,tmp,0,3,mesh,iex)/pi
        call intquc(zjr,tmp,r,nblock,nsize)
        zljl(m)=tmp(jmt)/pi
        zljlws(m)=ylagc(rws,r,tmp,0,3,mesh,iex)/pi
c       write(6,'(''   m'',t10,''zlzl(mt)'',t44,''zlzl(ws)''/
c    >          1x,i3,2e16.8,2x,2e16.8)')
c    >       m,(2*l+1)*zlzl(m),(2*l+1)*zlzlws(m)
c       write(6,'(''   m'',t10,''zljl(mt)'',t44,''zljl(ws)''/
c    >          1x,i3,2e16.8,2x,2e16.8)')
c    >       m,(2*l+1)*zljl(m),(2*l+1)*zljlws(m)
c _____________________________________________________________________
        zztmzj=zlzl(m)*tmat(m)+zljl(m)
        zztzjw=zlzlws(m)*tmat(m)+zljlws(m)
c       write(6,'('' l='',i1,3x,''within  r='',f8.5,'' mt radius''
c    >     ,3x,''-(zzt-zj)/pi=('',e12.5,'','',e12.5,'')''/
c    >            '' l='',i1,3x,''within  r='',f8.5,'' ws radius''
c    >     ,3x,''-(zzt-zj)/pi=('',e12.5,'','',e12.5,'')'') ')
c    >             l,r(jmt),zztmzj,l,rws,zztzjw
c _____________________________________________________________________
c                                                  end of main loop
  500   continue
c _____________________________________________________________________
        if(isws.eq.1) return
        call intquc(chden,tmp,r,nblock,nsize)
        chgmt=tmp(jmt)
        chgws=ylagc(rws,r,tmp,0,3,mesh,iex)
c _____________________________________________________________________
c              write(6,*)' lmax=',lmax,'   energy=',e
c       write(6,'(7x,''within  r='',f8.5,'' mt radius''
c    >     ,3x,''charge density'',e12.5,'','',e12.5/
c    >            7x,''within  r='',f8.5,'' ws radius''
c    >     ,3x,''charge density'',e12.5,'','',e12.5) ')
c    >             r(jmt),chgmt,rws,chgws
c _____________________________________________________________________
        return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine ricbes (nn,y,bj,bh,bn,dj,dh,dn)
c =====================================================================
c calculates the riccati-bessel functions jl, hl and nl
c this version calculates the derivatives
c ---------------------------------------------------------------------
      implicit complex*16 (a-h,o-z)
      real*8 flm,big
      complex*16 bj(nn),bh(nn),bn(nn)
      complex*16 dj(nn),dh(nn),dn(nn)
      complex*16 one,sqrtm1,zero,x,y,aj
      parameter (zero=(0.,0.),one=(1.,0.), sqrtm1=(0.,1.),big=1.e6)
c _____________________________________________________________________
c
      x=y
      if (abs(x).eq.0.0) go to 200
      bj(1) =  sin(x)
      bn(1) = -cos(x)
      bj(2) = bj(1)/x + bn(1)
      bn(2) = bn(1)/x - bj(1)
      dj(1) = -bn(1)
      dn(1) =  bj(1)
      dj(2) = -(bj(1)/x - dj(1))/x + dj(1)
      dn(2) = -(bn(1)/x - dn(1))/x - dn(1)
c ---------------------------------------------------------------------
c  recursion relations
c ---------------------------------------------------------------------
      bh(1)=-sqrtm1*exp( sqrtm1*x )
      bh(2)=( one/x - sqrtm1 )*bh(1)
      dh(1)= sqrtm1*bh(1)
      dh(2)= ((one - sqrtm1*x)*dh(1) - bh(1)/x ) /x
      do 730 l=3,nn
        flm = l+l-3
        bj(l) = (flm*bj(l-1)/x-bj(l-2))
        bn(l) = (flm*bn(l-1)/x-bn(l-2))
        bh(l) = (flm*bh(l-1)/x-bh(l-2))
        dj(l) = (flm*(dj(l-1)-bj(l-1)/x)/x) - dj(l-2)
        dn(l) = (flm*(dn(l-1)-bn(l-1)/x)/x) - dn(l-2)
        dh(l) = (flm*(dh(l-1)-bh(l-1)/x)/x) - dh(l-2)
  730   continue
c ---------------------------------------------------------------------
      if( abs(x) .le. 0.90 ) then
c
c  power-series for j if abs(x) is smaller than 0.9
c
          do 700 l=3,nn
             ll=l
             bj(l) = zero
             aj = one
             k = 0
             do 500 i=1,l
             aj = aj*x/(i+i-1)
  500        continue
  600        bj(l) = bj(l)+aj
             k = k+1
             aj = -aj*x**2/((2*k)*(2*(l+k)-1))
             if (abs(big*aj).gt.abs(bj(l))) goto 600
             bh(l)=bj(l)+sqrtm1*bn(l)
  700     continue
c
      endif
      return
c _____________________________________________________________________
c
  200  write(6,99) x
   99  format(' trouble in bessel: argument= ',2f10.5)
        stop
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        subroutine findns(r,mesh,jmt,jws,rws,ic,is,vr,zed,
     >          numbs,npbs,lbs,ebs,em,deepe,thresh,iprt)
c =====================================================================
       implicit real*8 (a-h,o-z)
c
       parameter (ipmesh=641,nsize=40)
       parameter (ipspin=2,ipcomp=4,ipbst=20)
c
       dimension r(ipmesh),vr(ipmesh)
       dimension      numbs(ipcomp,ipspin)
       dimension lbs(ipbst,ipcomp,ipspin)
       dimension npbs(ipbst,ipcomp,ipspin)
       dimension ebs(ipbst,ipcomp,ipspin)
c _____________________________________________________________________
c
c   routine finds the total number of bounds states
c     for all l  and below em=emin and above deepe=edeep
c     useful when starting a calculation
c _____________________________________________________________________
c
            emin=em
            edeep=deepe
            nstate=0
            do 97 l=0,3
            lb=l
c .....................................................................
            nst2=-1
            call nbound(zed,r,vr,mesh,jmt,lb,emin,nst2,ajl)
            if(edeep.ge.emin) then
                nst1=0
            else
                nst1=-1
                call nbound(zed,r,vr,mesh,jmt,lb,edeep,nst1,ajl)
            endif
c .....................................................................
            if(iprt.ge.1)
     >      write(6,*)' l=',l,'  number of bound states=',nst2-nst1
            do 91 nst=nst1+1,nst2
               nstate=nstate+1
               lbs(nstate,ic,is)=l
               ebs(nstate,ic,is)=0.0
               npbs(nstate,ic,is)=l+nst
  91         continue
c
  97     continue
c
c .....................................................................
c
         write(6,'(1x,78(''-''))')
         write(6,*) '        below e=',emin,' (ryd)'
         if(edeep.lt.emin)
     >      write(6,*) '    and above e=',edeep,' (ryd)'
         write(6,*) ' total number of bound states=',nstate
         numbs(ic,is)=nstate
         if( nstate.gt.0 .and. iprt.gt.1 )  then
            write(6,'(/4x,''n'',4x,''l'',2x,''energy(ryd)''
     >                    /3x,27(''-'') )')
            do 99 i=1,nstate
              npbs1=npbs(i,ic,is)
              lbs1=lbs(i,ic,is)
              ebs1=ebs(i,ic,is)
              write(6,'(2i5,f13.4)') npbs1,lbs1,ebs1
   99       continue
         endif
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        subroutine nbound (zed,r,v,mesh1,jmt,lq,emn,nq,ajl)
c =====================================================================
        implicit real*8 (a-h,o-z)
        parameter (ipmesh=641,onepd=1.0001,eminup=-.0005,esrch=-.5)
        dimension r(ipmesh),v(ipmesh),vme(ipmesh)
        dimension y(ipmesh),dy(ipmesh),aj(ipmesh)
c _____________________________________________________________________
c
c          this routine finds the number of bound states below the
c              energy  emn  with the angular momentum l=lq
c
c _____________________________________________________________________
       l=lq
       mesh=mesh1
       emin=min(emn,eminup)
       en=emin
c ---------------------------------------------------------------------
c zero out arrays
        do 3 i=1,mesh1
        y(i)=0.0
        dy(i)=0.0
        aj(i)=0.0
    3   continue
c ---------------------------------------------------------------------
        vme(1)=0.0
        do 4 j=2,jmt
        vme(j)=(v(j)+(l*l+l)/r(j))/r(j)   - emin
    4   continue
        do 5 j=jmt+1,mesh
        vme(j)=         (l*l+l)/r(j)**2   - emin
    5   continue
        lll=mesh
c ---------------------------------------------------------------------
        call outsch(zed,emin,l,lll,vme,v,r,y,dy)
c ---------------------------------------------------------------------
c
c       the radial solution can be expressed in terms of the bessel and
c          hankel function,  the coefficient multiplying the bessel
c          function is aj, the number of zeroes of aj(r) counts the
c          bound states below the current energy, when aj at the
c          muffin-tin radius is zero the current energy corresponds
c          to a state which matches onto a hankel function outside
c          the muffin-tin, i.e. a bound state energy
c
c ---------------------------------------------------------------------
       if(nq.lt.0) then
           num=0
           if(emin.lt.esrch) then
                 do 6 i=3,jmt
                 if(y(i).eq.0.0) then
                      num=num+1
                      go to 6
                 else
                      is=sign(onepd, y(i) )
                      js=sign(onepd,y(i-1))
                      if(is.ne.js) num=num+1
                 endif
   6             continue
                 ajl=hcoef(emin,r(jmt),l,y(jmt),dy(jmt))
c
           else
c
                 aj(2)=hcoef(emin,r(2),l,y(2),dy(2))
                 do 7 i=3,jmt
                 aj(i)=hcoef(emin,r(i),l,y(i),dy(i))
                 if(aj(i).eq.0.0) then
                     num=num+1
                     go to 7
                 else
                     if(aj(i-1).eq.0.0) go to 7
                     is=sign(onepd, aj(i) )
                     js=sign(onepd,aj(i-1))
                     if(is.ne.js) num=num+1
                 endif
   7             continue
                 ajl=aj(jmt)
c
           endif
           nq=num
c
c
        else
c          write(6,*)' number of nodes at ',emin,'   is ',num
          ajl=hcoef(emin,r(jmt),l,y(jmt),dy(jmt))
        endif
c ---------------------------------------------------------------------
c number of nodes  n-l-1
c
c       write(6,*)' found ',num,' zeroes of j-coef'
        return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        subroutine bstate(r,mesh,jmt,jws,rws,ic,is,vr,rhoc,zed,
     >          numbs,npbs,lbs,ebs,occbs,em,thresh,qsum,iprt)
c =====================================================================
       implicit real*8 (a-h,o-z)
       character*7 abad
       parameter (ipmesh=641,nsize=40)
       parameter (ipspin=2,ipcomp=4,ipbst=20)
       parameter (eswdp=-25.)
       dimension r(ipmesh),vr(ipmesh),rhoc(ipmesh)
       dimension temp(ipmesh),rhoc1(ipmesh)
       dimension      numbs(ipcomp,ipspin)
       dimension npbs(ipbst,ipcomp,ipspin),ebs(ipbst,ipcomp,ipspin),
     >           lbs(ipbst,ipcomp,ipspin),occbs(ipbst,ipcomp,ipspin)
c _____________________________________________________________________
c
       do 91 k=1,ipmesh
       rhoc(k)=0.0
  91   continue
c ...................................
         nmbs=numbs(ic,is)
         do  97 ist=1,nmbs
           en = ebs(ist,ic,is)
           nb =npbs(ist,ic,is)
           lb = lbs(ist,ic,is)
           if(en.ge.0) en=-zed**2/nb**4
           if(en.lt.eswdp) then
c .....................................................................
              call deepst(abad,zed,r,vr,mesh,jmt,en,thresh,nb,lb,
     >                 rhoc1,kkk,em,temp,qi)
c .....................................................................
           else
c .....................................................................
              call shalst(abad,zed,r,vr,mesh,jmt,en,.01*thresh,nb,lb,
     >                 rhoc1,kkk,em,temp,qi)
c .....................................................................
           endif
c
           if(iprt.ge.1) write(6,'(1x,'' bound state: '',a7,3x,
     >            ''energy level at '',e16.9,'' ry'' )') abad,en
c
           twolp1=2*lb+1
           qsum=qsum-twolp1*qi
           occbs(ist,ic,is)=twolp1
           ebs(ist,ic,is)=en
           do 95 k=1,kkk
             rhoc(k)=rhoc(k)+rhoc1(k)*twolp1
  95       continue
  97     continue
         nblock=mesh/nsize
         call intqud(rhoc,rhoc1,r,nblock,nsize)
         chgws=ylag(rws,r,rhoc1,0,3,mesh,iex)
          write(6,'(''  charge:'',f11.8,5x,''within mt r:'',f11.8
     >       ,5x,''within ws r:'',f11.8)')
     >         rhoc1(mesh),rhoc1(jmt),chgws
c .....................................................................
         if(iprt.lt.2) return
         mpr=40
         if(iprt.lt.3) mpr=20
               write(6,'(1x,78(''=''))')
               write(6,'(1x,i3,''ch='',f8.5,2x,i3,''ch='',f8.5,
     >    2x,i3,''ch='',f8.5,2x,i3,''ch='',f8.5,2x,i3,''ch='',f8.5)')
     >             (k,rhoc(k),k=1,mesh,mpr)
               write(6,'(1x,78(''=''))')
       return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        subroutine deepst
     > (abad,zed,r,v,mesh1,jmt,en,thresh,nq,lq,chden,kkk,emn,qi2,qii)
c =====================================================================
        implicit real*8 (a-h,o-z)
        character*7 abad,fndt,tble
        parameter (ipmesh=641,onepd=1.0001)
        parameter (fndt='foundit',tble='trouble')
        dimension r(ipmesh),v(ipmesh),chden(ipmesh)
        dimension temp(ipmesh),qi2(ipmesh),vme(ipmesh)
        dimension y(ipmesh),dy(ipmesh),yin(ipmesh),dyin(ipmesh)
        data eminup/-0.0005/
c _____________________________________________________________________
c
       abad=tble
       l=lq
       n=nq
       mesh=mesh1
       if(en.ge.0.0) en=-zed**2/n**4
       emin=min(emn,eminup)
       energy=en
       en=emin
       eold=99.99
c ---------------------------------------------------------------------
c zero out arrays
        do 3 i=1,mesh1
        qi2(i)=0.0
        y(i)=0.0
        dy(i)=0.0
        chden(i)=0.0
    3   continue
c ---------------------------------------------------------------------
        bigkr=150.
        if(abs(emin).gt.400.) bigkr=125.
        if(abs(emin).gt.1000.) bigkr=100.
        bigkr=.5*bigkr
        vme(1)=0.0
        do 4 j=2,jmt
        vme(j)=(v(j)+(l*l+l)/r(j))/r(j)   - emin
    4   continue
        do 5 j=jmt+1,mesh
        vme(j)=         (l*l+l)/r(j)**2   - emin
    5   continue
        do 6 i=2,jmt
        j=i
        if(r(i)*sqrt(-emin) .gt. bigkr)go to 7
   6    continue
        j=mesh
   7    lll=j
c ---------------------------------------------------------------------
        call outsch(zed,emin,l,lll,vme,v,r,y,dy)
c ---------------------------------------------------------------------
        num=0
        do 8 i=3,lll
        if(y(i).eq.0.0) then
                 num=num+1
                 go to 8
        else
             if(y(i-1).eq.0) go to 8
             is=sign(onepd, y(i) )
             js=sign(onepd,y(i-1))
             if(is.ne.js) num=num+1
        endif
   8    continue
c ---------------------------------------------------------------------
c number of nodes  n-l-1
c
c            write(6,*)' found ',num,' nodes'
       nnode=n-l-1
       if(nnode.ge.num) then
               write(6,1740) emin,nq,lq
 1740          format('  no bound state below energy=',f8.4,
     >         ', with quantum numbers: n=',i2,', l=',i2 )
               return
       endif
c ---------------------------------------------------------------------
c
c  initialize parameters for search
c
       maxit=40
       maxjt=20
       iext=0
       iter=0
       jter=0
c
       v(1)=-2*zed
c ---------------------------------------------------------------------
c ....................  main loop .....................................
c ---------------------------------------------------------------------
  10   continue
       iter=iter+1
       if(iter.gt.maxit) go to 800
       if(energy.gt.0) energy=eold/2.
       bigkr=150.
       if(abs(energy).gt.400.) bigkr=125.
       if(abs(energy).gt.1000.) bigkr=100.
       bigkr=.5*bigkr
c
c find classical turning point ........................................
c
        do 40 j=2,jmt
        vme(j)=(v(j)+(l*l+l)/r(j))/r(j)   - energy
  40    continue
        do 41 j=jmt+1,mesh
        vme(j)=          (l*l+l)/r(j)**2  - energy
  41    continue
c ---------------------------------------------------------------------
c
        do 42 k=3,jmt
            j=jmt+3-k
            if(vme(j).le.0.0) go to 45
   42   continue
c ---------------------------------------------------------------------
        j=jmt-1
c ---------------------------------------------------------------------
  45    match=j+1
c ---------------------------------------------------------------------
c
        match=min0(match,jmt)
        do 60 i=match,mesh
        j=i
        if(r(i)*sqrt(-energy) .gt. bigkr ) go to 61
  60    continue
  61    last=j
        if(iext.eq.0) lll=last
        if(iext.ne.0) lll=match+1
        mmmm=match-1
c ---------------------------------------------------------------------
c    outward integration of the schrodinger equation
c ---------------------------------------------------------------------
        call outsch(zed,energy,l,lll,vme,v,r,y,dy)
c ---------------------------------------------------------------------
        gamout=dy(match)/(sqrt(-energy)*y(match))
c                                                 count number of nodes
        num=0
        do 80 i=3,lll
        if(y(i).eq.0.) then
             num=num+1
             go to 80
        else
             if(y(i-1).eq.0.0) go to 80
             is=sign(onepd, y(i) )
             js=sign(onepd,y(i-1))
             if(is.ne.js) num=num+1
        endif
  80    continue
c ---------------------------------------------------------------------
c            write(6,*)'  iter=',iter,'   energy=',energy
c            write(6,*)' found ',num,' nodes'
c ---------------------------------------------------------------------
        if(iter.eq.1) then
             iext=0
             if(num.eq.nnode+1) then
                  enext=1.15*energy
                  jter=0
                  go to 87
             endif
             if(num.eq.nnode) then
                  enext=.85*energy
                  jter=0
                  go to 87
             endif
             jter=jter+1
             if(jter.gt.maxjt)  then
                  write(6,*) ' cannot find right number of nodes '
                  return
             endif
c
             iter=0
             if(num.lt.nnode) then
                  energy=.5*energy
                  go to 10
             else
                  energy=1.5*energy
                  go to 10
             endif
         endif
c ---------------------------------------------------------------------
         if(iter.eq.2)  then
             if(num.eq.nold) then
                  iter=1
                  jter=jter+1
c
                  if(jter.gt.maxjt) then
                       write(6,*) ' trouble in second iteration'
                       return
                  endif
c
                  dir=+1.
                  if(num.eq.nnode) dir=-1.
                  enext=energy-.8*dir*abs(energy-eold)
                  if(enext.ge.0.0) enext=.5*energy
                  energy=enext
                  go to 10
             endif
             if(iabs(num-nold).ne.1) then
                  iter=1
                  jter=jter+1
                  if(jter.gt.maxjt) then
                       write(6,*) ' trouble in second iteration'
                       return
                  endif
                  enext=.5*energy+.5*eold
                  energy=enext
                  go to 10
             endif
             go to 81
        endif
c ---------------------------------------------------------------------
         if(abs(eold-energy).lt.thresh) go to 200
c ---------------------------------------------------------------------
         if(num.ne.nnode) go to 83
   81    dir=+1
         go to 85
   83    dir=-1.
   85    enext=energy+.5*dir*abs(energy-eold)
   87    jter=0
   88    nold=num
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       if(abs(1.-eold/energy) .lt. 0.1 .or. iext.eq.1) then
c ---------------------------------------------------------------------
c   inward integration of schrodinger equation
c ---------------------------------------------------------------------
       if(match.lt.jmt) then
         call inwsch(zed,energy,mmmm,last,jmt,mesh,l,vme,r,v,yin,dyin)
         gamin=dyin(match)/(sqrt(-energy)*yin(match))
c ---------------------------------------------------------------------
       else
         write(6,*)' ++++++++++ calling gamhnk +++++++++++++'
         gamin=gamhnk(energy,r(jmt),l)
       endif
c ---------------------------------------------------------------------
             rat=(1.+gamin)/(1.-gamin)
             rat=rat-((1.+gamout)/(1.-gamout))
c
c            write(6,*) 'gamin,gamout  ',gamin,gamout
c            write(6,*) ' ratin=',rat
c            write(6,*)' ratdf=',rat
c
             if(rat.eq.0) go to 200
c
             if(iext.eq.1) then
                 if( abs(1.-ratold/rat).lt.1.e-30) then
                    enext=.5*(eold+energy)
                 else
                    enext=energy+(eold-energy)/( 1. - ratold/rat)
                 endif
             endif
c
             iext=1
             ratold=rat
         endif
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
         eold=energy
         energy=enext
         go to 10
c .....................................................................
c
  800    continue
c ---------------------------------------------------------------------
c ..........  can not find state  .....................................
c ---------------------------------------------------------------------
         write(6,*) ' too many tries for n=',nq,' l=',l
         return
c
  200    continue
c ---------------------------------------------------------------------
c .................... found state ............................
c ---------------------------------------------------------------------
       if(match.lt.jmt) then
         call inwsch(zed,energy,mmmm,last,jmt,mesh,l,vme,r,v,yin,dyin)
c ---------------------------------------------------------------------
         gamin=dyin(match)/(sqrt(-energy)*yin(match))
c ---------------------------------------------------------------------
       else
         write(6,*) 'calling hank,gamhnk'
         gamin=gamhnk(energy,r(jmt),l)
         sqe=sqrt(-energy)
         do 210 i=match,last
          yin(i)=hank(energy,r(i),l,dh)
          dyin(i)=dh*sqe
  210    continue
       endif
c ---------------------------------------------------------------------
c
         abad=fndt
         kkk=last
         do 225 i=2,match
         chden(i)=y(i)**2
  225    continue
         anorm=y(match)/yin(match)
         do 226 i=match,last
         y(i)=anorm*yin(i)
         dy(i)=anorm*dyin(i)
         chden(i)=y(i)**2
  226    continue
         do 224 i=last+1,mesh
         chden(i)=0.0
         dy(i)=0.0
         y(i)=0.0
  224    continue
c
         en=energy
c ---------------------------------------------------------------------
         nsize=40
         nblock=mesh1/nsize
         call intqud(chden,temp,r,nblock,nsize)
         qi2(1)=0.0
         do 227 i=2,last
              chden(i)=chden(i)/temp(mesh1)
              qi2(i)=(dy(i)-(l+1)*y(i)/r(i))**2/temp(mesh1)
  227    continue
         do 228 i=last+1,mesh1
              chden(i)=0.0
              qi2(i)=0.0
  228    continue
         call intqud(qi2,temp,r,nblock,nsize)
         qii=temp(jmt)
c       mpr=20
c       write(6,2234) (i,chden(i),i=1,mesh1,mpr)
c2234   format( 4(i4,'.den=',f9.5,2x) )
         return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        subroutine shalst
     > (abad,zed,r,v,mesh1,jmt,en,thresh,nq,lq,chden,kkk,emn,qi2,qii)
c =====================================================================
        implicit real*8 (a-h,o-z)
        character*7 abad,fndt,tble
c
        parameter (ipmesh=641)
        parameter (fndt='foundit',tble='trouble')
c
        dimension r(ipmesh),v(ipmesh),chden(ipmesh)
        dimension temp(ipmesh),qi2(ipmesh),vme(ipmesh)
        dimension y(ipmesh),dy(ipmesh),yin(ipmesh),dyin(ipmesh)
c
        data eminup/-0.0005/,onepd/1.000001/
c _____________________________________________________________________
c
       abad=tble
       l=lq
       n=nq
       mesh=mesh1
       if(en.ge.0.0) en=-zed**2/n**4
       emin=min(emn,eminup)
       energy=en
       en=emin
       eold=99.99
c ---------------------------------------------------------------------
c
c zero out arrays
        do 3 i=1,mesh1
        qi2(i)=0.0
        y(i)=0.0
        dy(i)=0.0
        chden(i)=0.0
    3   continue
c ---------------------------------------------------------------------
       num =-1
       call nbound(zed,r,v,mesh,jmt,l,emin,num,ajl)
c
       nnode=n-l-1
       if(nnode.ge.num) then
               write(6,1740) emin,nq,lq
 1740          format('  no bound state below energy=',f8.4,
     >         ', with quantum numbers: n=',i2,', l=',i2 )
               return
       endif
c ---------------------------------------------------------------------
c
c  initialize parameters for search
c
       maxit=100
       maxjt=20
       iext=0
       iter=0
       jter=0
c
       v(1)=-2*zed
c ---------------------------------------------------------------------
c ....................  main loop .....................................
c ---------------------------------------------------------------------
  10   continue
c
       iter=iter+1
c       write(6,*)' iteration= ',iter,'   energy= ',energy
       if(iter.gt.maxit) go to 800
       if(energy.gt.0) energy=eold/2.
       num =-1
c ---------------------------------------------------------------------
       call nbound(zed,r,v,mesh,jmt,l,energy,num,ajl)
c ---------------------------------------------------------------------
        if(iter.eq.1) then
             iext=0
               nod1=num
               e1=energy
ccc            iasgn1=sign(onepd,apl)
c
             if(num.eq.nnode+1) then
                  enext=1.15*energy
                  jter=0
                  go to 87
             endif
c
             if(num.eq.nnode) then
                  enext=.85*energy
                  jter=0
                  go to 87
             endif
             jter=jter+1
             if(jter.gt.maxjt)  then
                  write(6,*) ' cannot find right number of nodes '
                  return
             endif
             iter=0
            if(num.lt.nnode) then
                enext=.5*energy
                go to 10
            else
                enext=1.5*energy
                go to 10
            endif
c
         endif
c ---------------------------------------------------------------------
         if(iter.eq.2)  then
             nod2=num
             e2=energy
ccc          iasgn2=sign(onepd,apl)
             if(num.eq.nold) then
                  iter=1
                  jter=jter+1
                  if(jter.gt.maxjt) then
                       write(6,*) ' trouble in second iteration'
                       return
                  endif
                  enext=energy-1.05*ajl*(energy-eold)/(ajl-ajold)
                  if(jter.lt.9) enext=2*energy-eold
                  if(enext.ge.0.0) enext=.5*energy
                      ajold=ajl
                      eold=energy
                  energy=enext
                  go to 10
             endif
             if(iabs(num-nold).ne.1) then
                  iter=1
                  jter=jter+1
                  if(jter.gt.maxjt) then
                       write(6,*) ' trouble in second iteration'
                       return
                  endif
                  enext=.5*energy+.5*eold
                  energy=enext
                  go to 10
             endif
             go to 81
        endif
c ---------------------------------------------------------------------
         if(abs(eold-energy).lt.thresh) go to 200
c ---------------------------------------------------------------------
   81    continue
c .....................................................................
            enext=energy-1.1*ajl*(energy-eold)/(ajl-ajold)
            if(iter.ge.5) then
                   onp=1. + .1/iter
                   enext=energy-onp*ajl*(energy-eold)/(ajl-ajold)
                   if(iter.ge.8) then
                       edf=(enext-energy)*(enext-eold)
                       if(edf.ge.0) enext=.5*(eold+energy)
                   endif
             endif
             ndf=(num-nod1)*(num-nod2)
             if(ndf.ne.0) then
                   write(6,*) ' node trouble: num=',num
                   write(6,*) '              nod1=',nod1,
     >                                   '     nod2=',nod2
             endif
c .....................................................................
   87    jter=0
   88    nold=num
         ajold=ajl
         eold=energy
         energy=enext
         go to 10
c .....................................................................
c
  800    continue
c .....................................................................
c ..........  can not find state  .....................................
c ..........  can not find state  .....................................
c ..........  can not find state  .....................................
c ..........  can not find state  .....................................
c .....................................................................
         write(6,*) ' too many tries for n=',nq,' l=',l
         return
c
  200    continue
c .....................................................................
c ....................     found state     ............................
c ....................     found state     ............................
c ....................     found state     ............................
c ....................     found state     ............................
c .....................................................................
        do 205 j=2,jmt
        vme(j)=(v(j)+(l*l+l)/r(j))/r(j) - energy
  205   continue
        sqe=sqrt(-energy)
        call outsch(zed,energy,l,jmt,vme,v,r,y,dy)
c ---------------------------------------------------------------------
        do 210 i=jmt,mesh
        yin(i)=hank(energy,r(i),l,dh)
        dyin(i)=sqe*dh
  210   continue
c
c ---------------------------------------------------------------------
c          tmp=gamhnk(energy,r(jmt),l)
c          gamout=dy(jmt)/(sqrt(-energy)*y(jmt))
c          gamin=dyin(jmt)/(sqrt(-energy)*yin(jmt))
c           write(6,*)'lgder:gamhnk,hank,outsch ',
c    >        tmp,gamin,gamout
c ---------------------------------------------------------------------
c
         abad=fndt
         kkk=mesh
         do 225 i=2,jmt
         chden(i)=y(i)**2
  225    continue
         anorm=y(jmt)/yin(jmt)
         do 226 i=jmt,mesh
            y(i)=anorm*yin(i)
            dy(i)=anorm*dyin(i)
            chden(i)=y(i)**2
  226    continue
c
         en=energy
c ---------------------------------------------------------------------
         nsize=40
         nblock=mesh1/nsize
         call intqud(chden,temp,r,nblock,nsize)
         qi2(1)=0.0
         do 227 i=2,mesh
            chden(i)=chden(i)/temp(mesh1)
            qi2(i)=(dy(i)-(l+1)*y(i)/r(i))**2/temp(mesh1)
  227    continue
         call intqud(qi2,temp,r,nblock,nsize)
         qii=temp(jmt)
c ---------------------------------------------------------------------
c        mpr=20
c        write(6,2234) (i,chden(i),i=1,mesh1,mpr)
c2234    format( 4(i4,'.den=',f9.5,2x) )
c ---------------------------------------------------------------------
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        subroutine outsch(zed,energy,l,match,vme,v,r,y,dy)
c =====================================================================
        implicit real*8 (a-h,o-z)
        parameter (ipmesh=641)
        dimension vme(ipmesh),v(ipmesh),r(ipmesh)
        dimension y(ipmesh),dy(ipmesh)
        dimension fac(5)
        data fac/1.,3.,15.,105.,945./
c _____________________________________________________________________
c
        if(l.lt.0 .or. l.gt.4) then
             write(6,*) ' angular momentum quatum number is ',l
             write(6,*) ' stopping in outsch: l out of bounds'
             stop
        endif
c ---------------------------------------------------------------------
        last=match
        e=energy
        b1=v(1)
        f1=(v(2)-v(1))/r(2)
        f2=(v(3)-v(1))/r(3)
        f3=(v(4)-v(1))/r(4)
        bx1=f1
        bx2=(f2-f1)/(r(3)-r(2))
        bx3=(f3-f1)/(r(4)-r(2))
        bx3=(bx3-bx2)/(r(4)-r(3))
        b4=bx3
        b3=bx2-bx3*(r(2)+r(3))
        b2=bx1-bx2*r(2)+bx3*r(2)*r(3)
        a=sqrt(-e)**(l+1) / fac(l+1)
        b=b1*a/(2*l+2)
        c=(a*(b2-e)+b*b1)/(4*l+6)
        d=(a*b3+b*(b2-e)+c*b1)/(6*l+12)
        f=(a*b4+b*b3+c*(b2-e)+d*b1)/(8*l+20)
c
        do 10 i=2,5
        x=r(i)**l
        p1=a+r(i)*( b+ r(i)*( c+ r(i)*( d+ f*r(i) )))
        p2=b+r(i)*(2.*c+r(i)*(3.*d+4.*f*r(i) ))
        y(i)=x*r(i)*p1
        dy(i)=(l+1)*x*p1+x*r(i)*p2
   10   continue
c ---------------------------------------------------------------------
      delold=r(6)-r(5)
      ityp=0
      do 40 n=6,last
         del=r(n)-r(n-1)
c ---------------------------------------------------------------------
         if(abs(del-delold).gt.1.e-10)  then
                 ityp=1
                 delold=del
                 m4=n-7
                 m3=n-5
                 m2=n-3
                 m1=n-1
         else
c
                 if(ityp.eq.0) then
                     m4=n-4
                     m3=n-3
                     m2=n-2
                     m1=n-1
                 endif
                 if(ityp.eq.2) then
                     ityp=0
                     m4=n-5
                     m3=n-3
                     m2=n-2
                     m1=n-1
                 endif
                 if(ityp.eq.1) then
                     ityp=2
                     m4=n-6
                     m3=n-4
                     m2=n-2
                     m1=n-1
                 endif
         endif
c ---------------------------------------------------------------------
        h24=del/24.
c
        dy(n)=dy(m1)+h24*(55.*vme(m1)*y(m1) -59.*vme(m2)*y(m2)
     >                  + 37.*vme(m3)*y(m3) - 9.*vme(m4)*y(m4))
        do 37 ijk=1,5
        y(n)=y(m1)+h24*(   9.*dy(n)  +19.*dy(m1) -5.*dy(m2) + dy(m3))
        dy(n)=dy(m1)+h24*( 9.*vme(n)*y(n) +19.*vme(m1)*y(m1)
     >                   - 5.*vme(m2)*y(m2) +  vme(m3)*y(m3))
   37   continue
        y(n)=y(m1)+h24*(  9.*dy(n) +19.*dy(m1) -5.*dy(m2) + dy(m3))
c ---------------------------------------------------------------------=
   40  continue
c
             return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine inwsch(zed,energy,matc1,last1,jmt,mesh,l,vme,r,v,y,dy)
c ======================================================================
        implicit real*8 (a-h,o-z)
        parameter (ipmesh=641)
        dimension r(ipmesh),v(ipmesh),vme(ipmesh),y(ipmesh),dy(ipmesh)
c _____________________________________________________________________
c
        match=matc1
        last=last1
        if(mod(last-2,40).lt.4) last=min0(last+6,mesh)
c
        n=last
        p=sqrt(-energy)
c ---------------------------------------------------------------------=
        do 10 i=1,4
        xpr=exp(-p*r(n))
        y(n)=xpr
        dy(n)=-p*xpr
            if(i.eq.1) then
                   ym3=y(n)
                   dym3=dy(n)
                   vmem3=vme(n)
            endif
            if(i.eq.2) then
                   ym2=y(n)
                   dym2=dy(n)
                   vmem2=vme(n)
            endif
            if(i.eq.3) then
                   ym1=y(n)
                   dym1=dy(n)
                   vmem1=vme(n)
            endif
        n=n-1
  10    continue
c ---------------------------------------------------------------------
        nend=n
        n=n+1
        delold=r(n+1)-r(n)
        do 80 nn=match,nend
        n=n-1
        del=r(n+1)-r(n)
        m1=n+1
        m2=n+2
        m3=n+3
        m4=n+4
c ---------------------------------------------------------------------
        if(abs(del-delold).gt.1.e-10) then
             a1=dy(m1)
             a2=dy(m2)
             a3=dy(m3)
c
             call qrfit(delold,a1,a2,a3,f1,f2)
c
             dym1=a1
             dym2=f1
             dym3=a2
             dym4=f2
c
             a1=y(m1)
             a2=y(m2)
             a3=y(m3)
c
             call qrfit(delold,a1,a2,a3,f1,f2)
c
             ym1=a1
             ym2=f1
             ym3=a2
             ym4=f2
c
             a1=vme(m1)
             a2=vme(m2)
             a3=vme(m3)
c
             call qrfit(delold,a1,a2,a3,f1,f2)
c
             vmem1=a1
             vmem2=f1
             vmem3=a2
             vmem4=f2
c
             delold=del
c ---------------------------------------------------------------------
      else
c ---------------------------------------------------------------------
             ym4=ym3
             ym3=ym2
             ym2=ym1
             ym1=y(n+1)
c
             dym4=dym3
             dym3=dym2
             dym2=dym1
             dym1=dy(n+1)
c
             vmem4=vmem3
             vmem3=vmem2
             vmem2=vmem1
             vmem1=vme(n+1)
      endif
c ---------------------------------------------------------------------
c
      h24=-del/24.
      dy(n)=dym1 + h24 * (55.*vmem1*ym1-59.*vmem2*ym2
     >                   +37.*vmem3*ym3- 9.*vmem4*ym4)
      do 70 ijk=1,5
      y(n)=ym1 + h24*(9.*dy(n)+19.*dym1-5.*dym2+dym3)
      dy(n)=dym1 + h24 * ( 9.*vme(n)*y(n)+19.*vmem1*ym1
     >                    -5.*vmem2*ym2+vmem3*ym3)
   70 continue
      y(n)=ym1 + h24*(9.*dy(n)+19.*dym1-5.*dym2+dym3)
c ---------------------------------------------------------------------
   80 continue
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine qfit(del,a1,a2,a3,f1,f2)
c =====================================================================
       implicit complex*16 (a-h,o-z)
       real*8 del,del2,delh
c _____________________________________________________________________
c
c  given a function evaluated at three equally spaced points
c     returns values at the two midpoints using
c           quadratic interpolation
c                needs     del,  the spacing
c                       a1,a2,a3, the 3 function evaluations
c              returns    f1,f2,  the interpolated values
c        f(x)= a + bx + cx(x-d)
c        f(0)=a1
c        f(d)=a2
c        f(2d)=a3
c        f1=f(d/2)
c        f2=f(3d/2)
c _____________________________________________________________________
c
          del2=2*del
          delh=del*.5
          a=a1
          b=(a2-a1)/del
          c=(((a3-a1)/del2)-b)/del
          f1=a+delh*(b-c*delh)
          f2=a+3.*delh*(b+c*delh)
          return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine qrfit(del,a1,a2,a3,f1,f2)
c =====================================================================
       implicit real*8 (a-h,o-z)
       real*8 del,del2,delh
c _____________________________________________________________________
c
c  given a function evaluated at three equally spaced points
c     returns values at the two midpoints using
c         quadratic interpolation
c               needs     del,  the spacing
c                       a1,a2,a3, the 3 function evaluations
c              returns    f1,f2,  the interpolated values
c        f(x)= a + bx + cx(x-d)
c        f(0)=a1
c        f(d)=a2
c        f(2d)=a3
c        f1=f(d/2)
c        f2=f(3d/2)
c ------------------------------------------------------------------==
          del2=2*del
          delh=del*.5
          a=a1
          b=(a2-a1)/del
          c=(((a3-a1)/del2)-b)/del
          f1=a+delh*(b-c*delh)
          f2=a+3.*delh*(b+c*delh)
          return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
          function hank(e,r,l,dh)
c =====================================================================
c
c     hankel function   ......  expects e to be real and negative
c
c .....................................................................
          implicit real*8 (a-h,o-z)
c
          x=sqrt(-e)*r
          xi=1./x
          xp=exp(-x)
c
          go to (100,101,102,103) , l+1
c
          write(6,*)' stopping in hank, l.gt.3'
          stop
c
c .....................   l=0
c
  100     hank=xp
          dh=-xp
          return
c
c .....................   l=1
c
  101     hank=(1.+xi)*xp
          dh=-hank-xi**2*xp
          return
c
c .....................   l=2
c
  102     hank=(1.+3*xi*(1.+xi))*xp
          dh=-hank-xi**2*(3.+6.*xi)*xp
          return
c
c .....................   l=3
c
  103     hank=(1.+xi*(6.+xi*(15.+3.*xi)))*xp
          dh=-hank-xi**2*(6+xi*(30+9*xi))*xp
          return
c
c .....................................................................
c
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
         function hcoef(e,r,l,y,dy)
c =====================================================================
         implicit real*8 (a-h,o-z)
c
         bh=hank(e,r,l,dh)
         hcoef=dh*y-bh*dy/sqrt(-e)
c
         return
      end
c ______________________________________________________________________
c //////////////////////////////////////////////////////////////////////
          function gamhnk(e,r,l)
c =====================================================================
          implicit real*8 (a-h,o-z)
c .....................................................................
c
c  logarithmic derivative of the hankel functions
c     gammma = dh/dx  /  h
c
c .....................................................................
c
          x=sqrt(-e)*r
          go to (100,101,102,103) , l+1
c
          write(6,*)' stopping in gamhnk, l.gt.3'
          stop
c
c .....................   l=0
c
  100     gamhnk=-1.
          return
c
c .....................   l=1
c
  101     gamhnk=-1. - 1./(x*x+x)
          return
c
c .....................   l=2
c
  102     gamhnk=-1. - 3.*(x+2.)/(x*(3.+x*(3.+x)))
          return
c
c .....................   l=3
c
  103     gamhnk=-1. - (9.+x*(30.+6.*x))/(x*(3.+x*(15.+x*(6.+x))))
          return
c
c .....................................................................
c
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
          function bessj(e,r,l,djl)
c =====================================================================
          implicit real*8 (a-h,o-z)
c
          dimension bj(4),dj(4)
c
          if(l.gt.3) then
            write(6,*)' stopping in bessj: l.gt.3'
            stop
           endif
c
          x=sqrt(-e)*r
c _____________________________________________________________________
          bj(1)=sin(x)
          dj(1)=cos(x)
          bj(2)=bj(1)/x-dj(1)
          dj(2)=bj(1)-bj(2)/x
c
          do 10 ll=3,4
             fl=2*ll-3
             bj(ll)= fl*bj(ll-1)/x-bj(ll-2)
             dj(ll)=-fl*(bj(ll-1)/x-dj(ll-1))/x-dj(ll-2)
  10      continue
c
         bessj=bj(l+1)
           djl=dj(l+1)
c
         return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine ptable(iz,symb,name)
c =====================================================================
        character*2  symb,sym(103)
        character*15 name,nam(103)
c _____________________________________________________________________
c
       sym(  1)='h '
       nam(  1)='hydrogen       '
       sym(  2)='he'
       nam(  2)='helium         '
       sym(  3)='li'
       nam(  3)='lithium        '
       sym(  4)='be'
       nam(  4)='berylium       '
       sym(  5)='b '
       nam(  5)='boron          '
       sym(  6)='c '
       nam(  6)='carbon         '
       sym(  7)='n '
       nam(  7)='nitrogen       '
       sym(  8)='o '
       nam(  8)='oxygen         '
       sym(  9)='f '
       nam(  9)='fluorine       '
       sym( 10)='ne'
       nam( 10)='neon           '
       sym( 11)='na'
       nam( 11)='sodium         '
       sym( 12)='mg'
       nam( 12)='mangnesium     '
       sym( 13)='al'
       nam( 13)='aluminum       '
       sym( 14)='si'
       nam( 14)='silicon        '
       sym( 15)='p '
       nam( 15)='phosphorus     '
       sym( 16)='s '
       nam( 16)='sulfur         '
       sym( 17)='cl'
       nam( 17)='chlorine       '
       sym( 18)='ar'
       nam( 18)='argon          '
       sym( 19)='k '
       nam( 19)='potassium      '
       sym( 20)='ca'
       nam( 20)='calcium        '
       sym( 21)='sc'
       nam( 21)='scandium       '
       sym( 22)='ti'
       nam( 22)='titanium       '
       sym( 23)='v '
       nam( 23)='vanadium       '
       sym( 24)='cr'
       nam( 24)='chromium       '
       sym( 25)='mn'
       nam( 25)='manganese      '
       sym( 26)='fe'
       nam( 26)='iron           '
       sym( 27)='co'
       nam( 27)='cobalt         '
       sym( 28)='ni'
       nam( 28)='nickel         '
       sym( 29)='cu'
       nam( 29)='copper         '
       sym( 30)='zn'
       nam( 30)='zinc           '
       sym( 31)='ga'
       nam( 31)='gallium        '
       sym( 32)='ge'
       nam( 32)='germanium      '
       sym( 33)='as'
       nam( 33)='arsenic        '
       sym( 34)='se'
       nam( 34)='selenium       '
       sym( 35)='br'
       nam( 35)='bromine        '
       sym( 36)='kr'
       nam( 36)='krypton        '
       sym( 37)='rb'
       nam( 37)='rubidium       '
       sym( 38)='sr'
       nam( 38)='strontium      '
       sym( 39)='y '
       nam( 39)='yttrium        '
       sym( 40)='zr'
       nam( 40)='zirconium      '
       sym( 41)='nb'
       nam( 41)='niobium        '
       sym( 42)='mo'
       nam( 42)='molybdenum     '
       sym( 43)='tc'
       nam( 43)='technetium     '
       sym( 44)='ru'
       nam( 44)='ruthenium      '
       sym( 45)='rh'
       nam( 45)='rhodium        '
       sym( 46)='pd'
       nam( 46)='palladium      '
       sym( 47)='ag'
       nam( 47)='silver         '
       sym( 48)='cd'
       nam( 48)='cadmium        '
       sym( 49)='in'
       nam( 49)='indium         '
       sym( 50)='sn'
       nam( 50)='tin            '
       sym( 51)='sb'
       nam( 51)='antimony       '
       sym( 52)='te'
       nam( 52)='tellurium      '
       sym( 53)='i '
       nam( 53)='iodine         '
       sym( 54)='xe'
       nam( 54)='xenon          '
       sym( 55)='cs'
       nam( 55)='cesium         '
       sym( 56)='ba'
       nam( 56)='barium         '
       sym( 57)='la'
       nam( 57)='lanthanum      '
       sym( 58)='ce'
       nam( 58)='cerium         '
       sym( 59)='pr'
       nam( 59)='praseodymium   '
       sym( 60)='nd'
       nam( 60)='neodymium      '
       sym( 61)='pm'
       nam( 61)='promethium     '
       sym( 62)='sm'
       nam( 62)='samarium       '
       sym( 63)='eu'
       nam( 63)='europium       '
       sym( 64)='gd'
       nam( 64)='gadolinium     '
       sym( 65)='tb'
       nam( 65)='terbium        '
       sym( 66)='dy'
       nam( 66)='dysprosium     '
       sym( 67)='ho'
       nam( 67)='holmium        '
       sym( 68)='er'
       nam( 68)='erbium         '
       sym( 69)='tm'
       nam( 69)='thulium        '
       sym( 70)='yb'
       nam( 70)='ytterbium      '
       sym( 71)='lu'
       nam( 71)='lutetium       '
       sym( 72)='hf'
       nam( 72)='hafnium        '
       sym( 73)='ta'
       nam( 73)='tantalum       '
       sym( 74)='w '
       nam( 74)='tungsten       '
       sym( 75)='re'
       nam( 75)='rhenium        '
       sym( 76)='os'
       nam( 76)='osmium         '
       sym( 77)='ir'
       nam( 77)='iridium        '
       sym( 78)='pt'
       nam( 78)='platinum       '
       sym( 79)='au'
       nam( 79)='gold           '
       sym( 80)='hg'
       nam( 80)='mercury        '
       sym( 81)='tl'
       nam( 81)='thallium       '
       sym( 82)='pb'
       nam( 82)='lead           '
       sym( 83)='bi'
       nam( 83)='bismuth        '
       sym( 84)='po'
       nam( 84)='polonium       '
       sym( 85)='at'
       nam( 85)='astatine       '
       sym( 86)='rn'
       nam( 86)='radon          '
       sym( 87)='fr'
       nam( 87)='francium       '
       sym( 88)='ra'
       nam( 88)='radium         '
       sym( 89)='ac'
       nam( 89)='actinium       '
       sym( 90)='th'
       nam( 90)='thorium        '
       sym( 91)='pa'
       nam( 91)='protactinium   '
       sym( 92)='u '
       nam( 92)='uranium        '
       sym( 93)='np'
       nam( 93)='neptunium      '
       sym( 94)='pu'
       nam( 94)='plutonium      '
       sym( 95)='am'
       nam( 95)='americium      '
       sym( 96)='cm'
       nam( 96)='curium         '
       sym( 97)='bk'
       nam( 97)='berkelium      '
       sym( 98)='cf'
       nam( 98)='californium    '
       sym( 99)='es'
       nam( 99)='einsteinium    '
       sym(100)='fm'
       nam(100)='fermium        '
       sym(101)='md'
       nam(101)='mendelevium    '
       sym(102)='no'
       nam(102)='nobelium       '
       sym(103)='lr'
       nam(103)='lawrencium     '
c _____________________________________________________________________
c
      if(iz.lt.0 .or. iz.gt.103) then
          write(6,*)' atomic number is incorrect, iz=',iz
          stop
      endif
c _____________________________________________________________________
c
      name=nam(iz)
      symb=sym(iz)
c _____________________________________________________________________
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine intqud(y,yi,x,nb,nin)
c =====================================================================
      implicit real*8 (a-h,o-z)
      parameter (ipmesh=641)
      dimension x(ipmesh),y(ipmesh),yi(ipmesh)
c _____________________________________________________________________
c
      nbig=nin*nb+1
      yi(1)=0.0
      do 10 i=2,nbig
          top=x(i)
          bot=x(i-1)
          j=max0(i-2,1)
          j=min0(j,nbig-3)
          y1=y(j)
          y2=y(j+1)
          y3=y(j+2)
          y4=y(j+3)
          x1=x(j)
          x2=x(j+1)
          x3=x(j+2)
          x4=x(j+3)
          a1=y1
          a2=(y2-y1)/(x2-x1)
          a3=(((y3-y1)/(x3-x1))-a2)/(x3-x2)
          a4=((((((y4-y1)/(x4-x1))-a2)/(x4-x2)))-a3)/(x4-x3)
          b1=a1-a2*x1+a3*x1*x2-a4*x1*x2*x3
          b2=a2-a3*(x1+x2)+a4*(x2*x3+x1*x2+x1*x3)
          b3=a3-a4*(x1+x2+x3)
          b4=a4
          xint=b1+(.5*b2+.25*b4*(top**2+bot**2))*(top+bot)
     >        +b3*(top**2+top*bot+bot**2)/3.
          yi(i)=yi(i-1)+xint*(top-bot)
   10  continue
c
        return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine smpold(y,yi,r,nb,nin)
c =====================================================================
      implicit real*8(a-h,o-z)
      parameter (ipmesh=641)
      dimension r(ipmesh),y(ipmesh),yi(ipmesh)
c _____________________________________________________________________
c
      bsum=0.
      asum=0.
      yi(1)=0.
      i=1
      kmax=nin/2
c     nin must be an even integer
      do 50 j=1,nb
      dx=r(i+1)-r(i)
      h=dx/3.
      g=h/4.
      do 40 k=1,kmax
      i=i+2
      bsum=bsum+h*(y(i-2)+4.*y(i-1)+y(i))
      yi(i)=asum+bsum
      yi(i-1)=yi(i-2)+g*(5.*y(i-2)+8.*y(i-1)-y(i))
  40  continue
      asum=yi(i)
      bsum=0.
  50  continue
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine intquc(y,yi,x,nb,nin)
c =====================================================================
      implicit complex*16 (a-h,o-w,y-z)
      implicit real*8 (x)
      real*8 top,bot
      parameter (ipmesh=641)
      dimension x(ipmesh),y(ipmesh),yi(ipmesh)
c _____________________________________________________________________
c
      nbig=nin*nb+1
      yi(1)=0.0
      do 10 i=2,nbig
          top=x(i)
          bot=x(i-1)
          j=max(i-2,1)
          j=min(j,nbig-3)
          y1=y(j)
          y2=y(j+1)
          y3=y(j+2)
          y4=y(j+3)
          x1=x(j)
          x2=x(j+1)
          x3=x(j+2)
          x4=x(j+3)
          a1=y1
          a2=(y2-y1)/(x2-x1)
          a3=(((y3-y1)/(x3-x1))-a2)/(x3-x2)
          a4=((((((y4-y1)/(x4-x1))-a2)/(x4-x2)))-a3)/(x4-x3)
          b1=a1-a2*x1+a3*x1*x2-a4*x1*x2*x3
          b2=a2-a3*(x1+x2)+a4*(x2*x3+x1*x2+x1*x3)
          b3=a3-a4*(x1+x2+x3)
          b4=a4
          cint=b1+(.5*b2+.25*b4*(top**2+bot**2))*(top+bot)
     >        +b3*(top**2+top*bot+bot**2)/3.
          yi(i)=yi(i-1)+cint*(top-bot)
   10  continue
c
        return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
         subroutine derv5(y,dy,x,nn)
c =====================================================================
         implicit real*8 (a-h,o-z)
         parameter (ipmesh=641)
         dimension x(ipmesh),y(ipmesh),dy(ipmesh)
c _____________________________________________________________________
c
c                                      lagrangian fit of y=f(x)
c                                      derivative is then calculated
c _____________________________________________________________________
         mesh=nn
         do 10 i=2,mesh
         i1=max0(i-2,2)
         i1=min0(i1,mesh-4)
         a=y(i1)
         b=(y(i1+1)-a)/(x(i1+1)-x(i1))
         c=(y(i1+2)-a)/(x(i1+2)-x(i1))
         c=(c-b)/(x(i1+2)-x(i1+1))
         d=(y(i1+3)-a)/(x(i1+3)-x(i1))
         d=(d-b)/(x(i1+3)-x(i1+1))
         d=(d-c)/(x(i1+3)-x(i1+2))
         e=(y(i1+4)-a)/(x(i1+4)-x(i1))
         e=(e-b)/(x(i1+4)-x(i1+1))
         e=(e-c)/(x(i1+4)-x(i1+2))
         e=(e-d)/(x(i1+4)-x(i1+3))
         dy(i)=b+c*( (x(i)-x(i1)) + (x(i)-x(i1+1)) )
     >        + d*(  (x(i)-x(i1))*(x(i)-x(i1+1))
     >             + (x(i)-x(i1+1))*(x(i)-x(i1+2))
     >             + (x(i)-x(i1+2))*(x(i)-x(i1))      )
     >      +e*( (x(i)-x(i1))*(x(i)-x(i1+1))*(x(i)-x(i1+2))
     >      +  (x(i)-x(i1+3))*(x(i)-x(i1+1))*(x(i)-x(i1+2))
     >      +  (x(i)-x(i1))*(x(i)-x(i1+1))*(x(i)-x(i1+3))
     >      +  (x(i)-x(i1))*(x(i)-x(i1+2))*(x(i)-x(i1+3)) )
   10    continue
c
         dy(1)=dy(2)
c _____________________________________________________________________
         return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
         subroutine derv4(y,dy,x,nn)
c =====================================================================
         implicit real*8 (a-h,o-z)
         parameter (ipmesh=641)
         dimension x(ipmesh),y(ipmesh),dy(ipmesh)
c _____________________________________________________________________
c
c                                      lagrangian fit of y=f(x)
c                                      derivitave is then calculated
c _____________________________________________________________________
         mesh=nn
         do 10 i=2,mesh
         i1=max0(i-1,2)
         i1=min0(i1,mesh-3)
         a=y(i1)
         b=(y(i1+1)-a)/(x(i1+1)-x(i1))
         c=(y(i1+2)-a)/(x(i1+2)-x(i1))
         c=(c-b)/(x(i1+2)-x(i1+1))
         d=(y(i1+3)-a)/(x(i1+3)-x(i1))
         d=(d-b)/(x(i1+3)-x(i1+1))
         d=(d-c)/(x(i1+3)-x(i1+2))
         dy(i)=b+c*( (x(i)-x(i1)) + (x(i)-x(i1+1)) )
     >        + d*(  (x(i)-x(i1))*(x(i)-x(i1+1))
     >             + (x(i)-x(i1+1))*(x(i)-x(i1+2))
     >             + (x(i)-x(i1+2))*(x(i)-x(i1))      )
   10    continue
c
         dy(1)=dy(2)
c
         return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
         subroutine derv3(y,dy,x,nn)
c =====================================================================
         implicit real*8 (a-h,o-z)
         parameter (ipmesh=641)
         dimension x(ipmesh),y(ipmesh),dy(ipmesh)
c _____________________________________________________________________
c
c                                      lagrangian fit of y=f(x)
c                                      derivative is then calculated
c _____________________________________________________________________
         mesh=nn
         do 10 i=2,mesh
         i1=max0(i-1,2)
         i1=min0(i1,mesh-2)
         a=y(i1)
         b=(y(i1+1)-a)/(x(i1+1)-x(i1))
         c=(y(i1+2)-a)/(x(i1+2)-x(i1))
         c=(c-b)/(x(i1+2)-x(i1+1))
         dy(i)=b+c*( (x(i)-x(i1)) + (x(i)-x(i1+1)) )
   10    continue
c
         dy(1)=dy(2)
         return
      end
c ______________________________________________________________________
c //////////////////////////////////////////////////////////////////////
      function ylag(xi,x,y,ind1,n1,imax,iex)
c =====================================================================
c
c  program authors  a. a. brooks and e.c. long,
c  computing technology center, union carbide corp., nuclear div.,
c  oak ridge, tenn.
c
c     lagrangian interpolation
c     xi is intepolated entry into x-array
c     n is the order of lagrangran interpolation
c     y is array from which ylag is obtained by interpolation
c     ind is the min-i for x(i).gt.xi
c     if ind=0,x-array will be searched
c     imax is max index of x-and y-arrays
c     extrapolation can occur,iex=-1 or +1
c
c
      implicit real*8 (a-h,o-z)
c
      dimension x(imax),y(imax)
c
      ind=ind1
      n=n1
      iex=0
      if (n.le.imax) go to 10
      n=imax
      iex=n
   10 if (ind.gt.0) go to 40
      do 20 j = 1,imax
         if(abs(xi-x(j)).lt.1.0e-10) go to 130
c        if (xi-x(j)) 30,130,20
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
   40 if (ind.gt.1) go to 50
      iex=-1
   50 inl=ind-(n+1)/2
      if (inl.gt.0) go to 60
      inl=1
   60 inu=inl+n-1
      if (inu.le.imax) go to 80
   70 inl=imax-n+1
      inu=imax
   80 s=0.
      p=1.
      do j=inl,inu
         p=p*(xi-x(j))
         d=1.
         do i=inl,inu
            if (i.ne.j) go to 90
            xd=xi
            go to 100
   90       xd=x(j)
  100       d=d*(xd-x(i))
         enddo    
         s=s+y(j)/d
      enddo
      ylag=s*p
      return
  130 ylag=y(j)
      return
      end
c ______________________________________________________________________
c //////////////////////////////////////////////////////////////////////
      function ylagc(xi,x,y,ind1,n1,imax,iex)
c =====================================================================
c
c  program authors  a. a. brooks and e.c. long,
c  computing technology center, union carbide corp., nuclear div.,
c  oak ridge, tenn.
c
c     lagrangian interpolation
c     xi is intepolated entry into x-array
c     n is the order of lagrangran interpolation
c     y is array from which ylag is obtained by interpolation
c     ind is the min-i for x(i).gt.xi
c     if ind=0,x-array will be searched
c     imax is max index of x-and y-arrays
c     extrapolation can occur,iex=-1 or +1
c
c
      implicit complex*16 (a-h,o-w,y-z)
      implicit real*8 (x)
      real*8 d,p
c
      dimension x(imax),y(imax)
c
      ind=ind1
      n=n1
      iex=0
      if (n.le.imax) go to 10
      n=imax
      iex=n
   10 if (ind.gt.0) go to 40
      do 20 j = 1,imax
         if(abs(xi-x(j)).lt.1.0e-10) go to 130
c        if (xi-x(j)) 30,130,20
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
   40 if (ind.gt.1) go to 50
      iex=-1
   50 inl=ind-(n+1)/2
      if (inl.gt.0) go to 60
      inl=1
   60 inu=inl+n-1
      if (inu.le.imax) go to 80
   70 inl=imax-n+1
      inu=imax
   80 s=(0.0,0.0)
      p=1.0
      do j=inl,inu
         p=p*(xi-x(j))
         d=1.
         do i=inl,inu
            if (i.ne.j) go to 90
            xd=xi
            go to 100
   90       xd=x(j)
  100       d=d*(xd-x(i))
         enddo
         s=s+y(j)/d
      enddo
      ylagc=s*p
      return
  130 ylagc=y(j)
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine dgamix(jiter,jws,ncomp,nspins,rhoins,rhoots,rhomix,
     >                         alpha)
c =====================================================================
        implicit real*8 (a-h,o-z)
        parameter (ipmesh=641,ipcomp=4,ipspin=2,ipits=4)
        parameter (ipsize=ipits+1,one=1.)
        dimension rhoins(ipits,ipmesh+1,ipcomp,ipspin)
        dimension rhoots(ipits,ipmesh+1,ipcomp,ipspin)
        dimension rhomix(ipmesh+1,ipcomp,ipspin)
        dimension a(ipsize,ipsize),b(ipsize)
c
       iter=min(jiter,ipits)
       idim=iter+1
       do 2 i=1,iter
         a(i,idim)=one
         a(idim,i)=one
         b(i)=0
   2   continue
       a(idim,idim)=0
       b(idim)=one
c
       do 5 jj=1,iter
         j=jj
         do ii=1,jj
         i=ii
         a(ii,jj)=dpdif(i,j,jws,ncomp,nspins,rhoins,rhoots)
         if(ii.ne.jj) a(jj,ii)=a(ii,jj)
         enddo
   5  continue
c
      call mxlneq(a,b,idim)
c
      if(nspins.eq.1) then
         is=1
         do 15 ic=1,ncomp
           do i=1,jws+1
           rhomix(i,ic,is)=sumand(rhoins(1,i,ic,is),rhoots(1,i,ic,is),
     >                                 b,iter,alpha)
           enddo
   15    continue
c
      else
           do 30 ic=1,ncomp
           do i=1,jws+1
              call summag(rhomix(i,ic,1),rhomix(i,ic,2),
     >                  rhoins(1,i,ic,1),rhoins(1,i,ic,2),
     >                  rhoots(1,i,ic,1),rhoots(1,i,ic,2),
     >                           b,iter,alpha)
           enddo
  30       continue
      endif
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        function dpdif(i,j,jws,ncomp,nspins,rhoins,rhoots)
c =====================================================================
        implicit real*8 (a-h,o-z)
        parameter (ipmesh=641,ipcomp=4,ipspin=2,ipits=4)
        dimension rhoins(ipits,ipmesh+1,ipcomp,ipspin)
        dimension rhoots(ipits,ipmesh+1,ipcomp,ipspin)
c
         sum=0
         do 10 is=1,nspins
         do ic=1,ncomp
         do 5 k=1,jws+1
         sum=sum+(rhoins(i,k,ic,is)-rhoots(i,k,ic,is))
     >          *(rhoins(j,k,ic,is)-rhoots(j,k,ic,is))
   5     continue
         enddo
  10   continue
c
         dpdif=sum
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        subroutine savit(iter,rho,jws,totz,ncomp,nspins,rhosav)
c =====================================================================
       implicit real*8 (a-h,o-z)
       parameter (ipmesh=641,ipcomp=4,ipspin=2,ipits=4)
       dimension rho(ipmesh,ipcomp,ipspin),totz(ipcomp,ipspin)
       dimension rhosav(ipits,ipmesh+1,ipcomp,ipspin)
c
       n=mod(iter-1,ipits)+1
       do 10 is=1,nspins
       do ic=1,ncomp
       rhosav(n,jws+1,ic,is)=totz(ic,is)
       do i=1,jws
       rhosav(n,i,ic,is)=rho(i,ic,is)
       enddo
       enddo
  10   continue
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       function sumand(cin,cot,b,niter,alpha)
c =====================================================================
       implicit real*8 (a-h,o-z)
       parameter (ipits=4)
       dimension cin(ipits),cot(ipits),b(ipits)
c
       iter=min(niter,ipits)
c
       optin=0
       optot=0
       do 10 i=1,iter
       optin=optin+b(i)*cin(i)
       optot=optot+b(i)*cot(i)
  10   continue
c
       sumand=optin+alpha*(optot-optin)
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine mxlneq(a,y,nn)
c =====================================================================
c
       implicit real*8 (a-h,o-z)
       parameter (ipits=4,ipsize=ipits+1,one=1.)
       dimension a(ipsize,ipsize),y(ipsize)
       n=nn
       do 10 i=2,n
            f1=-one/a(i-1,i-1)
            do  2 j=i,n
            f2=f1*a(j,i-1)
            do 5 k=1,n
            a(j,k)=a(j,k)+f2*a(i-1,k)
    5       continue
            y(j)=y(j)+f2*y(i-1)
    2       continue
   10  continue
c
c ..........................  calculate the determinant ...............
c
          det=one
          do 19 i=1,n
          det=det*a(i,i)
  19      continue
c
          write(6,1235) n,det
 1235     format(' in mxlneq: order, det',2x,i4,e17.9)
c
      y(n)=y(n)/a(n,n)
      do 20 ij=1,n-1
          i=n-ij
          do 15 j=1,ij
          y(i)=y(i)-y(i+j)*a(i,i+j)
   15     continue
      y(i)=y(i)/a(i,i)
   20 continue
c
c
       return
      end
c
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine summag(up,dn,cupin,cdnin,cupot,cdnot,b,niter,alpha)
c =====================================================================
      implicit real*8 (a-h,o-z)
      parameter (ipits=4)
      dimension b(ipits)
      dimension cupin(ipits),cdnin(ipits)
      dimension cupot(ipits),cdnot(ipits)
c
      iter=min(niter,ipits)
      ain=0
      aot=0
      xin=0
      xot=0
      do 10 i=1,iter
         aini=cupin(i)+cdnin(i)
         xini=cupin(i)-cdnin(i)
         aoti=cupot(i)+cdnot(i)
         xoti=cupot(i)-cdnot(i)
         ain=ain+b(i)*aini
         xin=xin+b(i)*xini
         aot=aot+b(i)*aoti
         xot=xot+b(i)*xoti
   10 continue
c
      c=ain+alpha*(aot-ain)
      up=.5*(c+xot)
      dn=.5*(c-xot)
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
c
c                                          end of program
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
c =====================================================================
