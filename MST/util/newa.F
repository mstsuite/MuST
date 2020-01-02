c
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
c  fjp version of the herman-skillman program       jan/feb 1987
c             spin-polarized     and      total energy
c /////////////////////////////////////////////////////////////////////
           implicit real*8 (a-h,o-z)
c _____________________________________________________________________
c
c                      order of subroutines
c
c        subroutine    .........   1  ....  redin
c        subroutine    .........   2  ....  xsetup
c        subroutine    .........   3  ....  config
c        subroutine    .........   4  ....  thofer
c        subroutine    .........   5  ....  intqud
c        subroutine    .........   6  ....  smpold
c        subroutine    .........   7  ....  potgen
c        subroutine    .........   8  ....  alpha2
c        subroutine    .........   9  ....  derv5
c        subroutine    .........  10  ....  bound
c        subroutine    .........  11  ....  outsch
c        subroutine    .........  12  ....  inwsch
c        subroutine    .........  13  ....  qrfit
c        subroutine    .........  14  ....  savout
c        subroutine    .........  15  ....  savin
c        subroutine    .........  16  ....  dpdif
c        subroutine    .........  17  ....  dgand
c        subroutine    .........  18  ....  mxlneq
c _____________________________________________________________________
c
           parameter (three=3.,third=1./three,twoth=2./three)
           character*2 sym
           character*15 name
           character*7 abad,fndt,tble
           dimension x(1001),r(1001),rho(1001,2),rv(1001,2),temp(1001)
           dimension q2(1001),tmp(1001)
           dimension occ(20,2),ncore(20,2),lcore(20,2),numcor(2),
     >                    elcnum(2),ecore(20,2),snlo(1001)
           dimension rvnew(1001,2),rhoold(1001,2)
           dimension occup(20),ncoup(20),lcoup(20),ecoup(20)
           dimension occdn(20),ncodn(20),lcodn(20),ecodn(20)
           common/e1/fndt,tble
           common/drvdr/drv(1001,2)
           common/excorr/excpo(1001,2),excen(1001)
c _____________________________________________________________________
c
               equivalence (occup,occ(1,1))
               equivalence (occdn,occ(1,2))
               equivalence (ncoup,ncore(1,1))
               equivalence (ncodn,ncore(1,2))
               equivalence (lcoup,lcore(1,1))
               equivalence (lcodn,lcore(1,2))
               equivalence (numcup,numcor(1))
               equivalence (numcdn,numcor(2))
               equivalence (ecoup,ecore(1,1))
               equivalence (ecodn,ecore(1,2))
c _____________________________________________________________________
c
         data etop/1.e-6/,thresh/1.e-7/,fndt/'foundit'/,tble/'trouble'/
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
           mesh=1001
           nsize=40
           nblock=mesh/nsize
c /////////////////////////////////////////////////////////////////////
c      icon  is zero  for starting from thomas-fermi atom
c      icon  is one   for continuing for potential stored in oldfl
c      icon  is two   for changing the electron configuration
c
c      zed is the atomic number
c
c      xmom is the moment;   up spins minus down spins
c
c      numsp is one   for a paramagnetic calculation
c
c      maxit is the number of iterations (negative for atom)
c
c      alfa is the mixing parameter ( i have been using .15 )
c /////////////////////////////////////////////////////////////////////
c
c _____________________________________________________________________
      call redin(icon,zed,xmom,numsp,maxit,alfa,iexc,iunin,iunot)
c _____________________________________________________________________
c
      write(6,1050)
      write(6,1503)
c
 1050  format(1x,73('=')/1x,'atomic code: local density (spherical)'
     >   ,3x,'written jan/feb 87 by f.j.pinski')
 1503 format(1x,73('.') )
c
       zup=.5*(zed+xmom)
c
c /////////////////////////////////////////////////////////////////////
c         xsetup contructs r grid
c         thofer sets up an intial approximation to the potential
c               and charge density using latter's fit of the
c               thomas-fermi atom
c =====================================================================
           call xsetup(x,r,mesh,nsize,zed,delta,sym,name)
c =====================================================================
c
c =====================================================================
           if(icon.eq.0) call thofer(r,rv,rhoold,zed,mesh,numsp)
c =====================================================================
c
c /////////////////////////////////////////////////////////////////////
              do 6 i=1,numsp
              zd=(i-1)*zed+(3-2*i)*zup
c =====================================================================
              call config(zd,numcor(i),ncore(1,i),lcore(1,i),occ(1,i)
     >                          ,ecore(1,i),icon,numsp,iunin)
c =====================================================================
    6         continue
c /////////////////////////////////////////////////////////////////////
c
           write(6,1005) zed,sym,zup,zed-zup
 1005  format(1x,73('=')/'  atomic number ',f6.2,8x,a2,2x,6('='),
     >      '====  atom:    up electrons=',f5.2/1x,39('='),
     >      '===========  down electrons=',f5.2/7x,a15,8x,a2,2x,
     >   /1x,73('=')  )
c
c /////////////////////////////////////////////////////////////////////
         if(icon.le.0) go to 17
c ----------------------------------- read old potential if  icon.gt.0
              read(iunin,3330) rv
              read(iunin,3330) (dummy,i=1,1001)
              read(iunin,3330) rhoold
              close (33)
 3330         format(4e20.13)
c ---------------------------------------------------------------------
              if(abs(rv(1,1)+2.*zed).gt.0.000001) then
                  write(6,*)' renormalizing potential:  to  ',zed
     >                                    ,'  from  ',-.5*rv(1,1)
                  if(icon.ne.2)
     >            write(6,*) ' ====== warning did not reoccupy levels'
                  renorm=-2.*zed/rv(1,1)
                  do 230 ns=1,numsp
                      do j=1,mesh
                      rv(j,ns)=renorm*rv(j,ns)
                      enddo
  230             continue
              else
c ====================================================================
                  if(alfa.lt.0) call potgen
     >                   (r,rv,rhoold,zed,mesh,nsize,numsp,iexc)
c ====================================================================
              endif
   17      continue
c ---------------------------------------  end of reading old potential
c /////////////////////////////////////////////////////////////////////
c =====================================================================
c --------------> main loop <--------------> main loop <---------------
c --------------> main loop <--------------> main loop <---------------
c --------------> main loop <--------------> main loop <---------------
c --------------> main loop <--------------> main loop <---------------
c --------------> main loop <--------------> main loop <---------------
c =====================================================================
c /////////////////////////////////////////////////////////////////////
        iter=0
   20   continue
        iter=iter+1
              if(icon.eq.0) then
                    itermx=max(1,iter-2)
              else
                    itermx=iter
              endif
c _____________________________________________________________________
        call savin(itermx,rv,numsp)
c =====================================================================
        write(6,1777) iter,abs(alfa)
 1777  format(1x,73('=')/ 1x,' iteration number:',i4,19x,
     >  'mixing parameter: ',f6.3 / 1x,73('=') )
c =====================================================================
        qsum=0.0
        do 35 ns=1,numsp
        elcnum(ns)=0.0
             do 21 k=1,1001
             rho(k,ns)=0.
   21        continue
c =====================================================================
             do 30 n=1,numcor(ns)
             en=ecore(n,ns)
             lambda=lcore(n,ns)
             nfl=ncore(n,ns)
             em=etop
c =====================================================================
c   bound -- solves the schrodinger equation: finds bound states
c =====================================================================
             call bound(abad,zed,r,rv(1,ns),mesh,en,thresh,
     >                        nfl,lambda,temp,kkk,em,q2,qi)
c =====================================================================
             ecore(n,ns)=en
c =====================================================================
             if(abad.eq.fndt) then
                  elcnum(ns) = elcnum(ns) + occ(n,ns)
                  qsum=qsum-(3-numsp)*occ(n,ns)*qi
                  do 22 k=1,kkk
                  snlo(k)=temp(k)
   22             continue
                  if(kkk.eq.mesh) go to 24
                  kk1=kkk+1
                  do 23 k=kk1,mesh
                  snlo(k)=0.0
   23             continue
   24             continue
                  xxx=occ(n,ns)
                  do 25 k=1,kkk
                  rho(k,ns)=rho(k,ns)+xxx*snlo(k)
   25             continue
             endif
c ---------------------------------------------------------------------
             tmo=(3-numsp)*occ(n,ns)
             write(6,1788) abad,ns,n,nfl,lambda,en,tmo
 1788        format(1x,a7,', spin=',i1,', state=',i2,', n=',i1,
     >         ', l=',i1,', energy=',f18.11,', oc=',f6.3)
c ---------------------------------------------------------------------
   30        continue
c ---------------------------------------------------------------------
   35   continue
c ---------------------------------------------------------------------
        elcore=elcnum(1)+elcnum(numsp)
        write(6,1200) name,zed,elcnum(1),elcnum(numsp)
 1200   format(1x,73('=')/1x,a15,'at.num=',f5.1,5x,f7.3,
     >  ' up   and',f7.3,'  down   electrons'/1x,73('='))
c =====================================================================
        call potgen(r,rvnew,rho,zed,mesh,nsize,numsp,iexc)
c =====================================================================
c
c                                       total energy section
c                                       total energy section
c                                       total energy section
c                                       total energy section
c                                       ====================
c
        enxc=0.0
        ecorr=0.0
        do 79 ns=1,numsp
c ---------------------------------------------------------------------
            do 32 i=1,mesh
            temp(i)=rho(i,ns)*(4.*excen(i)-3.*excpo(i,ns))
   32       continue
            call intqud(temp,tmp,r,nblock,nsize)
            enxc=enxc+(3-numsp)*tmp(mesh)
c
c ----------------------------- correction terms to make --------------
c ------------------------- the total energy variational --------------
c
            temp(1)=0.0
            do 56 i=2,mesh
            temp(i)=(rv(i,ns)-rvnew(i,ns))/r(i)
   56       continue
            call derv5(temp,tmp,r,mesh)
            do 57 i=2,mesh
            temp(i)=r(i)*rho(i,ns)*tmp(i)
   57       continue
            call intqud(temp,tmp,r,nblock,nsize)
            ecorr=ecorr+(3-numsp)*tmp(mesh)
   79   continue
c ---------------------------------------------------------------------
        call savout(itermx,rvnew,numsp)
c ---------------------------------------------------------------------
        tote=qsum+enxc
        write(6,1188) qsum,enxc,tote
 1188   format(1x,73('/')/1x,' qsum=',f12.6,3x,' excor=',f12.6,5x,
     >     '      energy=',f15.7)
c
        tote=tote+ecorr
        write(6,1189) ecorr,tote
 1189   format(1x,' corr=',f12.6,3x,'       ',12x,
     >      5x,'total energy=',f15.7/1x,73('/')/1x,73('_') )
c =====================================================================
c
c  ---------------      mix potentials     ------------------------
c  ---------------      mix potentials     ------------------------
c  ---------------      mix potentials     ------------------------
c  ---------------      mix potentials     ------------------------
c
            alfa=abs(alfa)
            rms=0.0
            vrms=0.0
            do 80 ns=1,numsp
            do k=1,mesh
               rms=rms+(rhoold(k,ns)-rho(k,ns))**2
               vrms=vrms+(rvnew(k,ns)-rv(k,ns))**2
               rhoold(k,ns)=rho(k,ns)
            enddo
   80       continue
            rms=sqrt(rms/(mesh*numsp))
            vrms=sqrt(vrms/(mesh*numsp))
c =====================================================================
            ndim=min(4,itermx)
            call dgand(rv,ndim,numsp,alfa)
c =====================================================================
            write(6,1787) iter,rms,iter,vrms
 1787       format(1x,9('='),' iter=',i3,1x,
     >      '==== rms difference in ch. density=',g11.4,1x,'======='
     >        /1x,9('='),' iter=',i3,1x,
     >      '==== rms difference in r*potential=',g11.4,1x,'=======' )
c =====================================================================
c =====   write info onto file   ======================================
c =====   write info onto file   ======================================
c =====   write info onto file   ======================================
c =====   write info onto file   ======================================
c =====================================================================
c
                write(iunot,*) numcup,occup
                write(iunot,*) ncoup,lcoup,ecoup
                if(numsp.eq.2)write(iunot,*) numcdn,occdn
                if(numsp.eq.2)write(iunot,*) ncodn,lcodn,ecodn
                write(iunot,3330) rv
                write(iunot,3330) r
                write(iunot,3330) rho
                rewind iunot
c =====================================================================
c
            if(iter.le.maxit) go to 20
c =====================================================================
c ============   end of main   ========================================
c ============   end of main   ========================================
c ============   end of main   ========================================
c ============   end of main   ========================================
c ============   end of main   ========================================
c =====================================================================
c
      stop 'Ok'
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine redin(icon,zed,xmom,nspins,nscfit,alpa,iexc,
     >             iunin,iunot)
c =====================================================================
        implicit real*8 (a-h,o-z)
        character*20 infile,otfile,blank,ptflin,ptflot
        character*78 title
        data blank/'                    '/
c
c       call dropfile(0)
c
        isave=6
c
c                     ntract=1   interactive mode for input
c                           =0   batch mode
c
c
        write(6,*)' atomic program, spherical atom, scf-total energy'
ccc     write(6,*)' do you wish to run this program in an interactive'
ccc  >       ,' mode? (1=yes,0=no)'
ccc     read(5,*) ntract
ccc     write(6,'('' ntract '',i5)') ntract
c
        read(5,'(a)')title
        ntract = 0
c
        if(ntract.eq.0) then
c
ccc         write(6,*)' input for this program will be read from a file'
ccc         write(6,*)' please type the name of this file. (format:a20)'
ccc         read(5,'(a20)') infile
ccc         write(6,'(a20)') infile
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
c
ccc      write(6,*) 'please type name of output file. (format:a20)'
ccc      write(6,*) 'if blank: the program will continue typing to the'
ccc  >             ,' screen'
         read(5,'(a20)') otfile
         write(6,'(2a)')'Output file name: ',otfile
c
        if(ntract.eq.0) then
ccc        open(unit=5,file=infile,form='formatted',status='old')
           if(otfile.ne.blank) then
              open(unit=6,file=otfile,form='formatted',status='unknown')
           endif
ccc        read(5,*) icon
           icon = 0
           read(5,*) zed
           read(5,*) xmom
           read(5,*) nspins
           if(nspins.eq.1) xmom=0
           read(5,*) iexc
           read(5,*) nscfit
           read(5,*) alpa
           read(5,'(a20)') ptflot
           if(icon.gt.0) read(5,'(a20)') ptflin
           close(5)
        else
            write(6,*) ' input: title of run.  (format:a78)'
            read(5,'(a78)') title
            write(6,*)' input: starting condition:  0 for fresh start'
            write(6,*)'                             1 for continuing '
            write(6,*)'                          or 2 for continuing '
     >          ,         'but changinging occupations'
            read(5,*) icon
            write(6,*)' input: atomic number'
            read(5,*) zed
            write(6,*)' input: moment'
            read(5,*) xmom
            write(6,*)' input: ex-corr. potential (1=vb-h,2=vosko)'
            read(5,*) iexc
            write(6,*)' input: number of spins (1 or 2)'
            read(5,*) nspins
            if(nspins.eq.1) xmom=0
            if(icon.gt.0) then
            write(6,*)' input: beginning potential'
            write(6,*)' please type the name of this file. (format:a20)'
            read(5,'(a20)') ptflin
            endif
            write(6,*)' input:     output potential'
            write(6,*)' please type the name of this file. (format:a20)'
            read(5,'(a20)') ptflot
            write(6,*) ' input: number of scf iterations'
            read(5,*) nscfit
            write(6,*) ' input: mixing parameter (typically .2)'
            read(5,*) alpa
        endif
c .....................................................................
        write(6,'(a78)') title
        write(6,*)' zed=',zed
        write(6,*)' icon=',icon
        write(6,*) ' nspins=',nspins
        if(nspins.eq.2) write(6,*)' xmom=',xmom
        write(6,*) 'iexc=',iexc
        write(6,*) ' number of scf iterations is ',nscfit,
     >            '     with mixing parameter alpha=',alpa
        if(icon.gt.0)
     >     write(6,'(1x,'' input potential file is '',a20)') ptflin
        write(6,'(1x,''output potential file is '',a20)') ptflot
c .....................................................................
c
        write(isave,'(a78)') title
        write(isave,
     &     '(1x,i4,t26,''continuing run (0=no,1=yes,2=yes+ch.occ)'')') 
     &     icon
        write(isave,'(1x,f9.5,t26,''zed: atomic number'')') zed
        write(isave,'(1x,f9.5,t26,''moment'')')xmom
        write(isave,'(i5,t26,''number of spins'')') nspins
        write(isave,'(i5,t26,''ex-cor type (1=vb-hedin,2=vosko)'')')
     >                             iexc
        write(isave,'(i5,t26,''number of iterations'')') nscfit
        write(isave,'(1x,f9.5,t26,''mixing parameter'')')alpa
        write(isave,'(a20,5x,''pot file (output)'')') ptflot
        if(icon.gt.0)
     >      write(isave,'(a20,5x,''pot file  (input)'')') ptflin
        if(isave.ne.6) close(isave)
c .....................................................................
        if(ntract.eq.1) then
              write(6,*)' do you wish to continue?  (1=yes,0=no)'
              read(5,*) icon1
              if(icon1.eq.0) stop
        endif
c .....................................................................
        if(icon.gt.0) then
           iunin=33
           open(unit=iunin,file=ptflin,form='formatted',status='old')
        endif
        iunot=34
        open(unit=iunot,file=ptflot,form='formatted',status='unknown')
c       write(6,*) iunin,iunot
c .....................................................................
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        subroutine xsetup(x,r,mesh,nsize,zed,delta,symb,name)
c =====================================================================
        implicit real*8 (a-h,o-z)
        character*2 symb,sym(103)
        character*15 name,nam(103)
        dimension x(1001),r(1001)
        parameter(three=3.,third=1./three)
c
c =====================================================================
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
c
c =====================================================================
c
        ized = zed+.1
        name=nam(ized)
        symb=sym(ized)
        write(6,1705) zed,symb,name
        write(6,1706)
 1705   format(1x,73('=')/5x,'atomic number ',f4.0,
     >     3x,'chemical symbol: ',a2,
     >     3x,'element: ',a15  )
 1706   format(1x,73('=')  )
c
c =====================================================================
 1002   format(11x,' mesh=',i5,8x,'nsize=',i4,8x,'nblock=',i4)
        nblock=mesh/nsize
        write(6,1002) mesh,nsize,nblock
        c=0.75
        dx=.00005
        i=1
        x(1)=0.
        r(1)=0.
        do 50 j=1,nblock
           dr=c*dx
           imk=(j-1)*nsize+1
           ximk=x(imk)
           do 40 k=1,nsize
                i=i+1
                x(i)=ximk+k*dx
                r(i)=c*x(i)
   40      continue
c          write(6,*) 'n,r(n),dr=',i,r(i),dr
           if(dr.lt.0.20)dx=dx+dx
   50   continue
c =====================================================================
        m3=mesh/150
        m3=150*m3
        write(6,1001) (k,r(k),k=1,m3,50)
 1001   format(3x,'r(',i4,')=',f13.8
     >        ,4x,'r(',i4,')=',f13.8
     >        ,4x,'r(',i4,')=',f13.8)
c =====================================================================
        if(i.le.mesh) return
        write(6,*)' trouble in setup:  mesh problem'
      stop
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine config(zed,nstate,nprin,l,fill,e,icon,numsp,iunit)
c =====================================================================
       implicit real*8 (a-h,o-z)
       character*1 nc(2,20)
       dimension nprnc(20),ll(20),nprin(20),l(20),fill(20),e(20)
       data nc/'1','s', '2','s', '2','p', '3','s',
     >         '3','p', '4','s', '3','d', '4','p',
     >         '5','s', '4','d', '5','p', '6','s',
     >         '4','f', '5','d', '6','p', '7','s',
     >         '5','f', '6','d', '7','p', '8','s'/
       data nprnc/1,2,2,3,3,4,3,4,5,4,5,6,4,5,6,7,5,6,7,8/
       data  ll  /0,0,1,0,1,0,2,1,0,2,1,0,3,2,1,0,3,2,1,0/
c
       nm=0
       nstate=0
       do 10 i=1,20
       l(i)=ll(i)
       nprin(i)=nprnc(i)
       e(i)=-3.*zed*zed/nprin(i)**6
       fill(i)=0.
   10  continue
   20  nstate=nstate+1
c =====================================================================
       if(nstate.gt.20) then
          write(6,*) ' stopping in subroutine config:  nstate problem'
          stop
       endif
c =====================================================================
       ndeg=2*ll(nstate)+1
       nm=nm+ndeg
       fill(nstate)=ndeg
       if(zed-nm.gt.0.0) go to 20
       dif=nm-zed
       fill(nstate)=ndeg-dif
c =====================================================================
           sum=0.0
           do 19 i=1,nstate
              sum=sum+fill(i)
   19      continue
c =====================================================================
           if(icon.ge.1) then
                read(iunit,*) nstate,fill
                read(iunit,*) nprin,l,e
c3330           format(4e20.13)
c =====================================================================
                sum=0.0
                jstate=nstate
                if(icon.gt.1)jstate=min(nstate+2,20)
                do 88 j=1,jstate
                if(icon.gt.1) then
                    write(6,*)' change occupancy?  state:',j,'   ',
     >              nc(1,j),nc(2,j),' n,l=(',nprnc(j),
     >              ',',ll(j),')  occ=',fill(j)
                    write(6,*) ' input occ:'
                    read(5,*) fill(j)
                    write(6,*) fill(j)
                    if(j.gt.nstate .and. fill(j).gt.0) nstate=j
                endif
c =====================================================================
                emmm=-.001
                e(j)=min(e(j),emmm)
                sum=sum+fill(j)
   88           continue
          endif
c =====================================================================
       zed=sum
       write(6,1000)   (3-numsp)*sum
       write(6,1001)(i,nc(1,i),nc(2,i),nprnc(i),ll(i),
     >    e(i),(3-numsp)*fill(i),i=1,nstate)
       write(6,1002)
c =====================================================================
c
 1000 format(/1x,23('+'),' number of electrons:',f7.3,1x,21('+')
     >   /1x,73('=')/29x,'state configuration'/1x,73('-') )
 1001  format(1x,i2,'.',2a1,':n=',i1,',l=',i1,',e=',f9.3,',occ=',f5.2
     > ,      2x,i2,'.',2a1,':n=',i1,',l=',i1,',e=',f9.3,',occ=',f5.2)
 1002 format(1x,73('=') )
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine thofer(r,v,rho,zed,mesh,numsp)
c =====================================================================
       implicit real*8 (a-h,o-z)
       parameter (pi=3.14159265358979,rf=4./(3.*pi) )
       parameter(three=3.,oneth=1./three)
       dimension r(1001)
       dimension v(1001,2),rho(1001,2)
c
c     v   is r*potential
c    rho  is 4pi r**2 * density
c =====================================================================
c
       scale=zed**oneth/.88534138
c
       v(1,1)=-2.*zed
       rho(1,1)=0.
       v(1,numsp)=-2.*zed
       rho(1,numsp)=0.
c =====================================================================
       do 10 i=2,mesh
       rr=r(i)*scale
       root= sqrt(rr)
       rvozin=1.            + 0.02747 *root
     >      + 1.243  *rr    - 0.14860 *rr*root
     >      + .2302  *rr*rr + 0.007298*rr*rr*root
     >      + .006944*rr*rr*rr
       v(i,1)=2.*zed/rvozin
       rho(i,1)=rf*v(i,1)* sqrt(v(i,1)*r(i))
       v(i,1)=-v(i,1)
       if(numsp.le.1) go to 10
          v(i,2)=v(i,1)
          rho(i,2)=rho(i,1)
   10  continue
c =====================================================================
c
c       write(6,1001) (r(i),i,v(i,1),i=1,mesh,100)
c1001   format(1x,'r=',e11.5,1x,'v(',i4,')=',e11.5
c    > ,       3x,'r=',e11.5,1x,'v(',i4,')=',e11.5)
c       write(6,1002) (r(i),i,rho(i,1),i=1,mesh,100)
c1002   format(1x,'r=',e11.5,1x,'n(',i4,')=',e11.5
c    > ,       3x,'r=',e11.5,1x,'n(',i4,')=',e11.5)
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      subroutine intqud(y,yi,x,nb,nin)
c =====================================================================
      implicit real*8 (a-h,o-z)
      dimension x(1001),y(1001),yi(1001)
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
      implicit real*8 (a-h,o-z)
      dimension r(1001),y(1001),yi(1001)
      bsum=0.
      asum=0.
      yi(1)=0.
      i=1
      kmax=nin/2
c                             nin must be an even integer
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
       subroutine potgen(r,rv,rhos,zed,mesh,nsize,numsp,iexcor)
c =====================================================================
       implicit real*8 (a-h,o-z)
       parameter (onemd=.99999999999,epsm=1.e-19)
       parameter (three=3.,third=1./three)
c
       dimension r(1001),rhos(1001,2),rv(1001,2)
       dimension dcoul(1001),excpo(1001,2),dexcpo(1001,2)
       dimension q(1001),q1(1001),coul(1001),drv(1001,2)
       dimension rho(1001),xm(1001),temp(1001)
       common/drvdr/drv
       common/excorr/excpo,excen(1001)
c =====================================================================
c
       nblock=mesh/nsize
       rho(1)=0.
       temp(1)=0.
       xm(1)=0.
       do 30 i=2,mesh
         rho(i)=rhos(i,1)+rhos(i,numsp)
         temp(i)=rho(i)/r(i)
         xm(i)=rhos(i,1)-rhos(i,numsp)
   30  continue
c =====================================================================
c
       q(1)=0.
       call intqud(rho,q,r,nblock,nsize)
       write(6,*)' total number of electrons: ',q(mesh)
c
       q1(1)=0.
       call intqud(temp,q1,r,nblock,nsize)
c =====================================================================
       coul(1)=-2*zed
       dcoul(1)=0.0
       temp(1)=0.0
       do 62 ns=1,numsp
         sp=3-2*ns
         rv(1,ns)=coul(1)
         do 60 i=2,mesh
            rs=((3.*r(i)**2)/max(rho(i),epsm))**third
            dz=(numsp-1)*xm(i)/max(rho(i),epsm)
            dz=min(onemd,max(-onemd,dz))
            if(ns.eq.1) then
                        coul(i)=-2*(zed-q(i)+r(i)*(q1(i)-q1(mesh)))
                        dcoul(i)=-2*(q1(i)-q1(mesh))
            endif
            excpo(i,ns)=alpha2(rs,dz,sp,iexcor,edum)
            temp(i)=r(i)*excpo(i,ns)
            excen(i)=edum
            rv(i,ns)=coul(i)+temp(i)
   60    continue
         call derv5(temp,dexcpo(1,ns),r,mesh)
         do 61 i=1,mesh
            drv(i,ns)=dcoul(i)+dexcpo(i,ns)
   61    continue
   62  continue
c =====================================================================
c
        m2=mesh/200
        m2=m2*200
        do 77 ns=1,numsp
        write(6,1001) (r(i),i,rv(i,ns),i=1,m2,100)
 1001   format(1x,'r=',e14.6,1x,'v(',i4,')=',e14.6
     > ,       3x,'r=',e14.6,1x,'v(',i4,')=',e14.6)
c       write(6,1002) (r(i),i,rhos(i,ns),i=1,m2,100)
c1002   format(1x,'r=',e14.6,1x,'n(',i4,')=',e14.6
c    > ,       3x,'r=',e14.6,1x,'n(',i4,')=',e14.6)
   77    continue
c =====================================================================
c             if(numsp.eq.1) return
c =====================================================================
c     write(6,1008)(r(i),i,xm(i)/max(rho(i),epsm)
c    >                 ,i=1,m2,100)
c1008   format(1x,'r=',e12.6,1x,'p(',i4,')=',e12.6
c    > ,       3x,'r=',e12.6,1x,'p(',i4,')=',e12.6)
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
      function alpha2 (rs,dz,sp,iexcor,exchg)
c =====================================================================
c
      implicit real*8 (a-h,o-z)
      parameter(one=1.,two=2.,three=3.,thrd=1./three,for3=4./three)
      parameter (rsmax=250.)
c
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
      data x0/-0.0047584,-0.10498,-0.3250/
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
      alpha2=0
      exchg=0
      if(rs.gt.rsmax) return
c
      go to (10,20) iexcor
c
c     von barth-hedin  exch-corr potential
c     j. phys. c5,1629(1972)
c
  10  continue
      fm=two**for3-two
      fdz=( (one+dz)**for3 + (one-dz)**for3 - two )/fm
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
   1  continue
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
         subroutine derv5(y,dy,x,nn)
c =====================================================================
         implicit real*8 (a-h,o-z)
c        parameter (ipmesh=641)
         dimension x(nn),y(nn),dy(nn)
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
        subroutine bound(abad,zed,r,v,mesh1,en,thresh,
     >                      nq,lq,chden,kkk,emn,qi2,qii)
c =====================================================================
        implicit real*8 (a-h,o-z)
        character*7 abad,fndt,tble
        parameter (onepd=1.00001,eminup=-.0005)
        common/e1/fndt,tble
        dimension r(1001),v(1001),chden(1001),temp(1001),qi2(1001)
        dimension vme(1001),y(1001),dy(1001),yin(1001),dyin(1001)
c =====================================================================
       abad=tble
       l=lq
       n=nq
       mesh=mesh1
       if(en.ge.0.0) en=-zed**2/n**4
       emin=min(emn,eminup)
       energy=en
       en=emin
c      write(6,*) ' in bound: emin=',emin
       eold=99.99
c =====================================================================
c zero out chden
        do 3 i=1,mesh1
        qi2(i)=0.0
        chden(i)=0.0
    3   continue
c =====================================================================
        vme(1)=0.0
        do 4 j=2,mesh
        vme(j)=(v(j)+(l*l+l)/r(j))/r(j)   - emin
    4   continue
        do 6 i=2,mesh
        j=i
        if(r(i)*sqrt(-emin) .gt. 150.)go to 7
   6    continue
   7    lll=j
c =====================================================================
c       write(6,*)' in bound: limit=',lll
        call outsch(zed,emin,l,lll,vme,v,r,y,dy)
c =====================================================================
        num=0
        do 8 i=3,lll
           if(y(i).eq.00.0) then
                num=num+1
                go to 8
           else
                if(y(i-1).eq.0) go to 8
                is=sign(onepd, y(i) )
                js=sign(onepd,y(i-1))
                if(is.ne.js) num=num+1
           endif
   8    continue
c =====================================================================
c number of nodes  n-l-1
c
c      write(6,*)' found ',num,' nodes'
       nnode=n-l-1
       if(nnode.ge.num) then
               write(6,1740) emin,nq,lq
 1740          format('  no bound state below energy=',f8.4,
     >         ', with quantum numbers: n=',i2,', l=',i2 )
               return
       endif
c =====================================================================
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
c =====================================================================
c ....................  main loop .....................................
c =====================================================================
  10   continue
       iter=iter+1
       if(iter.gt.maxit) go to 800
       if(energy.gt.0) energy=eold/2.
       bigkr=150.
       if(abs(energy).gt.400.) bigkr=125.
       bigkr=.5*bigkr
c
c find classical turning point ........................................
c
        do 40 j=2,mesh
        vme(j)=(v(j)+(l*l+l)/r(j))/r(j)   - energy
  40    continue
c =====================================================================
        do 42 k=3,mesh
        j=mesh+1-k
        if(vme(j).le.0.00) go to 45
   42   continue
c =====================================================================
        write(6,*) ' cannot find classical turning point '
        if(iter.eq.1) then
               iter=0
               jter=jter+1
               if(jter.gt.maxjt) then
                    write(6,*) ' turning point trouble'
                    return
               endif
               energy=.75*energy
               go to 10
         else
               iter=iter-1
               jter=jter+1
               if(jter.gt.maxjt) then
                    write(6,*) ' turning point trouble'
                    return
               endif
               energy=.5*(energy+eold)
          endif
c =====================================================================
  45    match=j+1
c =====================================================================
        match=min(match,mesh-40)
        do 60 i=match,mesh
        j=i
        if(r(i)*sqrt(-energy) .gt. bigkr ) go to 61
  60    continue
  61    last=j
        if(iext.eq.0) lll=last
        if(iext.ne.0) lll=match+1
        mmmm=match-1
c =====================================================================
c    outward integration of the schrodinger equation
c =====================================================================
        call outsch(zed,energy,l,lll,vme,v,r,y,dy)
c =====================================================================
        gamout=dy(match)/(sqrt(-energy)*y(match))
c                                                 count number of nodes
        num=0
        do 80 i=3,lll
           if(y(i).eq.0.00) then
                num=num+1
                go to 80
           else
                if(y(i-1).eq.0) go to 80
                is=sign(onepd, y(i) )
                js=sign(onepd,y(i-1))
                if(is.ne.js) num=num+1
           endif
  80    continue
c =====================================================================
        if(iter.eq.1) then
             iext=0
             if(num.eq.nnode+1) then
                  enext=1.25*energy
                  jter=0
                  go to 87
             endif
             if(num.eq.nnode) then
                  enext=.75*energy
                  jter=0
                  go to 87
             endif
             jter=jter+1
             if(jter.gt.maxjt)  then
                  write(6,*) ' cannot find right number of nodes '
                  write(6,*) ' looking for ',nnode,' nodes'
                  write(6,*) '       found ',num,' nodes'
                  return
             endif
             iter=0
             if(num.lt.nnode) then
                  energy=.5*energy
                  go to 10
             else
                  energy=1.5*energy
                  go to 10
             endif
         endif
c =====================================================================
         if(iter.eq.2)  then
             if(num.eq.nold) then
                  iter=1
                  jter=jter+1
                  if(jter.gt.maxjt) then
                       write(6,*) ' trouble in second iteration'
                       return
                  endif
                  dir=+1.
                  if(num.eq.nnode) dir=-1.
                  enext=energy-.8*dir*abs(energy-eold)
                  if(enext.ge.0.00) enext=.5*energy
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
c =====================================================================
         if(abs(eold-energy).lt.thresh) go to 200
c =====================================================================
         if(num.ne.nnode) go to 83
   81    dir=+1
         go to 85
   83    dir=-1.
   85    enext=energy+.5*dir*abs(energy-eold)
   87    jter=0
   88    nold=num
         if( abs(1.-eold/energy) .lt. 0.1  .or. iext.eq.1) then
c =====================================================================
c   inward integration of schrodinger equation
c =====================================================================
             call inwsch(zed,energy,mmmm,last,l,vme,r,v,yin,dyin)
c =====================================================================
             gamin=dyin(match)/(sqrt(-energy)*yin(match))
             rat=(1.+gamin)/(1.-gamin)
             rat=rat-((1.+gamout)/(1.-gamout))
             if(rat.eq.0) go to 200
             if(iext.eq.1) then
                 enext=energy+(eold-energy)/( 1. - ratold/rat)
             endif
             iext=1
             ratold=rat
         endif
c =====================================================================
         eold=energy
         energy=enext
         go to 10
c
  800    continue
c =====================================================================
c ..........  can not find state  .....................................
c =====================================================================
         write(6,*) ' too many tries for n=',nq,' l=',l
         return
c
  200    continue
c =====================================================================
c .................... found state ............................
c =====================================================================
         call inwsch(zed,energy,mmmm,last,l,vme,r,v,yin,dyin)
c =====================================================================
         gamin=dyin(match)/(sqrt(-energy)*yin(match))
         abad=fndt
         kkk=last
         do 225 i=1,match
         chden(i)=y(i)**2
  225    continue
         anorm=y(match)/yin(match)
         do 226 i=match,last
         y(i)=anorm*yin(i)
         dy(i)=anorm*dyin(i)
         chden(i)=(anorm*yin(i))**2
  226    continue
         en=energy
c =====================================================================
         nsize=40
         nblock=mesh1/nsize
         call intqud(chden,temp,r,nblock,nsize)
         qi2(1)=0.0
         do 227 i=2,mesh1
         if(i.lt.last) then
                chden(i)=chden(i)/temp(mesh1)
                qi2(i)=(dy(i)-(l+1)*y(i)/r(i))**2/temp(mesh1)
         else
                chden(i)=0.0
                qi2(i)=0.0
         endif
  227    continue
         call intqud(qi2,temp,r,nblock,nsize)
         qii=temp(mesh1)
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
        subroutine outsch(zed,energy,l,match,vme,v,r,y,dy)
c =====================================================================
        implicit real*8 (a-h,o-z)
        dimension vme(1001),v(1001),r(1001)
        dimension y(1001),dy(1001)
        dimension fac(5)
        data fac/1.,3.,15.,105.,945./
c =====================================================================
        if(l.lt.0 .or. l.gt.4) then
             write(6,*) ' angular momentum quatum number is ',l
             write(6,*) ' stopping in outsch: l out of bounds'
             stop
        endif
c =====================================================================
        last=match
        e=energy
c
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
c =====================================================================
      delold=r(6)-r(5)
      ityp=0
      do 40 n=6,last
         del=r(n)-r(n-1)
c =====================================================================
         if(abs(del-delold).gt.1.e-7)  then
                 ityp=1
                 delold=del
                 m4=n-7
                 m3=n-5
                 m2=n-3
                 m1=n-1
         else
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
c =====================================================================
        h24=del/24.
c
        dy(n)=dy(m1)+h24*(55.*vme(m1)*y(m1) -59.*vme(m2)*y(m2)
     >                  + 37.*vme(m3)*y(m3) - 9.*vme(m4)*y(m4))
        do 37 ijk=1,5
        y(n)=y(m1)+h24*(   9.*dy(n)  +19.*dy(m1) -5.*dy(m2) + dy(m3))
        dy(n)=dy(m1)+h24*( 9.*vme(n)*y(n) +19.*vme(m1)*y(m1)
     >                   - 5.*vme(m2)*y(m2) +  vme(m3)*y(m3))
   37  continue
        y(n)=y(m1)+h24*(  9.*dy(n) +19.*dy(m1) -5.*dy(m2) + dy(m3))
c =====================================================================
   40  continue
c
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine inwsch(zed,energy,match1,last1,l,vme,r,v,y,dy)
c =====================================================================
        implicit real*8 (a-h,o-z)
        dimension r(1001),v(1001),vme(1001),y(1001),dy(1001)
c =====================================================================
        match=match1
        last=last1
        if(mod(last-2,40).lt.4) last=min(last+6,1001)
c
        n=last
        p=sqrt(-energy)
c =====================================================================
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
c =====================================================================
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
c =====================================================================
         if(abs(del-delold).gt.1.e-7) then
             a1=dy(m1)
             a2=dy(m2)
             a3=dy(m3)
             call qrfit(delold,a1,a2,a3,f1,f2)
             dym1=a1
             dym2=f1
             dym3=a2
             dym4=f2
c
             a1=y(m1)
             a2=y(m2)
             a3=y(m3)
             call qrfit(delold,a1,a2,a3,f1,f2)
             ym1=a1
             ym2=f1
             ym3=a2
             ym4=f2
c
             a1=vme(m1)
             a2=vme(m2)
             a3=vme(m3)
             call qrfit(delold,a1,a2,a3,f1,f2)
             vmem1=a1
             vmem2=f1
             vmem3=a2
             vmem4=f2
c
             delold=del
c =====================================================================
        else
c =====================================================================
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
c =====================================================================
c
         h24=-del/24.
         dy(n)=dym1 + h24 * (55.*vmem1*ym1-59.*vmem2*ym2
     >                      +37.*vmem3*ym3- 9.*vmem4*ym4)
         do 70 ijk=1,5
         y(n)=ym1 + h24*(9.*dy(n)+19.*dym1-5.*dym2+dym3)
         dy(n)=dym1 + h24 * ( 9.*vme(n)*y(n)+19.*vmem1*ym1
     >                       -5.*vmem2*ym2+vmem3*ym3)
   70    continue
         y(n)=ym1 + h24*(9.*dy(n)+19.*dym1-5.*dym2+dym3)
c =====================================================================
   80 continue
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine qrfit(del,a1,a2,a3,f1,f2)
c ======================================================================
       implicit real*8 (a-h,o-z)
       real*8 del,del2,delh
c=======================================================================
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
c =====================================================================
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
       subroutine savout(iter,vr,numsp)
c =====================================================================
       implicit real*8 (a-h,o-z)
       dimension vr(1001,2)
       common/vout/vrout(1001,2,4)
c
       n=mod(iter-1,4)+1
       do 10 ns=1,numsp
       do i=1,1001
       vrout(i,ns,n)=vr(i,ns)
       enddo
   10  continue
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine savin(iter,vr,numsp)
c =====================================================================
       implicit real*8 (a-h,o-z)
       dimension vr(1001,2)
       common/vin/vrin(1001,2,4)
c
       n=mod(iter-1,4)+1
       do 10 ns=1,numsp
       do i=1,1001
       vrin(i,ns,n)=vr(i,ns)
       enddo
   10  continue
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       function dpdif(i,j,numsp)
c =====================================================================
       implicit real*8 (a-h,o-z)
       common/vin/vrin(1001,2,4)
       common/vout/vrout(1001,2,4)
c
       sum=0.0
       do 10 ns=1,numsp
       do k=1,1001
           sum=sum+(vrin(k,ns,i)-vrout(k,ns,i))
     >      *k*(vrin(k,ns,j)-vrout(k,ns,j))
       enddo
   10     continue
       dpdif=sum
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine dgand(vrnew,ndim,numsp,alpa)
c =====================================================================
       implicit real*8 (a-h,o-z)
       common/vin/vrin(1001,2,4)
       common/vout/vrout(1001,2,4)
       dimension a(5,5),b(5),vrnew(1001,2)
c
       jdim=ndim+1
       do 2 i=1,ndim
         a(i,jdim)=1.
         a(jdim,i)=1.
         b(i)=0.
    2  continue
c
       a(jdim,jdim)=0.
       b(jdim)=1.
       do 5 j=1,ndim
         jj=j
         do i=1,j
           ii=i
           a(i,j)=dpdif(ii,jj,numsp)
           if(i.ne.j) a(j,i)=a(i,j)
         enddo
    5  continue
c
       call mxlneq(a,b,jdim)
c
       do 20 ns=1,numsp
       do k=1,1001
          sum=0.0
          do 10 i=1,ndim
          sum=sum+b(i)*(vrin(k,ns,i)+alpa*(vrout(k,ns,i)-vrin(k,ns,i)))
   10     continue
          vrnew(k,ns)=sum
       enddo
   20  continue
      return
      end
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
       subroutine mxlneq(a,y,nn)
c =====================================================================
       implicit real*8 (a-h,o-z)
       dimension a(5,5),y(5)
c
       n=nn
       do 10 i=2,n
            f1=-1./a(i-1,i-1)
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
          det=1.
          do 19 i=1,n
          det=det*a(i,i)
  19      continue
c
c         write(6,1235) n,det
c1235     format(' in mxlneq: order, det',2x,i4,e17.9)
c
      y(n)=y(n)/a(n,n)
      do 20 ij=1,n-1
          i=n-ij
          do 15 j=1,ij
          y(i)=y(i)-y(i+j)*a(i,i+j)
   15     continue
      y(i)=y(i)/a(i,i)
   20  continue
c
c
      return
      end
c _____________________________________________________________________
c ////////////////////////////////////////end of program///////////////
c _____________________________________________________________________
c /////////////////////////////////////////////////////////////////////
