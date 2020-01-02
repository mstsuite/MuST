      program murn
c
c     modified              : bk  3.1.93     
c     modified from murn2.f : r.s 19.04.91
c     modified from murn1.f : gps 6.12.90
c     least-squares fit of e vs v by murnaghan equation of state
c
      implicit double precision (a-h,o-z)
      double precision kb
      logical tzero
      dimension edata(50),vdata(50),x(40),aux(40)
      external fun

c conversion factor from rydberg to eV, Boltzmann constant in Ha:
      data convyy, kb /27.211616839, 3.16666225714e-6/
      common /c/ edata,vdata, volmin,volmax,b0min,b0max,b0pmin,b0pmax,
     $     ula,nnn

      in=5
      iout=6
      ula=0.52917715
      write(iout,'("# unit of length =",f15.6,"  angstroem")') ula

      read(in,*)iunit
      if(iunit .eq. 1) then 
        conv_inp = 27.2116168391d0/2.d0
        write(iout,1070)
      else if(iunit .eq. 2) then 
        conv_inp = 1.d0
        write(iout,1080)
      else if(iunit .eq. 3) then
        conv_inp = 27.2116168391d0
        write(iout,1081)
      else
        write(iout,1090)
        stop
      end if

 1070 format('# energies input in rydbergs')
 1080 format('# energies input in electronvolts')
 1081 format('# energies input in hartrees')
 1090 format('# only iunit = 1, 2 and 3 are recognized; stop.')

c     conversion factor from a**3 to the volume of unit cell
c     (e.g. 0.25, if the structure is such that v = a**3/4 ):
      read(in,*) convav
      write(iout,1020) convav
 1020 format('# unit cell volume = a**3 *',f15.6)
c
c     read in alat boundaries and number of points 
      read(in,*) alat_min, alat_max, nr_alat
      write(iout,
     &  '("# a(min)=",f20.7/,"# a(max)=",f20.7/,"# points:",i4)') 
     &  alat_min, alat_max, nr_alat
      if (nr_alat .lt. 2 .or. alat_max .le. alat_min) stop '# something w
     &  rong with this range'

c     Zero point lattice energy correction by Moruzzi-Janak-Williams
c     Some values:
c             Gamma        T(Debye)      Volume/atom for T(Debye)
c     Al:    2.19  428  109.9
c     Fe:    1.66  467  78.95
c     Cu:    2.00  343  78.92

c     read in the values for zero point vibration energy correction - apsi
      read(in,*) tzero, gamma, thetad, atvol

      if ( tzero ) then
        write(iout,'(a)')
     $       "#  Zero point lattice energy by Moruzzi-Janak-Williams:"
        write(iout,'(a,f7.3,/,a,f6.2,/,a,f8.2)')
     $       "#  Grueneisen constant = ", gamma,
     $       "#  T(Debye) = ", thetad,
     $       "#  Volume/at for T(Debye) = ", atvol
        write(iout,*)
        write(iout,'(a)')
     $       "#  Tested only for fcc lattice (Al)"
        write(iout,*)
      end if

c number of data points a, e(a):
      read(in,*) nnn
      write(iout,1030) nnn
 1030 format('# number of data points: n =',i5)
      if (nnn .gt. 50) then
        write(iout,*)' n too large, nmax=50'
        stop
      endif

      write(iout,1120)
      write(iout,1040)
 1040 format('#    i',5x,'alatt',9x,'volume',12x,'E (ry)')
      do 10 i=1,nnn
        read(in,*)alatt,eninp
        alatt=alatt
        vdata(i)=alatt**3*convav

c     Zero point vibrational correction - see above
        if ( tzero ) then
          eninporig = eninp
          eninp = eninp + 9.0/8.0 * kb*thetad * (atvol/vdata(i))**gamma
        end if

        edata(i)=eninp*conv_inp
        if ( tzero ) then
          write(iout,1115)i,alatt,vdata(i),eninp,eninporig
 1115     format('#',1x,i5,2f15.6,1f19.10,2x,1f19.10)
        else
          write(iout,1110)i,alatt,vdata(i),eninp
 1110     format('#',1x,i5,2f15.6,1f19.10)
        end if
   10 continue
      write(iout,1120)
 1120 format('#',10(1h-))
c
c     starting values for iteration:
c
      write(iout,1130)
 1130 format('#',5x,'alatt0',3x,'vol0 (ula**3)',2x,'b0 (mbar)',
     &  5x,'b0prime',6x,'e0 (ry)')

c mimimum in edata(i):
      e0min=1.d6
      volmin=1.d6
      do 100 i=1,nnn
        if(edata(i) .ge. e0min) go to 100
        e0min=edata(i)
        volmin=vdata(i)
  100 continue

      vol0=volmin
      e0ev=e0min
      e0ry=e0*convyy
      alatt0=(vol0/convav)**(1.d0/3.d0)
c     
c for other variables, we choose:
c     
      b0=1.d0
      b0prim=4.d0

      write(iout,1140)alatt0,vol0,b0,b0prim,e0ry
 1140 format('# ',f9.4,3f13.5,1f13.5)

      almin=0.5d0*alatt0
      almax=1.5d0*alatt0
      b0min=0.01d0
      b0max=10.d0
      b0pmin=0.01d0
      b0pmax=10.d0

      volmin=almin**3.d0*convav
      volmax=almax**3.d0*convav
      write(iout,1210)almin,volmin,b0min,b0pmin,almax,volmax,b0max,b0pm
     &  ax
 1210 format('# min:',f9.4,3f13.5/,'# max:',f9.4,3f13.5)
c     
c     a uniform shift of all the energies
c     (will not influence the b0, b0prime, vol0, only shifts e0):
c     
      shift=-e0ev-1.d0
      do 20 i=1,nnn
   20 edata(i)=edata(i)+shift

c murnaghan least-squares fit:
      write(iout,1120)
      write(iout,1280)
 1280 format('# fit by murnaghan equation equation of state')
      write(iout,1120)
      write(iout,1190)
 1190 format('# iter, alatt0, vol0 (ula**3),' 
     1  ,' b0 (mbar), b0prime, e0 (ev), sq sum, ierr')

      x(1)=e0ev+shift
      x(2)=b0
      x(3)=b0prim
      x(4)=vol0
      nvar=4
      lim=20

      do 70 i = 1, 75
        call dfmnd(fun,x,fff,nvar,lim,aux,ierr)
        write(iout,'("# ",i3,f9.4,5f12.5,i2)')
     &       20*i,ula*(x(4)/convav)**(1.d0/3.d0),x(4),x(2),x(3),
     &       (x(1)-shift)/convyy,fff,ierr
        if (ierr .eq. 0) go to 80
   70 continue
   80 continue

        write(iout,'("bb ",i3,f9.4,5(",",f12.5),i2)')
     &       20*i,ula*(x(4)/convav)**(1.d0/3.d0),x(4),x(2),x(3),
     &       (x(1)-shift)/convyy,fff,ierr

      write(iout,1120)
c     write(iout,1130)
      vol0=x(4)
      alatt0=(vol0/convav)**(1.d0/3.d0)
      b0=x(2)
      b0prim=x(3)
      e0ev=x(1)-shift
      e0ry=e0ev/convyy
c     write(iout,'("# alat=",f10.5\,"# b0=", f10.4, " b0p=", f10.3," e0=",
c     &  f11.6)') 
c     &  alatt0, b0, b0prim, e0ev / 27.212
c     write(iout,1140)alatt0,vol0,b0,b0prim,e0ev,e0ry
c     write(iout,1120)
c     
c     print out fitted energy values for nr_alat lattice
c     
      open (3,FILE='murn.gnu',STATUS='UNKNOWN',FORM='FORMATTED')
      do i=1, nnn
        write (3,'(f10.5,f18.6,f16.8)')
     &       (vdata(i)/convav)**(1.d0/3.d0),
     &       (edata(i)-shift)/conv_inp,
     &       (edata(i)-shift)/conv_inp-e0ev/conv_inp
      end do
      write(iout,500)
  500 format('# a (a.u.), E (ryd)')
      write(iout,1120)
      do i=1,nr_alat
c	vol=0.6d0*vol0 + (i-1)*0.05*vol0
c	alatt=(vol/convav)**(1./3.)
c	alatt=alatt*ula
	alatt=alat_min + (i-1)*(alat_max-alat_min)/(nr_alat-1)
        vol = alatt**3 * convav
	alatt=alatt*ula
	call murng1(ula,vol,vol0,b0,b0prim,e0ev,etot)
	etotry=etot/convyy
	write(iout,'(f10.5,f18.6,f10.6)') 
     &    alatt/ula,etot/conv_inp,(etot-e0ev)/conv_inp
c	write(iout,502) alatt, etotry, alatt/ula, etotry/2
c 502 format(2(f10.5,f20.7))
      enddo
c
      stop
      end
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      double precision function fun(x)
c     function to be minimized in the least-squares fit (by subroutine dfmn
c     d)
c     function
c     calculates the sum of the  a b s o l u t e  deviations
c            (e(theor)-e(exp))**2
c     divided by the number of exp. points,
c     assuming for equation of state the murnaghan expression.
c     
c     meaning of variables:
c     x(1) .... e0
c     x(2) .... b0
c     x(3) .... b0prim
c     x(4) .... vol0
c     
c
      implicit double precision (a-h,o-z)
      dimension edata(50),vdata(50),x(40)
      common /c/ edata,vdata, volmin,volmax,b0min,b0max,b0pmin,b0pmax,
     $     ula,nnn
c
c     *         *         *         *         *        *
      e0=x(1)
      b0=x(2)
      b0prim=x(3)
      vol0=x(4)
c
      sum=0.d0
c     the sum of squares:
      do 10 i=1,nnn
        volact=vdata(i)
        call murng1(ula,volact,vol0,b0,b0prim,e0,etot)
        sum=sum+(etot-edata(i))**2
   10 continue
      fun=sum/dfloat(nnn)
      return
      end
c     -----------------------------------------------------------------
      subroutine murng1(ula,vol,vol0,b0,b0prim,e0,etot)
c     evaluation of the murnaghan expression for energy as a function
c     of volume.
c     
c     input data:
c     
c     ula ..... unit of length, in angstroems, used here.
c     vol ..... volume, in the above units of length cubed.
c     vol0 .... volume at the energy minimum, in the same units.
c     b0 ...... bulk modulus, in units megabar.
c     b0prim .. pressure derivative of the bulk modulus.
c     should be close neither to 0 nor to 1.
c     e0 ...... an additive constant (in electronvolts), added
c     to the energy expression.
c     (see,  pr b28, p 5484: fu and ho)
c     
c     output data:
c     
c     etot .... the energy, including the e0 constant, in electronvolts
c     
c     if b0 differs from 0. or from 1. by less than 1.d-6, then
c     etot is set at +111111.111111 electronvolts.
c
c     *      *      *      *      *      *      *      *
c     
      implicit double precision (a-h,o-z)
c     
c     conversion factor from erg to electronvolt:
      parameter( convxx = 1.60209d0 )
c     1 electronvolt = 1.60209 d-12 erg
c     
c     
      if(dabs(b0prim)-1.d0 .lt. 1.d-6   .or.
     1  dabs(b0prim) .lt. 1.d-6) then
        etot=+111111.111111d0
        return
      end if
c     
      if(vol .lt. 0.d0   .or.   vol0 .lt. 0.d0) then
        etot=+111111.111111d0
        return
      end if
c     
      etot = e0 - b0*vol0/b0prim *
     1  (((vol/vol0)**(1.d0-b0prim)-1.d0)/(1.d0-b0prim)-vol/vol0+1.d0)
     2  *ula**3/convxx
c     
      return
      end
c     --------
      subroutine dfmnd(f,x,y,n,lim,aux,ier)
c     
c     ******************************************************************
c     *   minimization of a function of several variables              *
c     *   using powell's algorithm without derivatives                 *
c     ******************************************************************
c     
      double precision f,x,y,aux,eps,eta,tol,
     1  da,db,dc,dm,dq,fa,fb,fc,fl,fm,fs,hd,hq,hx
      dimension x(1),aux(1)
c     
c     subroutines required: the external function f.
c     
c     input data:
c     f .... the function of n variables x(1)...x(n) to be minimized
c     x .... x(1) ... x(n) = starting guess for the variables;
c     beware: x has to be dimensioned in the main program
c     to considerably more than n.
c     n .... number of variables; the dimension of x and aux in the
c     calling program must be, however, much higher:
c     - perhaps 10 times?
c     lim .. maximum number of iterations
c     aux .. auxiliary array of the same dimension as x.
c     output data:
c     x .... x(1) ... x(n) the resulting minimum
c     y .... value of the function at the minimum
c     ier .. some error indication 
c     ierr=0 means 'convergence achieved'.
c     
c     *      *      *      *      *      *      *      *
c     
      isw  =ier
      ier  =0
      if (n) 1,1,2
    1 ier  =1000
      goto 109
    2 if (lim) 3,3,4
    3 ier  =2000
      goto 109
c     
c     set initial values and save initial argument
c
    4 n1   =n+1
      n2   =n+n
      n3   =n+n2
      n4   =n+n3
      n5   =n*n+n4
      eps  =1.d-15
      eta  =n*eps
      do 5 k=1,n
        aux(k)=x(k)
        j    =n3+k
    5 aux(j)=1.d0
      fs   =f(aux)
      fb   =fs
      i    =1
      ic   =1
      it   =0
      m    =n4
      mf   =0
      is   =1
c     
c     start iteration cycle
c     
    6 it   =it+1
      fl   =fs
      do 7 k=n1,n2
    7 aux(k)=0.d0
      id   =i
      i    =ic
      iw   =0
c     
c     start minimization along next direction
c     
    8 ns   =0
      ip   =0
      db   =0.d0
      if (iw) 10,9,10
    9 hx   =aux(i)
   10 if (it-1) 11,11,14
   11 if (is-1) 14,12,14
   12 dm   =.1d0
      if (dabs(hx)-1.d0) 38,38,13
   13 dm   =-dm*hx
      goto 38
   14 if (is-2) 18,15,18
   15 if (it-1) 17,16,17
   16 dm   =hq
      goto 38
   17 dm   =dq
      goto 38
c     
c     interpolation using estimate of second derivative
c     
   18 if (iw-1) 20,19,20
   19 j    =n2+i
      goto 21
   20 j    =n3+i
   21 hd   =aux(j)
      dc   =1.d-2
      if (it-2) 23,23,22
   22 dc   =hq
   23 dm   =dc
      mk   =1
      goto 51
   24 dm   =dc*hd
      if (dm) 26,25,26
   25 dm   =1.d0
   26 dm   =.5d0*dc-(fm-fb)/dm
      mk   =2
      if (fm-fb) 27,29,29
   27 fc   =fb
      fb   =fm
      db   =dc
      if (dm-db) 28,67,28
   28 dc   =0.d0
      goto 51
   29 if (dm-db) 31,30,31
   30 da   =dc
      fa   =fm
      goto 37
   31 fc   =fm
      goto 51
c     
c     analyse interpolated function value
c     
   32 if (fm-fb) 34,33,33
   33 da   =dm
      fa   =fm
      goto 35
   34 da   =db
      fa   =fb
      db   =dm
      fb   =fm
   35 if ((dc-da)/(db-da)) 36,36,50
   36 if (db) 67,37,67
   37 ns   =1
      dm   =-dc
c     
c     linear search for smaller function values
c     along current direction
c
   38 if (ns-15) 43,43,39
   39 if (fs-fm) 41,40,41
   40 mf   =n+2
      db   =0.d0
      goto 67
   41 if (dabs(dm)-1.d6) 43,43,42
   42 ier  =100
      goto 67
   43 ns   =ns+1
      mk   =3
      goto 51
   44 if (fm-fb) 45,46,47
   45 da   =db
      fa   =fb
      db   =dm
      fb   =fm
      dm   =dm+dm
      goto 38
   46 if (fs-fb) 47,45,47
   47 if (ns-1) 48,48,49
   48 da   =dm
      fa   =fm
      dm   =-dm
      goto 38
   49 dc   =dm
      fc   =fm
c
c     refine minimum using quadratic interpolation
c
   50 hd   =(fc-fb)/(dc-db)+(fa-fb)/(db-da)
      dm   =.5d0*(da+dc)+(fa-fc)/(hd+hd)
      ip   =ip+1
      mk   =4
c
c     step argument vector and calculate function value
c
   51 if (iw-1) 54,52,54
   52 do 53 k=1,n
         l    =m+k
   53    aux(k)=x(k)+dm*aux(l)
      goto 55
   54 aux(i)=hx+dm
   55 fm   =f(aux)
      goto (24,32,44,56),mk
c
c     analyse interpolated function value
c
   56 if (fm-fb) 61,61,57
   57 if (ip-3) 58,62,62
   58 if ((dc-db)/(dm-db)) 60,60,59
   59 dc   =dm
      fc   =fm
      goto 50
   60 da   =dm
      fa   =fm
      goto 50
   61 db   =dm
      fb   =fm
c
c     calculate new estimate of second derivative
c     along the current direction
c
   62 hd=(hd+hd)/(dc-da)
      if (iw-1) 64,63,64
   63 j    =n2+i
      goto 65
   64 j    =n3+i
   65 aux(j)=hd
      if (fb-fs) 67,66,67
   66 db   =0.d0
c
c     save argument vector with smallest function value found
c
   67 if (iw-1) 70,68,70
   68 do 69 k=1,n
         l    =m+k
         j    =n+k
         hd   =db*aux(l)
         aux(j)=aux(j)+hd
         hd   =x(k)+hd
         aux(k)=hd
   69    x(k) =hd
      goto 71
   70 j    =n+i
      aux(j)=aux(j)+db
      hd   =hx+db
      aux(i)=hd
      x(i) =hd
   71 if (ier-100) 72,108,72
c
c     determine direction for next linear search
c
   72 fs   =fb
      if (i-n) 74,73,73
   73 i    =0
   74 i    =i+1
      if (is) 75,75,80
   75 if (db) 77,76,77
   76 if (i-ic) 8,77,8
   77 ic   =i
      is   =1
      if (it-n) 79,79,78
   78 iw   =1
   79 i    =id
      goto 8
   80 m    =m+n
      if (m-n5) 82,81,81
   81 m    =n4
   82 if (is-1) 83,83,94
   83 if (i-1) 84,84,85
   84 iw   =1
   85 if (i-id) 8,86,8
   86 hq=0.d0
      do 87 k=n1,n2
   87    hq=hq+aux(k)*aux(k)
      if (hq) 90,88,90
   88 if (mf-n1) 108,108,89
   89 ier  =200
      goto 108
   90 dq   =dsqrt(hq)
      hq   =dq
      if (hq-1.d0) 92,92,91
   91 hq   =1.d0
   92 do 93 k=n1,n2
         l    =m+k-n
   93    aux(l)=aux(k)/dq
      is   =2
      goto 8
c
c     end of iteration cycle
c     test for termination of minimization
c
   94 is   =0
      tol  =eps
      if (dabs(fs)-1.d0) 96,96,95
   95 tol  =eps*dabs(fs)
   96 if (fl-fs-tol) 100,100,97
   97 if (it-lim) 99,99,98
   98 ier  =10
      goto 108
   99 mf   =0
      goto 6
  100 if (mf-n1) 102,102,101
  101 ier  =200
      goto 108
  102 mf   =mf+1
      dq   =0.d0
      do 105 k=1,n
         j    =n+k
         if (dabs(aux(k))-1.d0) 103,103,104
  103    dq   =dq+dabs(aux(j))
         goto 105
  104    dq   =dq+dabs(aux(j)/aux(k))
  105    continue
      if (dq-eta) 108,108,106
  106 if (mf-n) 6,6,107
  107 ier  =1
  108 y    =fb
      if (ier) 111,111,109
  109 if (isw+12345) 110,111,110
c 110 call wier(ier,20212)
  110 continue
  111 return
      end
