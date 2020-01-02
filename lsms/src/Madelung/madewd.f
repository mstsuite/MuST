c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine madewd(lastrs,lastkn,myatom,num_atoms,alat,
     >                  rslat_x,rslat_y,rslat_z,rslatmd,
     >                  knlat_x,knlat_y,knlat_z,knlatsq,
     >                  atom_posi_x,atom_posi_y,atom_posi_z,
     >                  eta,tau,madmat,
     >                  iprint,istop)
c     ================================================================
c
c     ****************************************************************
c     program calculates the madelung constants
c     written by w.a.s.,jr. 2-21-88
c     modified by y.w., 06-20-94
c     modified by y.w., 06-06-95
c
c     ****************************************************************
c     input :
c     =====
c         1.) ewald card
c             lastrs = number real-space vectors
c             rslat= real-space vectors (in the units of a0=1)
c             rslatmd= rslat module
c             lastkn = number reciprocal-space vectors
c             knlat = reciprocal-space vectors (in the units of a0=1)
c             knlatsq= knlat square
c             eta = ewald parameter
c             num_atoms = number of atoms in the unit cell
c             tau = Wigner-Seitz volume
c
c     calling sequence
c     ================
c       |
c     intewld                      ! evaluates integral over eta
c       |
c       |
c     madsum                       ! real & recip space sums
c
c     ****************************************************************
c
      implicit   none
c
      character  istop*32
      character  sname*32
      parameter (sname='madewd')
c
      integer    lastrs
      integer    lastkn
      integer    num_atoms
      integer    iprint
      integer    ibegin
      integer    myatom
      integer    n
      integer    i
c
      real*8     atom_posi_x(num_atoms)
      real*8     atom_posi_y(num_atoms)
      real*8     atom_posi_z(num_atoms)
      real*8     rslat_x(lastrs)
      real*8     rslat_y(lastrs)
      real*8     rslat_z(lastrs)
      real*8     rslatmd(lastrs)
      real*8     knlat_x(lastkn)
      real*8     knlat_y(lastkn)
      real*8     knlat_z(lastkn)
      real*8     knlatsq(lastkn)
      real*8     madmat(num_atoms)
      real*8     alat
      real*8     aij(3)
      real*8     r0tm
      real*8     term0
      real*8     tau
      real*8     eta
      real*8     fnpi
      real*8     pi
      real*8     pi4
      real*8     zero
      real*8     two
      real*8     four
c
      parameter (zero=0.0d0)
      parameter (two=2.0d0)
      parameter (four=4.0d0)
c
c     ****************************************************************
c     calculate Madelung's constant matrix:
c
c     for i <> j,
c                                   2   2       -> ->
c            4*pi          1    -eta *Kq /4 - i*Kq*aij
c     M   =  ---- * sum  ----- e
c      ij    tau    q<>0    2
c                         Kq
c
c                          ->   ->                  2
c                 1 - erf(|Rn + aij|/eta)     pi*eta
c          + sum ------------------------- - ---------
c             n         ->   ->                 tau
c                      |Rn + aij|
c
c     for i = j,
c                                   2   2
c            4*pi          1    -eta *Kq /4
c     M   =  ---- * sum  ----- e
c      ii    tau    q<>0    2
c                         Kq
c
c                                             2       
c                  1 - erf(Rn/eta)      pi*eta           2
c          + sum  ----------------- - ---------- - --------------
c            n<>0         Rn             tau        sqrt(pi)*eta
c
c     eta is the Ewald's parameter;
c     tau, atom_posi_*, rslat_* and knlat_* are in the units of a0=1;
c     madmat is in the units of a0=alat;
c     ****************************************************************
c
      pi=fnpi()
      pi4=four*pi
c
      term0=-pi*eta*eta/tau
c
c     ================================================================
c     start madmat calculation
c     ================================================================
      do n =1,num_atoms
c        =============================================================
c        aij is in the units of a0 = 1................................
c        =============================================================
         aij(1) = atom_posi_x(myatom) - atom_posi_x(n)
         aij(2) = atom_posi_y(myatom) - atom_posi_y(n)
         aij(3) = atom_posi_z(myatom) - atom_posi_z(n)
c
         if( n .eq. myatom ) then
            ibegin=2
            r0tm=-two/sqrt(pi)/eta
         else
            ibegin=1
            r0tm=zero
         endif
c        =============================================================
c        subrtact aij from rslat and calculate rslatmd which is used in 
c        calculating the real-space integral
c        rslatmd, and aij are in the units of a0 = 1
c        =============================================================
         do i = 1,lastrs
            rslatmd(i)=sqrt( (rslat_x(i)-aij(1))*(rslat_x(i)-aij(1))
     >                      +(rslat_y(i)-aij(2))*(rslat_y(i)-aij(2))
     >                      +(rslat_z(i)-aij(3))*(rslat_z(i)-aij(3)) )
         enddo
c        =============================================================
c        perform the reciprocal-space sum and real-space sum
c        -------------------------------------------------------------
         call madsum(lastkn,lastrs,ibegin,aij,
     >               rslatmd,knlat_x,knlat_y,knlat_z,knlatsq,
     >               eta,tau,madmat(n),
     >               pi4,iprint,istop)
c        -------------------------------------------------------------
         madmat(n)=madmat(n)+r0tm+term0
c        =============================================================
c        The unit of madmat is finally resumed to a0 = alat
c        =============================================================
         madmat(n)=madmat(n)/alat
      enddo                                  ! end do loop over n
c
      if( istop .eq. sname ) then
        call fstop(sname)
      endif
c     ================================================================
c
      return
      end
