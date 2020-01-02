c
c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine lattice(vbra,vcut,nm1,nm2,nm3,
     >                   vlat_1,vlat_2,vlat_3,vsq,nv,ipmax,
     >                   iprint,istop)
c     ================================================================
c
      implicit   none
c
      character  sname*32
      character  istop*32
c
      integer    nv
      integer    ipmax
      integer    iprint
      integer    nm1
      integer    nm2
      integer    nm3
      integer    tn1p1
      integer    tn2p1
      integer    tn3p1
      integer    nt12
      integer    nt23
      integer    nt123
      integer    n1
      integer    n2
      integer    n3
      integer    i
c
      real*8     vlat_1(ipmax)
      real*8     vlat_2(ipmax)
      real*8     vlat_3(ipmax)
      real*8     vsq(ipmax)
      real*8     vn(3)
      real*8     vbra(3,3)
      real*8     vnsq
      real*8     vcut
      real*8     vcut2
      real*8     tol
      real*8     zero
c
      parameter (tol=1.0d-06)
      parameter (zero=0.0d0)
      parameter (sname='lattice')
c
c     ****************************************************************
c     Calculates lattice vectors .....................................
c     input:
c       vbra     : basis vectors for bravais lattice
c       vcut     : cut off radius for lattice vectors
c       nm1      : number of repeats of vbra(1)
c       nm2      : number of repeats of vbra(2)
c       nm3      : number of repeats of vbra(3)
c       ipmax    : dimension of array that will contain vectors
c     output:
c       vlat     : bravais lattice vectors
c       vsq      : square of length of lattice vectors
c       nv       : number of lattice vectors
c
c     ****************************************************************
c
c     ================================================================
      if(iprint.ge.1) then
         write(6,'(''  lattice:: nm1,nm2,nm3 = '',3i3)') nm1,nm2,nm3
         write(6,'(''  lattice:: Basis vectors'')')
         do n1=1,3
            write(6,'(5x,3f10.5)')  (vbra(n2,n1),n2=1,3)
         enddo
      endif
c     ================================================================
c     generate lattice vectors........................................
      nv=0
      vcut2=vcut*vcut+tol
      tn1p1=2*nm1+1
      tn2p1=2*nm2+1
      tn3p1=2*nm3+1
      nt12=tn1p1*tn2p1
      nt23=tn2p1*tn3p1
      nt123=nt12*tn3p1
      do i=1,nt123
	 n1=i-1
         n1=mod(n1,tn1p1)-nm1
         n2=(i-1)/tn1p1
         n2=mod(n2,tn2p1)-nm2
         n3=(i-1)/nt12
         n3=mod(n3,tn3p1)-nm3
         vn(1) = n1*vbra(1,1)+n2*vbra(1,2)+n3*vbra(1,3)
         vn(2) = n1*vbra(2,1)+n2*vbra(2,2)+n3*vbra(2,3)
         vn(3) = n1*vbra(3,1)+n2*vbra(3,2)+n3*vbra(3,3)
         vnsq=vn(1)*vn(1)+vn(2)*vn(2)+vn(3)*vn(3)
         if(vnsq.le.vcut2) then
            if(nv+1.gt.ipmax) then
               write(6,'('' lattice:: nv.gt.ipmax:'',2i5)')
     >                                nv,ipmax
               call fstop(sname)
            endif
c           ----------------------------------------------------------
            call ord3v(vlat_1,vlat_2,vlat_3,vsq,nv,vn,vnsq,istop)
c           ----------------------------------------------------------
         endif
      enddo
c     ================================================================
      if(iprint.ge.1) then
         write(6,'(/'' lattice: number of lattice vectors is'',
     >                i5)') nv
c        write(6,'('' n,vlat '',i5,4f10.5)')
c    >            (n1,vlat_1(n1),vlat_2(n1),vlat_3(n1),vsq(n1),n1=1,nv)
      endif
c
c     ================================================================
      if(sname.eq.istop) then
         call fstop(sname)
      else
         return
      endif
c
      end
