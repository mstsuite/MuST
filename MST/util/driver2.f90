program driver2
   use KindParamModule, only : IntKind, RealKind
   use MathParamModule, only : PI, THREE, THIRD, PI4
   use MPPModule, only : initMPP, endMPP
   use MPPModule, only : syncAllPEs
   use MPPModule, only : MyPE, NumPEs, AllPEs, AnyPE
   use MPPModule, only : bcastMessage, sendMessage, recvMessage
   use MPPModule, only : packMessage, unpackMessage, sendPackage, recvPackage
   use PolyhedraModule, only : initPolyhedra, genPolyhedron, endPolyhedra
   use PolyhedraModule, only : getNeighborIndex, getNumPlanes, getPlaneArea
   use PolyhedraModule, only : getInscrSphRadius, getOutscrSphRadius, getVolume
   use PolyhedraModule, only : getWignerSeitzRadius, getBoxVolume
   use PolyhedraModule, only : printPolyhedron, getNeighborDistance
   implicit none
!
   character (len=70) :: text
   character (len=70) :: path
   character (len=30) :: info_table
!
   integer (kind=IntKind) :: num_atoms
   integer (kind=IntKind) :: i, j, k, lmax, n
   integer (kind=IntKind) :: mtasa
   integer (kind=IntKind) :: ngaussr, ngaussq, np
   integer (kind=IntKind), allocatable :: alist(:)
   integer (kind=IntKind), allocatable :: anum(:)
!
   real (kind=RealKind) :: bravais(3,3)
   real (kind=RealKind) :: rins, rins0
   real (kind=RealKind) :: rcut, r1, r2, r3, r4, rmt, rws, rcirc
   real (kind=RealKind), allocatable :: position(:,:)
   real (kind=RealKind), allocatable :: rad(:)
   real (kind=RealKind), allocatable :: area(:), dist(:)
   real (kind=RealKind) :: omegmt, omegws, omegbra, v
!
   call initMPP()
!
   if (MyPE == 0) then
      do i=1,6
         read(5,'(a)')text
      enddo
      read(5,'(a)')path
      do i=1,13
         read(5,'(a)')text
      enddo
      read(5,*)num_atoms
      read(5,'(a)')text
      read(5,*)i,j,mtasa,rins
      do i=1,5
         read(5,'(a)')text
      enddo
      read(5,*)ngaussr, ngaussq
      read(5,'(a)')text
      read(5,'(a)')info_table
!     ----------------------------------------------------------------
      call packMessage(num_atoms)
      call packMessage(mtasa)
      call packMessage(ngaussr)
      call packMessage(ngaussq)
      call packMessage(rins)
      call sendPackage(1000,AllPEs)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call recvPackage(1000,0)
      call unpackMessage(num_atoms)
      call unpackMessage(mtasa)
      call unpackMessage(ngaussr)
      call unpackMessage(ngaussq)
      call unpackMessage(rins)
!     ----------------------------------------------------------------
   endif
!
   allocate( position(3,num_atoms) )
   allocate( rad(num_atoms) )
!
   if (MyPE == 0) then
      allocate( anum(num_atoms) )
      open(unit=10,file=trim(path)//info_table,status='old')
      read(10,'(a)')text
      read(10,'(a)')text
      read(10,'(a)')text
      read(10,*)bravais(1:3,1)
      read(10,*)bravais(1:3,2)
      read(10,*)bravais(1:3,3)
!     ----------------------------------------------------------------
      call bcastMessage(bravais,3,3,0)
!     ----------------------------------------------------------------
      do i=1,5
         read(10,'(a)')text
      enddo
      do i=1,num_atoms
         read(10,*)anum(i),k,position(1,i),position(2,i),position(3,i),lmax, &
                   rcut,r1,r2,r3,r4,rad(i)
      enddo
!     ----------------------------------------------------------------
      call bcastMessage(position,3,num_atoms,0)
      call bcastMessage(rad,num_atoms,0)
!     ----------------------------------------------------------------
      close(unit=10)
   else
!     ----------------------------------------------------------------
      call bcastMessage(bravais,3,3,0)
      call bcastMessage(position,3,num_atoms,0)
      call bcastMessage(rad,num_atoms,0)
!     ----------------------------------------------------------------
   endif
!
   n = 0
   do i=1,num_atoms
      if (mod(i-1,NumPEs) == MyPE) then
         n = n+1
      endif
   enddo
!
   call initPolyhedra(n,bravais,'none',-1)
!
   if (MyPE == 0) then
      open(unit=11,file=trim(path)//'voronoi.dat',status='unknown')
      open(unit=12,file=trim(path)//'neighbor.dat',status='unknown')
      write(11,'(80(''=''))')
      write(11,'(a)')            &
      ' Atom    No. Surf     Vcell       Vmt         Rmt       Rcirc       Rws'
      write(11,'(80(''-''))')
   endif
!
!  -------------------------------------------------------------------
   call syncAllPEs()
!  -------------------------------------------------------------------
!
   rins0 = rins
!  -------------------------------------------------------------------
   omegbra = getBoxVolume()
!  -------------------------------------------------------------------
   v = 0.0d0
   n = 0
   do i=1,num_atoms
      if (mod(i-1,NumPEs) == MyPE) then
         rins = rins0
         n = n+1
!        -------------------------------------------------------------
         call genPolyhedron(n,i,num_atoms,position,rad)
!        -------------------------------------------------------------
         np = getNumPlanes(n)
         rmt = getInscrSphRadius(n)
         rws = getWignerSeitzRadius(n)
         rcirc = getOutscrSphRadius(n)
         omegws = getVolume(n)
!        -------------------------------------------------------------
         allocate(alist(np), area(np), dist(np))
         do j=1,np
            alist(j) = getNeighborIndex(n,j)
            area(j)  = getPlaneArea(n,j)
            dist(j)  = getNeighborDistance(n,j)
         enddo
         omegmt = PI4*THIRD*rmt**3
!
         if (mtasa >= 1) then
            if(rins > 0.0d0) then
!              if rins is preset in the input then force all rws the same
               rws=(THREE*omegbra/(PI4*num_atoms))**THIRD
               omegws=PI4*rws**3/THREE
            else
!              otherwise rins is the real rmt
               rins=rmt
            endif
            rmt = rws
            omegmt = omegws
            if (mtasa == 1) then
               rins=rws
            endif
         else if (rins <= 0.0d0) then
!           otherwise rins is the real rmt
            rins=rmt
         endif
!
         if (MyPE /= 0) then
            call packMessage(omegws)
            call packMessage(omegmt)
            call packMessage(rmt)
            call packMessage(rcirc)
            call packMessage(rws)
            call sendPackage(3000+i,0)
            call sendMessage(np,4000+i,0)
            call sendMessage(alist,np,5000+i,0)
            call sendMessage(area,np,6000+i,0)
            call sendMessage(dist,np,7000+i,0)
         endif
      else if (MyPE == 0) then
         j = mod(i-1,NumPEs)
         call recvPackage(3000+i,j)
         call unpackMessage(omegws)
         call unpackMessage(omegmt)
         call unpackMessage(rmt)
         call unpackMessage(rcirc)
         call unpackMessage(rws)
         call recvMessage(np,4000+i,j)
         allocate( alist(np), area(np), dist(np) )
         call recvMessage(alist,np,5000+i,j)
         call recvMessage(area,np,6000+i,j)
         call recvMessage(dist,np,7000+i,j)
      endif
      if (MyPE == 0) then
         v = v + omegws
         write(11,'(i5,4x,i5,4x,5(1x,f10.5))')i,np,omegws,omegmt,rmt,rcirc,rws
         write(12,'(3i5))')i,anum(i),np
         write(12,'(3(i5,2f10.5,a))')(alist(j),area(j),dist(j),';',j=1,np)
      endif
      if (allocated(alist)) then
         deallocate( alist )
      endif
      if (allocated(area)) then
         deallocate( area )
      endif
      if (allocated(dist)) then
         deallocate( dist )
      endif
   enddo
   if (MyPE == 0) then
      close(unit=11)
      close(unit=12)
      deallocate( anum )
      print *,'v = ',v
      print *,'omegbra = ',omegbra
      print *,'Volume difference = ',abs(v-omegbra)
   endif
!
   deallocate( position, rad )
   call endPolyhedra()
   call endMPP()
end program driver2
