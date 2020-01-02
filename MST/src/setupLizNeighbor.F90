subroutine setupLizNeighbor(print_level)
   use KindParamModule, only : IntKind, RealKind
!
   use MathParamModule, only : ZERO, ONE, TWO, ten2m8, ten2m6
!
   use GroupCommModule, only : getGroupID, getMyPEinGroup, GlobalMaxInGroup
!
   use ErrorHandlerModule, only : WarningHandler, StopHandler, ErrorHandler
!
   use SortModule, only : HeapSort
!
   use NeighborModule, only : setNeighbor, printNeighbor
   use NeighborModule, only : getNumReceives, setMaxReceives
   use NeighborModule, only : printCommunicationTable
!
   use SystemModule, only : getBravaisLattice, getLmaxKKR, setSiteLIZ
   use SystemModule, only : getNumAtoms, getAtomPosition, getAtomicNumber, getScalingFactor
!
   use Atom2ProcModule, only : getLocalNumAtoms, getGlobalIndex
   use Atom2ProcModule, only : getLocalIndex, getAtom2ProcInGroup
!
   use PublicTypeDefinitionsModule, only : LizLmaxStruct, NeighborStruct
!
   use AtomModule, only : getLocalAtomPosition
   use AtomModule, only : getLizLmax
!
   implicit none
!
!   logical, intent(in) :: isScreenTau
   logical :: foundRepeat
!
   integer (kind=IntKind), parameter :: MaxLizAtoms = 2048
   integer (kind=IntKind), parameter :: MaxLizShells = 600
!
   integer (kind=IntKind), intent(in) :: print_level(*)
!
   integer (kind=IntKind) :: GlobalNumAtoms, MyPEinGroup, GroupID
   integer (kind=IntKind) :: LocalNumAtoms, na, ns, ms
   integer (kind=IntKind) :: index_mn_max, imn
   integer (kind=IntKind) :: i, j, k, m, id, js, ja, ig
   integer (kind=IntKind) :: MaxReceives
   integer (kind=IntKind) :: LizAtoms_max
   integer (kind=IntKind) :: aid(MaxLizAtoms)
   integer (kind=IntKind) :: ash(MaxLizAtoms)
   integer (kind=IntKind) :: lmax_ja, lmax_sh
   integer (kind=IntKind) :: indx(MaxLizShells)
   integer (kind=IntKind) :: nas(MaxLizShells)
   integer (kind=IntKind) :: n1min, n1max, n2min, n2max, n3min, n3max, num_ks
   integer (kind=IntKind) :: vecCount, vecCount_sav
   integer (kind=IntKind) :: maxIter, count
   integer (kind=IntKind), allocatable :: m1(:), m2(:), m3(:)
   integer (kind=IntKind), allocatable :: index_mn(:)
!
   real (kind=RealKind) :: posl(3), posg(3)
   real (kind=RealKind) :: x0, y0, z0, x, y, z, rcut, rcut_max
   real (kind=RealKind) :: r, r1, r2, r3, bk1, bk2, bk3, p, scaling
   real (kind=RealKind) :: rs(MaxLizShells)
   real (kind=RealKind) :: brav(3,3)
   real (kind=RealKind) :: apos(3,MaxLizAtoms)
   real (kind=RealKind), allocatable :: apos_mn(:,:,:)
!
   type (NeighborStruct) :: neighb
   type (LizLmaxStruct), pointer :: LizLmax
!
   brav(1:3,1:3) = getBravaisLattice()
   r1 = sqrt( brav(1,1)*brav(1,1) + brav(2,1)*brav(2,1) + brav(3,1)*brav(3,1) )
   r2 = sqrt( brav(1,2)*brav(1,2) + brav(2,2)*brav(2,2) + brav(3,2)*brav(3,2) )
   r3 = sqrt( brav(1,3)*brav(1,3) + brav(2,3)*brav(2,3) + brav(3,3)*brav(3,3) )
!
   GlobalNumAtoms = getNumAtoms()
   GroupID = getGroupID('Unit Cell')
   MyPEinGroup = getMyPEinGroup(GroupID)
   LocalNumAtoms = getLocalNumAtoms()
   scaling = getScalingFactor()
!
   maxIter = 10000
!
   allocate( index_mn(GlobalNumAtoms), apos_mn(3,MaxLizAtoms,GlobalNumAtoms) )
   MaxReceives = 0
   LOOP_id: do id = 1, LocalNumAtoms
      LizLmax => getLizLmax(id)
!
      neighb%NumAtoms = 0
      neighb%NumReceives = 0
!
      ig = getGlobalIndex(id)
      posl(1:3) = getLocalAtomPosition(id)
      rcut_max = LizLmax%rad
      rcut = TEN2m6
      LizAtoms_max = LizLmax%nmax
!
!     ================================================================
!     determine the number of shifts of the unit cell along brav(:,1),
!     brav(:,2), and brav(:,3) directions.
!     ================================================================
      p = abs( posl(1)*brav(1,1)+posl(2)*brav(2,1)+posl(3)*brav(3,1) )/r1
      n1min = -ceiling((rcut_max+p)/r1)
      n1max = -n1min               ! floor((LizLmax%rad+p)/r1)
      p = abs( posl(1)*brav(1,2)+posl(2)*brav(2,2)+posl(3)*brav(3,2) )/r2
      n2min = -ceiling((rcut_max+p)/r2)
      n2max = -n2min               ! floor((LizLmax%rad+p)/r2)
      p = abs( posl(1)*brav(1,3)+posl(2)*brav(2,3)+posl(3)*brav(3,3) )/r3
      n3min = -ceiling((rcut_max+p)/r3)
      n3max = -n3min               ! floor((LizLmax%rad+p)/r3)
!     ================================================================
!     Check to make sure that there are enough cell shifts (aka repeats)
!     to catch all of the neighbors closer than rcut_max
!     ================================================================
      vecCount = 0
      vecCount_sav = 0
      count = 0
      foundRepeat = .false.
      do while(.not.foundRepeat.and.count.lt.maxIter)
         do i=n1min,n1max
            do j=n2min,n2max
               do k=n3min,n3max
                  do m = 1, GlobalNumAtoms
                     posg(1:3) = getAtomPosition(m)
                     x0 = posg(1) - posl(1)
                     y0 = posg(2) - posl(2)
                     z0 = posg(3) - posl(3)
                     bk1  = brav(1,1)*i + brav(1,2)*j + brav(1,3)*k
                     bk2  = brav(2,1)*i + brav(2,2)*j + brav(2,3)*k
                     bk3  = brav(3,1)*i + brav(3,2)*j + brav(3,3)*k
                     x    = x0 + bk1
                     y    = y0 + bk2
                     z    = z0 + bk3
                     r = sqrt( x*x + y*y + z*z )
                     if(r <= rcut_max+ten2m8) then
                        vecCount = vecCount+1
                     endif
                  enddo
               enddo
            enddo
         enddo
         if(vecCount.eq.vecCount_sav) then
            foundRepeat = .true.
         else
            vecCount_sav = vecCount
            vecCount = 0
            n1min = n1min-1
            n2min = n2min-1
            n3min = n3min-1
            n1max = n1max+1
            n2max = n2max+1
            n3max = n3max+1
         endif
         count = count+1
      enddo
      num_ks = (n3max-n3min+1)*(n2max-n2min+1)*(n1max-n1min+1)
!
!     n1max = ceiling(LizLmax%rad/r1)
!     n1min = -n1max
!     n2max = ceiling(LizLmax%rad/r2)
!     n2min = -n2max
!     n3max = ceiling(LizLmax%rad/r3)
!     n3min = -n3max
!     num_ks = (n3max-n3min+1)*(n2max-n2min+1)*(n1max-n1min+1)
!
      if (print_level(id) >= 0) then
         write(6,*)"setupLIZ :: Number of cell shifts:",num_ks
         write(6,*)"                            nXmax:", n1max, n2max, n3max
         write(6,*)"                            nXmin:", n1min, n2min, n3min
      endif
!
      allocate( m1(num_ks), m2(num_ks), m3(num_ks) )
      ns = 0
      do i=n3min, n3max
         do j=n2min, n2max
            do k=n1min, n1max
               ns = ns + 1
               m1(ns) = k
               m2(ns) = j
               m3(ns) = i
            enddo
         enddo
      enddo
!
      if ( LizAtoms_max>0 ) then
         ns = 0
         nas = 0
         do i = 1, GlobalNumAtoms
            posg(1:3) = getAtomPosition(i)
            x0 = posg(1) - posl(1)
            y0 = posg(2) - posl(2)
            z0 = posg(3) - posl(3)
            LOOP_k0: do k = 1,num_ks
               bk1  = brav(1,1)*m1(k) + brav(1,2)*m2(k) + brav(1,3)*m3(k)
               bk2  = brav(2,1)*m1(k) + brav(2,2)*m2(k) + brav(2,3)*m3(k)
               bk3  = brav(3,1)*m1(k) + brav(3,2)*m2(k) + brav(3,3)*m3(k)
               x    = x0 + bk1
               y    = y0 + bk2
               z    = z0 + bk3
               r = sqrt( x*x + y*y + z*z )
               if (r < ONE .and. ig == i) then
                  cycle LOOP_k0
!              else if (r < ONE) then
!                 -------------------------------------------------
!                 call ErrorHandler('setupLizNeighbor','atoms are too close')
!                 -------------------------------------------------
!                 write(6,'(50(''=''))')
!                 write(6,'(a,f15.8)')'The following pair of atoms has small distance:',r/scaling
!                 write(6,'(a)')'Index         x                y                z'
!                 write(6,'(i5,3(2x,f15.8))')ig,posl(1:3)/scaling
!                 write(6,'(i5,3(2x,f15.8))')i,posg(1:3)/scaling
!                 if (m1(k) /= 0 .or. m2(k) /= 0 .or. m3(k) /= 0) then
!                    write(6,'(a,3i5,x,a)')'After the cell shifted by:',m1(k),m2(k),m3(k), &
!                                          'the following atom is at'
!                    write(6,'(i5,3(2x,f15.8))')i,(x+posl(1))/scaling,(y+posl(2))/scaling, &
!                                                 (z+posl(3))/scaling
!                 endif
!                 write(6,'(50(''-''))')
               else if ( r <= rcut_max+ten2m8 ) then
                  do js = 1, ns
                     if (abs(r - rs(js)) < TEN2m6) then
                        nas(js) = nas(js)+1
                        cycle LOOP_k0
                     endif
                  enddo
                  ns = ns + 1 ! count no. of neighboring shells within LizLmax%rad
                  if (ns > MaxLizShells) then
!                    -------------------------------------------------
                     call StopHandler('setupLizNeighbor','ns > MaxLizShells')
!                    ----------------------------------------------
                  endif
                  rs(ns) = r
                  nas(ns) = nas(ns)+1
               endif
            enddo LOOP_k0
         enddo
!        -------------------------------------------------------------
         call HeapSort(ns, rs, indx)
!        -------------------------------------------------------------
         na = 0
         Loop_js: do js = 1,ns
            na = na+nas(indx(js))
            if ( na <= LizAtoms_max ) then
               rcut = rs(js)
            else
               rcut = rcut+0.001d0
               exit Loop_js
            endif
         enddo Loop_js
      endif
!     ================================================================
!
      index_mn = 0
      index_mn_max = 0
      ns = 0
      na = 0
      do i = 1, GlobalNumAtoms
         posg(1:3) = getAtomPosition(i)
!        print *,num_ks,i,ig,"posL,posG::",posl,posg
         x0 = posg(1) - posl(1)
         y0 = posg(2) - posl(2)
         z0 = posg(3) - posl(3)
         LOOP_k: do k = 1,num_ks
            bk1  = brav(1,1)*m1(k) + brav(1,2)*m2(k) + brav(1,3)*m3(k)
            bk2  = brav(2,1)*m1(k) + brav(2,2)*m2(k) + brav(2,3)*m3(k)
            bk3  = brav(3,1)*m1(k) + brav(3,2)*m2(k) + brav(3,3)*m3(k)
            x    = x0 + bk1
            y    = y0 + bk2
            z    = z0 + bk3
            r = sqrt( x*x + y*y + z*z )
            if ( r < (LizLmax%rad_s + ten2m8) ) then
               index_mn(i) = index_mn(i) +1     ! count no. of i atoms
!                                               ! contributing to tau
               if ( index_mn(i) > MaxLizAtoms ) then
!                 ----------------------------------------------------
                  call StopHandler( 'setupLizNeighbor',               &
                                    'na > MaxLizAtoms' )
!                 ----------------------------------------------------
               endif
               apos_mn(1,index_mn(i),i) = x
               apos_mn(2,index_mn(i),i) = y
               apos_mn(3,index_mn(i),i) = z
            endif
            if (r < ONE .and. ig == i) then
               cycle LOOP_k
            else if (r < ONE) then
               write(6,'(50(''=''))')
               write(6,'(a,f15.8)')'The following pair of atoms has small distance:',r/scaling
               write(6,'(a)')'Index         x                y                z'
               write(6,'(i5,3(2x,f15.8))')ig,posl(1:3)/scaling
               write(6,'(i5,3(2x,f15.8))')i,posg(1:3)/scaling
               if (m1(k) /= 0 .or. m2(k) /= 0 .or. m3(k) /= 0) then
                  write(6,'(a,3i5,x,a)')'After the cell shifted by:',m1(k),m2(k),m3(k), &
                                        'the following atom is at'
                  write(6,'(i5,3(2x,f15.8))')i,(x+posl(1))/scaling,(y+posl(2))/scaling, &
                                               (z+posl(3))/scaling
               endif
               write(6,'(50(''-''))')
            else if (r <= rcut+ten2m8) then
               na = na + 1 ! count no. of atoms within LizLmax%rad
               if (na > MaxLizAtoms) then
!                 ----------------------------------------------------
                  call StopHandler('setupLizNeighbor','na > MaxLizAtoms')
!                 ----------------------------------------------------
               endif
               aid(na) = i ! stores global index of neighboring atom
               apos(1,na) = posg(1) + bk1
               apos(2,na) = posg(2) + bk2
               apos(3,na) = posg(3) + bk3
               do js = 1, ns
                  if (abs(r - rs(js)) < 0.01) then
                     ash(na) = js
                     cycle LOOP_k
                  endif
               enddo
               ns = ns + 1 ! count no. of neighboring shells within LizLmax%rad
               if (ns > MaxLizShells) then
!                 ----------------------------------------------------
                  call StopHandler('setupLizNeighbor','ns > MaxLizShells')
!                 ----------------------------------------------------
               endif
               rs(ns) = r
               ash(na) = ns ! stores shell index of neighboring atom
            endif
         enddo LOOP_k
         index_mn_max = max(index_mn_max,index_mn(i))
      enddo
!
!     if (ns == 0) then
!        -------------------------------------------------------------
!        call WarningHandler('setupLizNeighbor','no neighbor for this atom',ig)
!        -------------------------------------------------------------
!     else
      if (ns > 0) then
!        -------------------------------------------------------------
         call HeapSort(ns, rs, indx)
!        -------------------------------------------------------------
         neighb%NumAtoms = na
         neighb%NumShells = ns
         allocate( neighb%Z(na), neighb%Lmax(na) )
         allocate( neighb%ProcIndex(na), neighb%LocalIndex(na) )
         allocate( neighb%GlobalIndex(na), neighb%Position(3,na) )
         allocate( neighb%IndexMN(GlobalNumAtoms) )
         allocate( neighb%Rmunu(3,index_mn_max,GlobalNumAtoms))
         allocate( neighb%ShellIndex(na) ) 
         allocate( neighb%ShellRad(ns) )
         do i=1,ns
            neighb%ShellRad(i)=rs(i)
         enddo
         neighb%IndexMN = -1
         neighb%Rmunu = real(-10000)
         ms = LizLmax%NumShells
         do i = 1,GlobalNumAtoms
            neighb%IndexMN(i) = index_mn(i)
            do imn = 1,index_mn(i)
               neighb%Rmunu(1:3,imn,i) = apos_mn(1:3,imn,i)
            enddo
         enddo
         neighb%NumReceives = 0

         do i = 1, na
            ja = aid(i)
            neighb%GlobalIndex(i) = ja
            neighb%Z(i) = getAtomicNumber(ja)
            neighb%Position(1:3,i) = apos(1:3,i)
            neighb%ProcIndex(i) = getAtom2ProcInGroup(ja)
            neighb%LocalIndex(i) = getLocalIndex(ja)
            neighb%Lmax(i) = -1
            neighb%ShellIndex(i)=0
!
            do js = 1, ns
               if (ash(i) == indx(js)) then
                  neighb%ShellIndex(i)=js
                  lmax_ja = getLmaxKKR(ja)
                  if (js <= ms) then
                     lmax_sh = LizLmax%lmax_shell(js)
!                    if ( lmax_sh > lmax_ja .and. lmax_ja/=0 ) then
                     if ( lmax_sh > lmax_ja .and. lmax_ja>=0 ) then
                        neighb%Lmax(i) = lmax_ja
                     else
                        neighb%Lmax(i) = lmax_sh
                     endif
                  else
                     lmax_sh = LizLmax%lmax_shell(ms)
!                    if ( lmax_sh > lmax_ja .and. lmax_ja/=0 ) then
                     if ( lmax_sh > lmax_ja .and. lmax_ja>=0 ) then
                        neighb%Lmax(i) = lmax_ja
                     else
                        neighb%Lmax(i) = lmax_sh
                     endif
                  endif
               endif
            enddo
            if (neighb%Lmax(i) < 0) then
!               ------------------------------------------------------
                call WarningHandler('setupLizNeighbor','neighb%Lmax < 0')
!               ------------------------------------------------------
            endif
            if (neighb%ProcIndex(i) /= MyPEinGroup) then
               neighb%NumReceives = neighb%NumReceives + 1
            endif
         enddo
!        -------------------------------------------------------------
         call setNeighbor(id,neighb)
!        -------------------------------------------------------------
         deallocate( neighb%Z, neighb%Lmax, neighb%ProcIndex,          &
                     neighb%LocalIndex, neighb%GlobalIndex )
         deallocate( neighb%Position, neighb%IndexMN, neighb%Rmunu )
         deallocate( neighb%ShellIndex, neighb%ShellRad)
      endif
      deallocate( m1, m2, m3 )
!
      if (print_level(id) >= 0) then
!        -------------------------------------------------------------
         call printNeighbor(id)
!        -------------------------------------------------------------
      endif
      call setSiteLIZ(ig,neighb%NumAtoms)
      MaxReceives = max(MaxReceives, getNumReceives(id))
   enddo LOOP_id
   deallocate( index_mn, apos_mn )
!
!  -------------------------------------------------------------------
   call GlobalMaxInGroup(GroupID,MaxReceives)
   call setMaxReceives(MaxReceives)
!  -------------------------------------------------------------------
   call buildSendTable(print_level)
!  -------------------------------------------------------------------
   if (maxval(print_level(1:LocalNumAtoms)) >= 0) then
      call printCommunicationTable()
   endif
!
end subroutine setupLizNeighbor
