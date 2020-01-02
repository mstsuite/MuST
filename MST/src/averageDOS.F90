!!!ywg needs to fix nume  2-12-2010
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine averageDOS(LocalNumAtoms, npola, nume, egrd, dos)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use MathParamModule, only : ZERO, ONE
!
   use MPPModule, only : setCommunicator, sendMessage, recvMessage, getSource, resetCommunicator
!
   use GroupCommModule, only : getGroupID, getMyPEinGroup, getNumPEsInGroup, getGroupCommunicator
!
   use ScfDataModule, only : inputpath
!
   use SystemModule, only : getNumAtomTypes, getAtomType, getNumAtomsOfType
   use SystemModule, only : getAtomTypeName, getSystemID, getSystemTitle
!
   use Atom2ProcModule, only : getLocalNumAtoms, getMaxLocalNumAtoms
   use Atom2ProcModule, only : getGlobalIndex
!
   implicit none
!
   character (len=3) :: atname
   character (len=150) :: fname
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms, npola, nume
   integer (kind=IntKind) :: it, ie, i, j, funit, pe, nla, iglb, ksize
   integer (kind=IntKind) :: LocalMaxAtoms, NumAtomTypes, MessageSize
   integer (kind=IntKind) :: GroupID, MyPEinGroup, NumPEsInGroup, comm
   integer (kind=IntKind), parameter :: OutputNode = 0
   integer (kind=IntKind), parameter :: MessageID = 10010
!
   real (kind=RealKind), intent(in) :: dos(1:nume,1:npola,1:2,1:LocalNumAtoms)
   real (kind=RealKind), allocatable :: msgbuf(:)
   real (kind=RealKind), allocatable :: fspace(:,:,:,:)
   real (kind=RealKind) :: alp
!
   complex (kind=CmplxKind), intent(in) :: egrd(1:nume)
!
   LocalMaxAtoms = getMaxLocalNumAtoms()
   NumAtomTypes = getNumAtomTypes()
   MessageSize = 2*nume*npola*LocalMaxAtoms
   allocate( msgbuf(1:MessageSize) )
   msgbuf(1:MessageSize) = ZERO
!
   GroupID = getGroupID('Unit Cell')
   MyPEinGroup = getMyPEinGroup(GroupID)
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   comm = getGroupCommunicator(GroupID)
   call setCommunicator(comm,MyPEinGroup,NumPEsInGroup,sync=.true.)
!
   if ( MyPEinGroup /= OutputNode ) then
!     ================================================================
!     send the D.O.S. to node_output for calculating the average
!     potential....................................................
!     ----------------------------------------------------------------
      call dcopy(2*nume*npola*LocalNumAtoms,dos,1,msgbuf,1)
!     ----------------------------------------------------------------
      call sendMessage(msgbuf,MessageSize,MessageID,OutputNode)
!     ----------------------------------------------------------------
   else
      allocate( fspace(1:nume,1:npola,2,1:NumAtomTypes) )
      fspace(1:nume,1:npola,2,1:NumAtomTypes) = ZERO
!
      do it = 1, NumAtomTypes
         funit = 30+it
         atname = getAtomTypeName(it)
         fname=trim(inputpath)//trim(atname)//'_averdos_'//getSystemID()
         open(unit=funit,file=trim(fname),status='replace')
!
         write(funit,'(a)') getSystemTitle()
         write(funit,'(''Total spin indeces'',t40,''='',i5)') npola
         write(funit,'(''Num. of energies'',t40,''='',i5)') nume
         write(funit,'(a,a)')'Atom Type = ',atname
         write(funit,'(/,64(''=''))')
         if (npola == 1) then
            write(funit,'(7x,a,16x,a,16x,a)')'energy','DOS-MT', "DOS-WS"
         else
            write(funit,'(7x,a,15x,a,10x,a,15x,a,10x,a)')'energy', &
                  'DOS-MT Up','DOS-MT Down', 'DOS-WS Up','DOS-WS Down'
         endif
         write(funit,'(64(''-''))')
      enddo
!
      ksize = nume*npola*2
!
      do j = 1, LocalNumAtoms
         iglb = getGlobalIndex(j)
         it = getAtomType(iglb)
!        -------------------------------------------------------------
         call daxpy(ksize,ONE,dos(1:nume,1:npola,1:2,j),1,                &
                    fspace(1:nume,1:npola,1:2,it),1)
!        -------------------------------------------------------------
      enddo
!
      do i = 1, NumPEsInGroup-1
!        -------------------------------------------------------------
         call recvMessage(msgbuf,MessageSize,MessageID,-1)
!        -------------------------------------------------------------
         pe = getSource()
         nla = getLocalNumAtoms(pe)
         do j = 1, nla
            iglb = getGlobalIndex(j,pe)
            it = getAtomType(iglb)
!           ----------------------------------------------------------
            call daxpy(ksize,ONE,msgbuf(ksize*(j-1)+1:ksize*j),1,     &
                       fspace(1:nume,1:npola,1:2,it),1)
!           ----------------------------------------------------------
         enddo
      enddo
!
      do it = 1, NumAtomTypes
         alp=ONE/real(getNumAtomsOfType(it),kind=RealKind)
         fspace(1:nume,1:npola,1:2,it) = fspace(1:nume,1:npola,1:2,it)*alp
         if (npola == 1) then
            do ie = 1, nume
               write(funit,'(2x,2(1x,d15.8),4x,2(2x,d15.8))')           &
                          egrd(ie),fspace(ie,1,1,it),fspace(ie,1,2,it)
            enddo
         else
            do ie = 1, nume
               write(funit,'(2x,2(1x,d15.8),4x,2(2x,d15.8,2x,d15.8))') &
                          egrd(ie),fspace(ie,1,1,it),fspace(ie,2,1,it),  &
                                   fspace(ie,1,2,it),fspace(ie,2,2,it)
            enddo
         endif
      enddo
      close(unit=funit)
      deallocate( fspace )
   endif
   deallocate( msgbuf )
!
   call resetCommunicator()
!
   end subroutine averageDOS
