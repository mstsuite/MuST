module LSMSConductivityModule
   use KindParamModule, only : IntKind, QuadRealKind, QuadCmplxKind, RealKind, CmplxKind, LongIntKind
   use MathParamModule, only : PI, ZERO, CZERO, CONE, TEN2m6, TEN2m7, TEN2m8, HALF, SQRTm1
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
   use PublicTypeDefinitionsModule, only : NeighborStruct
   use NeighborModule, only : getNeighbor, sortNeighbors

public :: initLSMSConductivity, &
          endLSMSConductivity,  &
          calSigmaTildeLSMS

private
   integer (kind=IntKind) :: LocalNumAtoms, GlobalNumAtoms
   integer (kind=IntKind) :: dsize, n_spin_pola, total_pairs, total_neighbors
   integer (kind=IntKind) :: numspecies, maxnumneighbors
   integer (kind=IntKind), allocatable :: num_neighbors(:)
   real (kind=RealKind) :: omega
   complex (kind=CmplxKind) :: eval
   complex (kind=CmplxKind), pointer :: tau(:,:), tau00(:,:,:)

   type TauContainerType
      complex (kind=CmplxKind), allocatable :: taup(:,:)
      complex (kind=CmplxKind), allocatable :: taun(:,:)
   end type TauContainerType
  
   type TauNeighborType
      type (TauContainerType), allocatable :: neigh1j(:)
      type (TauContainerType), allocatable :: neighj1(:)
   end type TauNeighborType
 
   type (TauNeighborType), allocatable, target :: TauIJ(:)

   type (NeighborStruct), pointer :: Neighbor

   integer (kind=IntKind) :: num_dirs, num_es
   integer (kind=IntKind) :: NumPEsInGroup, MyPEinGroup, GroupID

   complex (kind=CmplxKind), allocatable, target :: remote_jmtx(:,:)
   complex (kind=CmplxKind), allocatable, target :: local_jmtx(:,:)
#ifdef OpenMPI
   complex (kind=CmplxKind), allocatable, target :: jmtx_buf(:,:,:)
#else
   complex (kind=CmplxKind), allocatable :: jmtx_buf(:,:)
#endif
!
   integer (kind=IntKind) :: print_level
   integer (kind=IntKind) :: NumNonLocalNeighbors
   integer (kind=IntKind) :: NumReceives
   integer (kind=IntKind) :: NumSends
   integer (kind=IntKind), pointer :: NonLocalNeighbors(:)
   integer (kind=IntKind), pointer :: TargetProc(:)
   integer (kind=IntKind), pointer :: SourceProc(:)
   integer (kind=IntKind), allocatable :: remote_jmsize(:)
   integer (kind=IntKind), allocatable :: mapping_jmtx(:,:)
   integer (kind=IntKind), allocatable :: recv_msgid(:), send_msgid(:)

contains
!
   include '../lib/arrayTools.F90'
!  ============================================================

!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initLSMSConductivity(ndir, npola, kmax, efermi, iprint)
!  ============================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
   use KuboDataModule, only : isFermiEnergyRealPart, &
      getFermiEnergyRealPart, useStepFunctionForSigma, getFermiEnergyImagPart
   use SystemModule, only : getNumAtoms
   use Atom2ProcModule, only : getLocalNumAtoms, getGlobalIndex
   use SSSolverModule, only : solveSingleScattering, getScatteringMatrix
   use ClusterMatrixModule, only : getTau, calClusterMatrix, &
              calClusterMatrixNonPeriodic, getClusterTau, getNeighborTau, checkIfNeighbor
   use CurrentMatrixModule, only : calCurrentMatrix, getJMatrix, getLindex, getJMatrixSize
   use WriteMatrixModule, only : writeMatrix
   use SystemModule, only : getAtomPosition
   use SystemVolumeModule, only : getSystemVolume
   use PolyhedraModule, only : getVolume
   use SystemModule, only : getAtomTypeName, getAtomType, getNumAtomTypes

   integer (kind=IntKind), intent(in) :: ndir, npola, kmax, iprint
   real (kind=RealKind), intent(in) :: efermi
   integer (kind=IntKind) :: shin, pot_type, id, i, j, L1, L2
   integer (kind=IntKind) :: jmsize, jmrank, ne
   real (kind=RealKind) :: delta
   real (kind=RealKind) :: posi(3), posj(3), dist
   complex (kind=CmplxKind) :: eneg, kappa
   complex (kind=CmplxKind), allocatable :: temp(:,:)
   complex (kind=CmplxKind), pointer :: Jx(:,:), Jy(:,:), Jz(:,:)

   numspecies = getNumAtomTypes()
   delta = getFermiEnergyImagPart()
   pot_type = useStepFunctionForSigma()
   LocalNumAtoms = getLocalNumAtoms()
   GlobalNumAtoms = getNumAtoms()
   n_spin_pola = npola
   omega = getSystemVolume()
   dsize = kmax
   print_level = iprint
   total_pairs = 0

   GroupID = getGroupID('Unit Cell')
   NumPEsInGroup = getNumPEsInGroup(GroupID)
   MyPEinGroup = getMyPEinGroup(GroupID)
  
   if (isFermiEnergyRealPart()) then
      eval = getFermiEnergyRealPart() + SQRTm1*delta
      eneg = getFermiEnergyRealPart() - SQRTm1*delta  
   else
      eval = efermi + SQRTm1*delta
      eneg = efermi - SQRTm1*delta
   endif

   kappa = -sqrt(eval)

   allocate(num_neighbors(LocalNumAtoms))
   do i = 1, LocalNumAtoms
      Neighbor => getNeighbor(i)
      num_neighbors(i) = Neighbor%NumAtoms + 1
   enddo
   maxnumneighbors = maxval(num_neighbors)

   allocate(TauIJ(LocalNumAtoms))
   do i = 1, LocalNumAtoms
      Neighbor => getNeighbor(i)
      allocate(TauIJ(i)%neigh1j(Neighbor%NumAtoms+1))
      allocate(TauIJ(i)%neighj1(Neighbor%NumAtoms+1))
      do j = 1, Neighbor%NumAtoms+1
         allocate(TauIJ(i)%neigh1j(j)%taup(kmax, kmax))
         allocate(TauIJ(i)%neigh1j(j)%taun(kmax, kmax))
         allocate(TauIJ(i)%neighj1(j)%taup(kmax, kmax))
         allocate(TauIJ(i)%neighj1(j)%taun(kmax, kmax))
    !  TauIJ(i,j)%taup = CZERO
    !  TauIJ(i,j)%taun = CZERO
      enddo
   enddo

!  Solving for eF - i*delta

!  do id = 1, LocalNumAtoms
!    print *, "Atom number ", id
!    --------------------------------------------------------------------------
!    call solveSingleScattering(spin=is, site=id, e=eneg, vshift=CZERO, kp=kappa)
!    --------------------------------------------------------------------------
!  enddo

!  call calClusterMatrix(energy=eneg, &
!     getSingleScatteringMatrix=getScatteringMatrix, tau_needed=.true., kp=kappa)

!  do i = 1, LocalNumAtoms
!    Neighbor => getNeighbor(i)
!    do j = 1, Neighbor%NumAtoms+1
!     if (checkIfNeighbor(i,j)) then
!        TauIJ(i)%neigh1j(j)%taun = getNeighborTau(i, j, 0)
!        TauIJ(i)%neighj1(j)%taun = getNeighborTau(i, j, 1)
!      endif
!    enddo
!  enddo

!  Solving for eF + i*delta

   do is = 1, n_spin_pola
      do id = 1, LocalNumAtoms
!        ----------------------------------------------------------------------
         call solveSingleScattering(spin=is, site=id, e=eval, vshift=CZERO)
!        ----------------------------------------------------------------------
         call calCurrentMatrix(n=id, is=is, eval=eval, pot_type=1, mode=2)
!        ----------------------------------------------------------------------
      enddo

      call calClusterMatrix(energy=eval, &
                            getSingleScatteringMatrix=getScatteringMatrix, tau_needed=.true.)

!     call calClusterMatrixNonPeriodic(energy=eval, &
!                                      getSingleScatteringMatrix=getScatteringMatrix, tau_needed=.true.)

!     do i = 1, LocalNumAtoms
!        do j = 1, LocalNumAtoms
!           TauIJ(i,j)%taup = getClusterTau(i, j)
!        enddo
!     enddo 

      do i = 1, LocalNumAtoms
         Neighbor => getNeighbor(i)
         do j = 1, Neighbor%NumAtoms+1
            TauIJ(i)%neigh1j(j)%taup = getNeighborTau(i,j,0)
            TauIJ(i)%neighj1(j)%taup = getNeighborTau(i,j,1)
         enddo
      enddo

!     call writeMatrix('Tau12', TauIJ(1,2)%taup, dsize, dsize)
      do i = 1, LocalNumAtoms
         Neighbor => getNeighbor(i)
         do j = 1, Neighbor%NumAtoms+1
            do L1 = 1, dsize
               do L2 = 1, dsize
                  TauIJ(i)%neigh1j(j)%taun(L1,L2) = (-1.0)**(getLindex(L1) - getLindex(L2))* &
                                                    conjg(TauIJ(i)%neighj1(j)%taup(L2,L1))
                  TauIJ(i)%neighj1(j)%taun(L1,L2) = (-1.0)**(getLindex(L1) - getLindex(L2))* &
                                                    conjg(TauIJ(i)%neigh1j(j)%taup(L2,L1))
               enddo
            enddo
         enddo
      enddo

      if (print_level >= 1) then
         do i = 1, LocalNumAtoms
            Neighbor => getNeighbor(i)
            do j = 1, Neighbor%NumAtoms+1
               if (j == 1) then
                  write(6,'(a,2i4,2x,2d15.8)')'tau1j = ',getGlobalIndex(i),getGlobalIndex(i),TauIJ(i)%neigh1j(j)%taup(1,1)
                  write(6,'(a,2i4,2x,2d15.8)')'tauj1 = ',getGlobalIndex(i),getGlobalIndex(i),TauIJ(i)%neighj1(j)%taup(1,1)
               else
                  write(6,'(a,2i4,2x,2d15.8)')'tau1j = ',getGlobalIndex(i),Neighbor%GlobalIndex(j-1),TauIJ(i)%neigh1j(j)%taup(1,1)
                  write(6,'(a,2i4,2x,2d15.8)')'tauj1 = ',getGlobalIndex(i),Neighbor%GlobalIndex(j-1),TauIJ(i)%neighj1(j)%taup(1,1)
               endif
            enddo
            write(6,'(a,i4,a,i4)')'Convergence Data For Atom ',getGlobalIndex(i),',  Species ',getAtomType(i)
            write(6,'(a)')'Index,  Species,  (1,1) element'
            do j = 2, Neighbor%NumAtoms+1
               shin = Neighbor%ShellIndex(j-1)
               write(6,'(2i4,f12.8,2x,d15.8)') Neighbor%GlobalIndex(j-1), getAtomType(Neighbor%GlobalIndex(j-1)), &
                                               Neighbor%ShellRad(shin), real(TauIJ(i)%neigh1j(j)%taup(1,1),kind=RealKind)
            enddo

            write(6,'(a)')'Index,  Species,  (2,2) element'
            do j = 2, Neighbor%NumAtoms+1
               shin = Neighbor%ShellIndex(j-1)
               write(6,'(2i4,f12.8,2x,d15.8)') Neighbor%GlobalIndex(j-1), getAtomType(Neighbor%GlobalIndex(j-1)), &
                                               Neighbor%ShellRad(shin), real(TauIJ(i)%neigh1j(j)%taup(2,2),kind=RealKind)
            enddo
     
            write(6,'(a)')'Index,  Species,  (3,3) element'
            do j = 2, Neighbor%NumAtoms+1
               shin = Neighbor%ShellIndex(j-1)
               write(6,'(2i4,f12.8,2x,d15.8)') Neighbor%GlobalIndex(j-1), getAtomType(Neighbor%GlobalIndex(j-1)), &
                                               Neighbor%ShellRad(shin), real(TauIJ(i)%neigh1j(j)%taup(3,3),kind=RealKind)
            enddo
         enddo
      endif
   enddo
!
   jmsize = getJMatrixSize(jmrank)
   num_dirs = ndir
   num_es = 4
!  -------------------------------------------------------------------
   call allocateExchJmat(n_spin_pola*num_dirs*num_es*jmsize) ! Allocate space for exchanging the current Matrices
                                                             ! between processes
!  -------------------------------------------------------------------
   call exchangeCurrentMatrix()
!  ----------------------------------------------------------------------------
  
   end subroutine initLSMSConductivity
!  ============================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endLSMSConductivity()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: i, j
!
   deallocate(num_neighbors)
   do i = 1, LocalNumAtoms
      Neighbor => getNeighbor(i)
      do j = 1, Neighbor%NumAtoms+1
         deallocate(TauIJ(i)%neigh1j(j)%taup)
         deallocate(TauIJ(i)%neigh1j(j)%taun)
         deallocate(TauIJ(i)%neighj1(j)%taup)
         deallocate(TauIJ(i)%neighj1(j)%taun)
      enddo
      deallocate(TauIJ(i)%neigh1j)
      deallocate(TauIJ(i)%neighj1)
   enddo
   deallocate(TauIJ)
   nullify(Neighbor)

   call cleanExchJmat()
!
   end subroutine endLSMSConductivity
!  ===================================================================
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTildeLSMS(is, dir1, dir2, caltype) result(sigmaval)
!  ============================================================

   use ClusterMatrixModule, only : getClusterTau
   use CurrentMatrixModule, only : getJMatrix, getJMatrixSize
   use WriteMatrixModule, only : writeMatrix
   use SystemModule, only : getAtomType
   use GroupCommModule, only : GlobalSumInGroup
!
   implicit none

   integer (kind=IntKind), intent(in) :: is, dir1, dir2, caltype
   integer (kind=IntKind) :: m, n, L, A, B, np, jmsize, jmrank, jid
   complex (kind=CmplxKind) :: coeff, sigmaval, trace, tmpsum
   complex (kind=CmplxKind) :: sigmamat(numspecies, numspecies)
   complex (kind=CmplxKind) :: sigmaatom(LocalNumAtoms)
   complex (kind=CmplxKind), pointer :: Jm(:,:), Jn(:,:), p_jmtx(:)
   complex (kind=CmplxKind), allocatable :: temp(:,:), temp2(:,:), temp3(:,:)
!
   allocate(temp(dsize, dsize), temp2(dsize, dsize), temp3(dsize, dsize))
   temp = CZERO; temp2 = CZERO; temp3 = CZERO

   jmsize = getJMatrixSize(jmrank)
   np = (is-1)*num_dirs*num_es*jmsize+(dir2-1)*num_es*jmsize+(caltype-1)*jmsize
!
   sigmamat = CZERO
   sigmaatom = CZERO
   sigmaval = CZERO
   coeff = -CONE/(PI*omega)
   do m = 1, LocalNumAtoms
      A = getAtomType(m)
      if (caltype == 2) then
         Jm => getJMatrix(n=m,ic=1,is=is,dir=dir1,en_type=3,tilde=0)
      else if (caltype == 3) then
         Jm => getJMatrix(n=m,ic=1,is=is,dir=dir1,en_type=2,tilde=0)
      else
         Jm => getJMatrix(n=m,ic=1,is=is,dir=dir1,en_type=caltype,tilde=0)
      endif
      Neighbor => getNeighbor(m)
      do n = 1, Neighbor%NumAtoms+1
         trace = CZERO; temp = CZERO; temp2 = CZERO; temp3 = CZERO
         if (n == 1) then
            B = A
            Jn => getJMatrix(n=m,ic=1,is=is,dir=dir2,en_type=caltype,tilde=0)
         else if (Neighbor%ProcIndex(n-1) == MyPEinGroup) then
            B = getAtomType(Neighbor%LocalIndex(n-1))
            Jn => getJMatrix(n=Neighbor%LocalIndex(n-1),ic=1,is=is,dir=dir2,en_type=caltype,tilde=0)
         else
            jid = mapping_jmtx(n-1,m)
            p_jmtx => remote_jmtx(np+1:np+jmsize,jid)
            Jn => aliasArray2_c(p_jmtx,jmrank,jmrank)
         endif
  !    call writeMatrix('Jn', Jn, dsize, dsize)
         if (caltype == 1) then
            call zgemm('n', 'n', dsize, dsize, dsize, CONE, Jn, dsize, &
                       TauIJ(m)%neighj1(n)%taup, dsize, CZERO, temp, dsize)
  !         call writeMatrix('TauMN', TauIJ(m,n)%taup, dsize, dsize)
  !         call writeMatrix('TauNM', TauIJ(n,m)%taup, dsize, dsize)
            call zgemm('n', 'n', dsize, dsize, dsize, CONE, &
                       TauIJ(m)%neigh1j(n)%taup, dsize, temp, dsize, CZERO, temp2, dsize)
         else if (caltype == 2) then
            call zgemm('n', 'n', dsize, dsize, dsize, CONE, Jn, dsize, &
                       TauIJ(m)%neighj1(n)%taun, dsize, CZERO, temp, dsize)
            call zgemm('n', 'n', dsize, dsize, dsize, CONE, &
                       TauIJ(m)%neigh1j(n)%taup, dsize, temp, dsize, CZERO, temp2, dsize)
         else if (caltype == 3) then
            call zgemm('n', 'n', dsize, dsize, dsize, CONE, Jn, dsize, &
                       TauIJ(m)%neighj1(n)%taup, dsize, CZERO, temp, dsize)
            call zgemm('n', 'n', dsize, dsize, dsize, CONE, &
                       TauIJ(m)%neigh1j(n)%taun, dsize, temp, dsize, CZERO, temp2, dsize)
         else if (caltype == 4) then
            call zgemm('n', 'n', dsize, dsize, dsize, CONE, Jn, dsize, &
                       TauIJ(m)%neighj1(n)%taun, dsize, CZERO, temp, dsize)
            call zgemm('n', 'n', dsize, dsize, dsize, CONE, &
                       TauIJ(m)%neigh1j(n)%taun, dsize, temp, dsize, CZERO, temp2, dsize)
         endif
         call zgemm('n', 'n', dsize, dsize, dsize, CONE, Jm, dsize, &
                    temp2, dsize, CZERO, temp3, dsize)
         do L = 1, dsize
            trace = trace + coeff*temp3(L,L)
         enddo
         sigmaval = sigmaval + trace
         sigmamat(A,B) = sigmamat(A,B) + trace
         sigmaatom(m) = sigmaatom(m) + trace
     ! print *, "Sigmaval for (m,n) ", m, n, " in direction (mu,nu) ", dir1, dir2, " and caltype", caltype, "is "
     ! print *, real(trace), "with total", real(sigmaval) 
      enddo
   enddo

   if (print_level >= 1) then
      do m = 1, LocalNumAtoms
         write(6,'(a,i4,a,2i4,a,i4)')'Conductivity = ', caltype, ', along dirs = ', dir1, dir2, ', for Atom = ', m
         write(6,'(a,d15.8)')'Re[sigmaatom] = ',real(sigmaatom(m),kind=RealKind)
      enddo

      do m = 1, numspecies
         do n = 1, numspecies
            write(6,'(a,i4,a,2i4,a,i4)')'Conductivity = ', caltype, ', along dirs = ', dir1, dir2, ', for Type = ', m, n
            write(6,'(a,d15.8)')'Re[sigmamat] = ', real(sigmamat(m,n),kind=RealKind)
         enddo
      enddo
   endif
!
!  ===================================================================
!  Notes: Will perform global sum of sigma here ...
!  -------------------------------------------------------------------
   call GlobalSumInGroup(GroupID,sigmaval)
!  -------------------------------------------------------------------
!
   deallocate(temp, temp2, temp3)
   nullify(Jn, Jm, p_jmtx)

   end function calSigmaTildeLSMS
!  ============================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine exchangeCurrentMatrix()
!  ===================================================================
   use MPPModule, only : AnyPE, MyPE
   use MPPModule, only : nbsendMessage, nbrecvMessage
   use MPPModule, only : recvMessage, isThereMessage, getSource, setMaxWaits
   use MPPModule, only : waitMessage, setCommunicator, resetCommunicator
!
   use GroupCommModule, only : getGroupCommunicator
!
   use DataServiceCenterModule, only : getDataStorage, ComplexMark
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use SpinRotationModule, only : rotateLtoG
!
   use CurrentMatrixModule, only : getJMatrixSize, getJMatrix
!
   implicit none
!
   integer (kind=IntKind) :: comm
   integer (kind=IntKind) :: jmsize, jmrank, nfac
   integer (kind=IntKind) :: i, is, je, dir, nr, np, ns, nm, mj, NumFills
   integer (kind=IntKind) :: j, k, proc, gid
!
   complex (kind=CmplxKind), pointer :: jm(:,:)
#ifdef OpenMPI
   complex (kind=CmplxKind), pointer :: p_trecv(:,:)
#endif
!
   jmsize = getJMatrixSize(jmrank)
   nfac = n_spin_pola*num_es*num_dirs
!
!  -------------------------------------------------------------------
   comm = getGroupCommunicator(GroupID)
   call setCommunicator(comm,MyPEinGroup,NumPEsInGroup,sync=.true.)
!  -------------------------------------------------------------------
!
   do i = 1, LocalNumAtoms
!     ================================================================
!     Here, we assume the local atoms have the same matrix size (= jmrank)
!     ================================================================
      local_jmtx(1,i) = cmplx(jmrank, getGlobalIndex(i), kind=CmplxKind)
      nm = 0
      do is = 1, n_spin_pola
!        =============================================================
!        Obtain the current matrix
!        =============================================================
         do dir = 1, num_dirs
            do je = 1, num_es
               jm => getJMatrix(n=i,ic=1,is=is,dir=dir,en_type=je,tilde=0)
!              =======================================================
!              Store the j-matrix into the local communication buffer
!              -------------------------------------------------------
               call zcopy(jmsize,jm,1,local_jmtx(nm+2,i),1)
!              -------------------------------------------------------
               nm = nm + jmsize
            enddo
         enddo
      enddo
   enddo
!
   nullify(jm)
!
   nr = 0; ns = 0; NumFills = 0
   remote_jmsize = 0
!
#ifdef OpenMPI
   nm = max(NumReceives, NumSends)
   call setMaxWaits(nm)
   do k = 1, nm
      if (k <= NumReceives) then
         p_trecv => jmtx_buf(1:nfac*jmsize+1,1:LocalNumAtoms,k)
!        -------------------------------------------------------------
         recv_msgid(k)=nbrecvMessage(p_trecv,nfac*jmsize+1,LocalNumAtoms,2010,AnyPE)
!        ------------------------------------------------------------
      endif
   enddo
!
   do k = 1, nm
      if (k < NumSends+1) then
         proc = TargetProc(k)
!        -------------------------------------------------------------
         send_msgid(k) = nbsendMessage(local_jmtx,nfac*jmsize+1,LocalNumAtoms,2010,proc)
!        -------------------------------------------------------------
      endif
   enddo
!
   do k = 1, nm
      if (nr < NumReceives) then
         p_trecv => jmtx_buf(1:nfac*jmsize+1,1:LocalNumAtoms,nr+1)
!        -------------------------------------------------------------
         call waitMessage(recv_msgid(k))
!        -------------------------------------------------------------
         nr = nr + 1
!
         do i = 1, LocalNumAtoms
            mj = int(real(jmtx_buf(1,i,nr),kind=RealKind))
            gid = int(aimag(jmtx_buf(1,i,nr)))
            LOOP_j: do j = 1, NumNonLocalNeighbors
               if (gid == NonLocalNeighbors(j)) then
!                 ----------------------------------------------------
                  call zcopy(nfac*jmsize,jmtx_buf(2,i,nr),1,remote_jmtx(1,j),1)
!                 ----------------------------------------------------
!                 remote_jmtx(1:nfac*jmsize,j) = jmtx_buf(2:nfac*jmsize+1,i,nr)
                  remote_jmsize(j) = mj
                  NumFills = NumFills + 1
                  exit LOOP_j
               endif
            enddo LOOP_j
         enddo
      else if (NumFills /= NumNonLocalNeighbors) then
!        -------------------------------------------------------------
         call ErrorHandler('exchangeSSSMatrix','NumFills /= NumNonLocalNeighbors', &
                           NumFills,NumNonLocalNeighbors)
!        -------------------------------------------------------------
      endif
   enddo
!
   do k = 1, nm
      if (k < NumSends+1) then
!        -------------------------------------------------------------
         call waitMessage(send_msgid(k))
!        -------------------------------------------------------------
         ns = ns + 1
      endif
   enddo
   nullify(p_trecv)
!
#else
   do while (nr < NumReceives .or. ns < NumSends)
      if (ns < NumSends) then
         ns = ns + 1
         proc = TargetProc(ns)
         send_msgid(ns) = nbsendMessage(local_jmtx,nfac*jmsize+1,LocalNumAtoms,2010,proc)
      endif
!
      if (nr < NumReceives) then
         if (isThereMessage(2010,AnyPE)) then
!           ----------------------------------------------------------
            call recvMessage(jmtx_buf,nfac*jmsize+1,LocalNumAtoms,2010,getSource())
!           ----------------------------------------------------------
            nr = nr + 1
!
            do i = 1, LocalNumAtoms
               mj = int(real(jmtx_buf(1,i),kind=RealKind))
               gid = int(aimag(jmtx_buf(1,i)))
               LOOP_j: do j = 1, NumNonLocalNeighbors
                  if (gid == NonLocalNeighbors(j)) then
!                    -------------------------------------------------
                     call zcopy(nfac*jmsize,jmtx_buf(2,i),1,remote_jmtx(1,j),1)
!                    -------------------------------------------------
!                    remote_jmtx(1:nfac*jmsize,j) = jmtx_buf(2:nfac*jmsize+1,i)
                     remote_jmsize(j) = mj
                     NumFills = NumFills + 1
                     exit LOOP_j
                  endif
               enddo LOOP_j
            enddo
         endif
      endif
   enddo
!
   if (NumFills /= NumNonLocalNeighbors) then
!     ----------------------------------------------------------------
      call ErrorHandler('exchangeSSSMatrix','NumFills /= NumNonLocalNeighbors', &
                        NumFills,NumNonLocalNeighbors)
!     ----------------------------------------------------------------
   endif
!
   do k = 1, NumSends
      call waitMessage(send_msgid(k))
   enddo
#endif
!
!  -------------------------------------------------------------------
   call resetCommunicator(sync=.true.)
!  -------------------------------------------------------------------
!
   end subroutine exchangeCurrentMatrix
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine allocateExchJmat(bsize)
!  ===================================================================
   use SendRecvTmatModule, only : getMaxNeighbors
   use SendRecvTmatModule, only : getNumNonLocalNeighbors, getNonLocalNeighbors
   use SendRecvTmatModule, only : getNumProcWiseSends, getNumProcWiseReceives
   use SendRecvTmatModule, only : getProcWiseTargetProcs, getProcWiseSourceProcs
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: bsize
   integer (kind=IntKind) :: m, n, i, j, iri, pe, gid
!
!  ===================================================================
!  We assume that bsize is the same on all processors. This needs to be
!  fixed in the future.
!  ===================================================================
   m = getMaxNeighbors()
   NumNonLocalNeighbors = getNumNonLocalNeighbors()
!
   NonLocalNeighbors => getNonLocalNeighbors()
   NumSends = getNumProcWiseSends()
   NumReceives = getNumProcWiseReceives()
   TargetProc => getProcWiseTargetProcs()
   SourceProc => getProcWiseSourceProcs()
!
   allocate( local_jmtx(bsize+1,LocalNumAtoms) )
   allocate(recv_msgid(NumReceives),send_msgid(NumSends))
#ifdef OpenMPI
   allocate( jmtx_buf(bsize+1,LocalNumAtoms,NumReceives) )
#else
   allocate( jmtx_buf(bsize+1,LocalNumAtoms) )
#endif
   allocate( mapping_jmtx(m,LocalNumAtoms) )
!
   mapping_jmtx = -1
   local_jmtx = CZERO
   jmtx_buf = CZERO
!
   do i = 1, LocalNumAtoms
      Neighbor => getNeighbor(i)
      do iri=1, Neighbor%NumAtoms
         pe = Neighbor%ProcIndex(iri)
         if (pe /= MyPEinGroup) then
            gid = Neighbor%GlobalIndex(iri)
            LOOP_j: do j = 1, NumNonLocalNeighbors
               if (gid == NonLocalNeighbors(j)) then
                  mapping_jmtx(iri,i) = j
                  exit LOOP_j
               endif
            enddo LOOP_j
         endif
      enddo
   enddo
   allocate( remote_jmsize(NumNonLocalNeighbors), remote_jmtx(bsize,NumNonLocalNeighbors) )
!
   end subroutine allocateExchJmat
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine cleanExchJmat()
!  ===================================================================
   implicit none
!
   NumNonLocalNeighbors = 0
   NumSends = 0
   NumReceives = 0
!
   deallocate(remote_jmsize, remote_jmtx, local_jmtx, jmtx_buf)
   deallocate(mapping_jmtx)
   deallocate(recv_msgid,send_msgid)
!
   nullify(NonLocalNeighbors, SourceProc, TargetProc)
!
   end subroutine cleanExchJmat
!  ===================================================================
end module LSMSConductivityModule
