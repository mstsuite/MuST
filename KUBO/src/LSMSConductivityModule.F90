module LSMSConductivityModule
   use KindParamModule, only : IntKind, QuadRealKind, QuadCmplxKind, RealKind, CmplxKind, LongIntKind
   use MathParamModule, only : PI, ZERO, CZERO, CONE, TEN2m6, TEN2m7, TEN2m8, HALF, SQRTm1
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler, StopHandler
   use PublicTypeDefinitionsModule, only : NeighborStruct
   use NeighborModule, only : getNeighbor, sortNeighbors

public :: initLSMSConductivity, &
          calSigmaTildeLSMS, &
          endLSMSConductivity

private
   integer (kind=IntKind) :: LocalNumAtoms, GlobalNumAtoms
   integer (kind=IntKind) :: dsize, spin_pola
   real (kind=RealKind) :: omega
   complex (kind=CmplxKind) :: eval
   complex (kind=CmplxKind), pointer :: tau(:,:), tau00(:,:,:)

   type TauContainerType
      complex (kind=CmplxKind), allocatable :: taup(:,:)
      complex (kind=CmplxKind), allocatable :: taun(:,:)
   end type TauContainerType
   
   type (TauContainerType), allocatable, target :: TauIJ(:,:)

contains
!  ============================================================

!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initLSMSConductivity(is, kmax, efermi)
!  ============================================================

   use KuboDataModule, only : isFermiEnergyRealPart, &
      getFermiEnergyRealPart, useStepFunctionForSigma, getFermiEnergyImagPart
   use SystemModule, only : getNumAtoms
   use Atom2ProcModule, only : getLocalNumAtoms
   use SSSolverModule, only : solveSingleScattering, getScatteringMatrix
   use ClusterMatrixModule, only : getTau, calClusterMatrix, getNeighborTau, checkIfNeighbor
   use CurrentMatrixModule, only : calCurrentMatrix, getJMatrix, getLindex
   use WriteMatrixModule, only : writeMatrix
   use SystemVolumeModule, only : getSystemVolume

   integer (kind=IntKind), intent(in) :: is, kmax
   real (kind=RealKind), intent(in) :: efermi
   integer (kind=IntKind) :: pot_type, id, i, j, L1, L2
   real (kind=RealKind) :: delta
   complex (kind=CmplxKind) :: eneg
   complex (kind=CmplxKind), pointer :: temp(:,:), temp2(:,:)
   
   delta = getFermiEnergyImagPart()
   pot_type = useStepFunctionForSigma()
   LocalNumAtoms = getLocalNumAtoms()
   GlobalNumAtoms = getNumAtoms()
   spin_pola = is
   omega = getSystemVolume()
   dsize = kmax
  
   if (isFermiEnergyRealPart()) then
      eval = getFermiEnergyRealPart() + SQRTm1*delta
   else
      eval = efermi + SQRTm1*delta
   endif
  
   eneg = eval - 2*SQRTm1*delta
  
   allocate(TauIJ(LocalNumAtoms, LocalNumAtoms))
   do i = 1, LocalNumAtoms
     do j = 1, LocalNumAtoms
       allocate(TauIJ(i,j)%taup(kmax, kmax))
       allocate(TauIJ(i,j)%taun(kmax, kmax))
       TauIJ(i,j)%taup = CZERO
       TauIJ(i,j)%taun = CZERO
     enddo
   enddo

!  Solving for eF - i*delta

!  do id = 1, LocalNumAtoms
!    print *, "Atom number ", id
!    --------------------------------------------------------------------------
!    call solveSingleScattering(spin=spin_pola, site=id, e=eneg, vshift=CZERO)
!    --------------------------------------------------------------------------
!  enddo

!  call calClusterMatrix(energy=eneg, &
!     getSingleScatteringMatrix=getScatteringMatrix, tau_needed=.true.)

!  do i = 1, LocalNumAtoms
 !   tau00 => getTau(i)
 !   call writeMatrix('tau00', tau00, kmax, kmax, 1)
!    do j = 1, LocalNumAtoms
!      if (checkIfNeighbor(i,j)) then
 !       print *, i, "and ", j, "are neighbors"
!        tau => getNeighborTau(i, j)
!        call zcopy(kmax*kmax, tau,1,TauIJ(j,i)%taun,1)
!      else
 !       print *, i, "and ", j, "are not neighbors!"
!      endif
!    enddo
!  enddo

!  Solving for eF + i*delta

   do id = 1, LocalNumAtoms
!    --------------------------------------------------------------------------
     call solveSingleScattering(spin=spin_pola, site=id, e=eval, vshift=CZERO)
!    --------------------------------------------------------------------------
     call calCurrentMatrix(n=id, is=spin_pola, eval=eval, pot_type=1, mode=2)
!    --------------------------------------------------------------------------
   enddo

   call calClusterMatrix(energy=eval, &
       getSingleScatteringMatrix=getScatteringMatrix, tau_needed=.true.)

   do i = 1, LocalNumAtoms
     do j = 1, LocalNumAtoms
       if (checkIfNeighbor(i,j)) then
         TauIJ(j,i)%taup = getNeighborTau(i,j)
   !     print *, "Tau matrix for (m,n) = ", i,j
   !     call writeMatrix('TauP', TauIJ(j,i)%taup, kmax, kmax)
        ! call zcopy(kmax*kmax, tau, 1, TauIJ(j,i)%taup,1)
         do L1 = 1, kmax
           do L2 = 1, kmax
             TauIJ(i,j)%taun(L1, L2) = (-1.0)**(getLIndex(L1) - getLIndex(L2))\
                      *conjg(TauIJ(j,i)%taup(L2, L1))
           enddo
         enddo
       endif
     enddo
   enddo
  ! do i = 1, LocalNumAtoms
  !  do j = 1, LocalNumAtoms
     ! print *, "Tau Matrix for (m,n) = ", i,j
     !  call writeMatrix('TauP',TauIJ(i,j)%taup,kmax,kmax)
     !  call writeMatrix('TauN',TauIJ(i,j)%taun,kmax,kmax)
     ! print *, "System Volume Per Atom", omega 
  !  enddo
  !enddo
  
   end subroutine initLSMSConductivity
!  ============================================================

!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTildeLSMS(dir1, dir2, caltype) result(sigmaval)
!  ============================================================

   use CurrentMatrixModule, only : getJMatrix
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: dir1, dir2, caltype
   integer (kind=IntKind) :: m, n, L, Jmcaltype
   complex (kind=CmplxKind) :: coeff, sigmaval, trace
   complex (kind=CmplxKind), pointer :: Jm(:,:), Jn(:,:)
   complex (kind=CmplxKind), allocatable :: temp(:,:), temp2(:,:), temp3(:,:)

   allocate(temp(dsize, dsize), temp2(dsize, dsize), temp3(dsize, dsize))
   temp = CZERO; temp2 = CZERO; temp3 = CZERO

   sigmaval = CZERO
   coeff = -CONE/(PI*omega)
   do m = 1, LocalNumAtoms
     if (caltype == 2 .or. caltype == 3) then
       Jmcaltype = 5 - caltype
       Jm => getJMatrix(n=m,ic=1,is=spin_pola,dir=dir1,en_type=Jmcaltype,tilde=0)
     else
       Jm => getJMatrix(n=m,ic=1,is=spin_pola,dir=dir1,en_type=caltype,tilde=0)
     endif
     do n = 1, LocalNumAtoms
       trace = CZERO
       Jn => getJMatrix(n=n,ic=1,is=spin_pola,dir=dir2,en_type=caltype,tilde=0)
       if (caltype == 1) then
         call zgemm('n', 'n', dsize, dsize, dsize, CONE, Jn, dsize, &
               TauIJ(n,m)%taup, dsize, CZERO, temp, dsize)
         call zgemm('n', 'n', dsize, dsize, dsize, CONE, &
           TauIJ(m,n)%taup, dsize, temp, dsize, CZERO, temp2, dsize)
       else if (caltype == 2) then
         call zgemm('n', 'n', dsize, dsize, dsize, CONE, Jn, dsize, &
               TauIJ(n,m)%taun, dsize, CZERO, temp, dsize)
         call zgemm('n', 'n', dsize, dsize, dsize, CONE, &
            TauIJ(m,n)%taup, dsize, temp, dsize, CZERO, temp2, dsize)
       else if (caltype == 3) then
         call zgemm('n', 'n', dsize, dsize, dsize, CONE, Jn, dsize, &
               TauIJ(n,m)%taup, dsize, CZERO, temp, dsize)
         call zgemm('n', 'n', dsize, dsize, dsize, CONE, &
            TauIJ(m,n)%taun, dsize, temp, dsize, CZERO, temp2, dsize)
       else if (caltype == 4) then
         call zgemm('n', 'n', dsize, dsize, dsize, CONE, Jn, dsize, &
               TauIJ(n,m)%taun, dsize, CZERO, temp, dsize)
         call zgemm('n', 'n', dsize, dsize, dsize, CONE, &
            TauIJ(m,n)%taun, dsize, temp, dsize, CZERO, temp2, dsize)
       endif
       call zgemm('n', 'n', dsize, dsize, dsize, CONE, Jm, dsize, &
             temp2, dsize, CZERO, temp3, dsize)
       do L = 1, dsize
         trace = trace + coeff*temp3(L,L)
       enddo
       sigmaval = sigmaval + trace
!      print *, "Sigmaval for (m,n) ", m, n, " in direction (mu,nu) ", dir1, dir2, " and caltype", caltype, "is "
!      print *, real(trace), "with total", real(sigmaval) 
     enddo
   enddo

   end function calSigmaTildeLSMS
!  ============================================================
end module LSMSConductivityModule
