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
   integer (kind=IntKind) :: dsize, spin_pola, total_pairs, total_neighbors
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
   use ClusterMatrixModule, only : getTau, calClusterMatrix, &
              calClusterMatrixNonPeriodic, getClusterTau, getNeighborTau, checkIfNeighbor
   use CurrentMatrixModule, only : calCurrentMatrix, getJMatrix, getLindex
   use WriteMatrixModule, only : writeMatrix
   use SystemModule, only : getAtomPosition
   use SystemVolumeModule, only : getSystemVolume
   use PolyhedraModule, only : getVolume

   integer (kind=IntKind), intent(in) :: is, kmax
   real (kind=RealKind), intent(in) :: efermi
   integer (kind=IntKind) :: pot_type, id, i, j, L1, L2
   real (kind=RealKind) :: delta
   real (kind=RealKind) :: posi(3), posj(3), dist
   complex (kind=CmplxKind) :: eneg, kappa
   complex (kind=CmplxKind), pointer :: temp(:,:,:), temp2(:,:)
   
   delta = getFermiEnergyImagPart()
   pot_type = useStepFunctionForSigma()
   LocalNumAtoms = getLocalNumAtoms()
   GlobalNumAtoms = getNumAtoms()
   spin_pola = is
   omega = getSystemVolume()
   dsize = kmax
   total_pairs = 0
  
   if (isFermiEnergyRealPart()) then
      eval = getFermiEnergyRealPart() + SQRTm1*delta
      eneg = getFermiEnergyRealPart() - SQRTm1*delta  
   else
      eval = efermi + SQRTm1*delta
      eneg = efermi - SQRTm1*delta
   endif

   kappa = -sqrt(eval)

 ! print *, "SQRT positive is ", sqrt(eval)
 ! print *, "SQRT negative is ", sqrt(eneg)
 ! print *, "What it should be", -sqrt(eval)  
 
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
!    call solveSingleScattering(spin=spin_pola, site=id, e=eneg, vshift=CZERO, kp=kappa)
!    --------------------------------------------------------------------------
!  enddo

!  call calClusterMatrix(energy=eneg, &
!     getSingleScatteringMatrix=getScatteringMatrix, tau_needed=.true., kp=kappa)

!  do i = 1, LocalNumAtoms
!    do j = 1, LocalNumAtoms
!      if (checkIfNeighbor(i,j)) then
!        TauIJ(j,i)%taun = getNeighborTau(i, j)
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

!  call calClusterMatrix(energy=eval, &
!      getSingleScatteringMatrix=getScatteringMatrix, tau_needed=.true.)

   call calClusterMatrixNonPeriodic(energy=eval, &
       getSingleScatteringMatrix=getScatteringMatrix, tau_needed=.true.)

   do i = 1, LocalNumAtoms
     do j = 1, LocalNumAtoms
       TauIJ(i,j)%taup = getClusterTau(i, j)
     enddo
   enddo 

!  do i = 1, LocalNumAtoms
!    do j = 1, LocalNumAtoms
!      if (checkIfNeighbor(i,j)) then
!        TauIJ(i,j)%taup = getNeighborTau(i,j)
!        total_pairs = total_pairs + 1
!      endif
!    enddo
!  enddo

   total_neighbors = total_pairs - LocalNumAtoms
   print *, "Total pairs ", total_pairs
   print *, "Total neighbors ", total_neighbors
 ! call writeMatrix('Tau12', TauIJ(1,2)%taup, dsize, dsize)
   do i = 1, LocalNumAtoms
    do j = 1, LocalNumAtoms
       do L1 = 1, dsize
         do L2 = 1, dsize
          TauIJ(i,j)%taun(L1,L2) = (-1.0)**(getLindex(L1) - getLindex(L2))*conjg(TauIJ(j,i)%taup(L2,L1))
         enddo
       enddo
     enddo
   enddo

!  call writeMatrix('Tau11n', TauIJ(1,1)%taun, dsize, dsize)
!  call writeMatrix('Tau12n', TauIJ(1,2)%taun, dsize, dsize)
!  call writeMatrix('Tau103LSMS', TauIJ(10, 3)%taup, dsize, dsize)
!  do i = 1, LocalNumAtoms
!    do j = 1, LocalNumAtoms
!      print *, "Tau Matrix for (m,n) = ", i,j
!      call writeMatrix('TauP',TauIJ(i,j)%taup,kmax,kmax)
!      call writeMatrix('TauN',TauIJ(i,j)%taun,kmax,kmax)
!      print *, "System Volume Per Atom", omega 
!    enddo
!  enddo
  
   end subroutine initLSMSConductivity
!  ============================================================

!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function calSigmaTildeLSMS(dir1, dir2, caltype) result(sigmaval)
!  ============================================================

   use CurrentMatrixModule, only : getJMatrix
   use WriteMatrixModule, only : writeMatrix

   integer (kind=IntKind), intent(in) :: dir1, dir2, caltype
   integer (kind=IntKind) :: m, n, L
   complex (kind=CmplxKind) :: coeff, sigmaval, trace, tmpsum
   complex (kind=CmplxKind) :: sigmamat(LocalNumAtoms,LocalNumAtoms)
   complex (kind=CmplxKind), pointer :: Jm(:,:), Jn(:,:)
   complex (kind=CmplxKind), allocatable :: temp(:,:), temp2(:,:), temp3(:,:)

   allocate(temp(dsize, dsize), temp2(dsize, dsize), temp3(dsize, dsize))
   temp = CZERO; temp2 = CZERO; temp3 = CZERO

   sigmamat = CZERO
   sigmaval = CZERO
   coeff = -CONE/(PI*omega)
   do m = 1, LocalNumAtoms
     if (caltype == 2) then
       Jm => getJMatrix(n=m,ic=1,is=spin_pola,dir=dir1,en_type=3,tilde=0)
     else if (caltype == 3) then
       Jm => getJMatrix(n=m,ic=1,is=spin_pola,dir=dir1,en_type=2,tilde=0)
     else
       Jm => getJMatrix(n=m,ic=1,is=spin_pola,dir=dir1,en_type=caltype,tilde=0)
     endif
     do n = 1, LocalNumAtoms
       trace = CZERO; temp = CZERO; temp2 = CZERO; temp3 = CZERO
       Jn => getJMatrix(n=n,ic=1,is=spin_pola,dir=dir2,en_type=caltype,tilde=0)
  !    call writeMatrix('Jn', Jn, dsize, dsize)
       if (caltype == 1) then
         call zgemm('n', 'n', dsize, dsize, dsize, CONE, Jn, dsize, &
               TauIJ(n,m)%taup, dsize, CZERO, temp, dsize)
  !      call writeMatrix('TauMN', TauIJ(m,n)%taup, dsize, dsize)
  !      call writeMatrix('TauNM', TauIJ(n,m)%taup, dsize, dsize)
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
       sigmamat(m,n) = trace
     ! print *, "Sigmaval for (m,n) ", m, n, " in direction (mu,nu) ", dir1, dir2, " and caltype", caltype, "is "
     ! print *, real(trace), "with total", real(sigmaval) 
     enddo
   enddo

   tmpsum = CZERO
   print *, "Total Single-site sum is "
   do m = 1, LocalNumAtoms
     tmpsum = tmpsum + sigmamat(m,m)
   enddo
   print *, real(tmpsum)
   tmpsum = CZERO
   print *, "Total Cross-site sum is "
   do m = 1, LocalNumAtoms
     do n = 1, LocalNumAtoms
       if (m .ne. n) then
         tmpsum = tmpsum + sigmamat(m,n)
       endif
     enddo
   enddo
   print *, tmpsum

   end function calSigmaTildeLSMS
!  ============================================================
end module LSMSConductivityModule
