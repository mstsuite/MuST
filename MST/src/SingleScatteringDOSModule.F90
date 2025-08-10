module SingleScatteringDOSModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, ONE, TWO, TEN2m6, Y0, PI, PI2, CZERO, SQRTm1, CONE
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
!
   use TimerModule, only : getTime
!
   use IntegerFactorsModule, only : lofk, lofj, mofj, m1m, mofk, jofk
!
   use MPPModule, only : MyPE
!
public :: setSScatteringDOSParam,       &
          getSScatteringPhaseShift,     &
          getSScatteringDOS,            & ! Returns results per spin
          getSScatteringDOSofCmplxE,    & ! Returns results per spin
          getRelSScatteringDOS,         & ! Returns results per spin
          getRelSScatteringDOSofCmplxE    ! Returns results per spin
!
   interface setSScatteringDOSParam
      module procedure setSSP_generic, setSSP_ef, setSSP_specific
   end interface setSScatteringDOSParam
!
private
!
   integer (kind=IntKind) :: n_spin_pola = 1
   integer (kind=IntKind) :: n_spin_cant = 1
   integer (kind=IntKind) :: node_print_level = 0
   integer (kind=IntKind) :: iharris = 0
!
   integer (kind=IntKind) :: NumRs = 0
   integer (kind=IntKind) :: jmax_green = 0
   integer (kind=IntKind) :: local_index = 0
!
   integer (kind=IntKind) :: NumPEsInEGroup, MyPEinEGroup, eGID
!
   logical :: rad_derivative
!
   real (kind=RealKind) :: chempot = ZERO
   real (kind=RealKind) :: Temperature = ZERO
!
   interface adjustEnergy
       function adjustEnergy_r(is,e) result(en)
          use KindParamModule, only : IntKind, RealKind, CmplxKind
          implicit none
          integer (kind=IntKind), intent(in) :: is
          real (kind=RealKind), intent(in) :: e
          real (kind=RealKind) :: en
       end function adjustEnergy_r
!
       function adjustEnergy_c(is,e) result(en)
          use KindParamModule, only : IntKind, RealKind, CmplxKind
          implicit none
          integer (kind=IntKind), intent(in) :: is
          complex (kind=CmplxKind), intent(in) :: e
          complex (kind=CmplxKind) :: en
       end function adjustEnergy_c
   end interface
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSSP_generic(ns,nc,dev,harris,iprint)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: ns, nc, harris, iprint
!
   logical, intent(in) :: dev
!
   n_spin_pola = ns
   n_spin_cant = nc
   rad_derivative = dev
   iharris = harris
   node_print_level = iprint
!
   eGID = getGroupID('Energy Mesh')
   NumPEsInEGroup = getNumPEsInGroup(eGID)
   MyPEinEGroup = getMyPEinGroup(eGID)
!
   end subroutine setSSP_generic
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSSP_specific(id,nr,jmax_in)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, nr, jmax_in
!
   local_index = id
   NumRs = nr
   jmax_green = jmax_in
!
   end subroutine setSSP_specific
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setSSP_ef(ef,temp)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: ef, temp
!
   chempot = ef
   Temperature = temp
!
   end subroutine setSSP_ef
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSScatteringPhaseShift(info,e) result(tps)
!  ===================================================================
   use SSSolverModule, only : solveSingleScattering, computePhaseShift, getPhaseShift
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: lmax_phi, kl, kmax_phi, is, id, atom
!
   real (kind=RealKind), intent(in) :: e
!
   complex (kind=CmplxKind) :: energy
!
   real (kind=RealKind) :: tps, t0
   real (kind=RealKind), pointer :: ps(:)
!
   is = info(1); id = info(2); atom = info(3); lmax_phi = info(5)
!
   if (id /= local_index) then
      call ErrorHandler('getSScatteringPhaseShift','Inconsistent local index',id,local_index)
   endif
!
   energy = adjustEnergy(is,e)
!
   if (atom < 1) then
      call ErrorHandler('getSScatteringPhaseShift','atom < 1',atom)
   endif
!
!  -------------------------------------------------------------------
   call solveSingleScattering(spin=is,site=id,atom=atom,e=energy,vshift=CZERO)
!  -------------------------------------------------------------------
   call computePhaseShift(spin=is,site=id,atom=atom)
   ps => getPhaseShift(spin=is,site=id,atom=atom)
!  -------------------------------------------------------------------
   tps = ZERO
   kmax_phi = (lmax_phi+1)**2
   do kl = 1, kmax_phi
      tps = tps + ps(kl)
   enddo
!
   end function getSScatteringPhaseShift
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSScatteringDOS(info,e,aux,rfac,redundant) result(ssdos)
!  ===================================================================
!
!  This function returns single site DOS for a given energy on the real 
!  energy axis.
!
!  The function also returns aux, a data set that can be integrated
!  on an energy grid
!
!  In the no-spin-polarized case, the returned dos is per spin.
!  ===================================================================
   use PhysParamModule, only : Boltzmann
!
   use GroupCommModule, only : GlobalSumInGroup
!
   use PotentialTypeModule, only : isASAPotential
!
   use RadialGridModule, only : getGrid, getRadialIntegration
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use SSSolverModule, only : solveSingleScattering, computeDOS, getDOS
   use SSSolverModule, only : getOutsideDOS, computePhaseShift, getPhaseShift
   use SSSolverModule, only : getDOSDerivative
!
   use AtomModule, only : getLocalAtomicNumber, getLocalAtomName
   use AtomModule, only : getLocalNumSpecies
!
   implicit none
!
   real (kind=RealKind), intent(in) :: e
   real (kind=RealKind), intent(in), optional :: rfac
!
   logical, intent(in), optional :: redundant
   logical :: red, InResonancePeak
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: jmax_dos, kmax_phi, iend, n, kl, jl, ir
   integer (kind=IntKind) :: is, id, print_dos, ia, i, atom, ib, ib0, lmax_phi
   integer (kind=IntKind) :: NumResonanceStates, NumBoundStates
!
   real (kind=RealKind), pointer :: dos(:), dos_mt(:), dos_out(:), tps(:)
   real (kind=RealKind) :: sfac, rmul, t0, ssdos
   real (kind=RealKind), pointer :: ps(:), p1(:)
   real (kind=RealKind), pointer :: msgbuf(:,:)
   real (kind=RealKind), allocatable, target :: wks_loc(:)
!
   complex (kind=CmplxKind), intent(out), target :: aux(:)
   complex (kind=CmplxKind) :: energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: der_dos_r_jl(:,:)
!
   type (GridStruct), pointer :: Grid
!
   interface 
      function getFermiDiracFunc(z,mu,kBT) result(fd)
         use KindParamModule, only : IntKind, RealKind
         implicit none
!
         real (kind=RealKind), intent(in) :: mu
         real (kind=RealKind), intent(in) :: kBT
         real (kind=RealKind), intent(in) :: z
         real (kind=RealKind) :: fd
      end function getFermiDiracFunc
   end interface
!
   is = info(1); id = info(2); atom = info(3);  print_dos = info(4); lmax_phi = info(5)
!
   if (id /= local_index) then
      call ErrorHandler('getSScatteringDOS','Inconsistent local index',id,local_index)
   endif
!
   sfac = ONE ! TWO/real(n_spin_pola,kind=RealKind)
!  ===================================================================
!  The following piece of code for multiplying the Fermi-Dirac function
!  needs to be carefully thought since for T = 0, if e > chempot, 
!  getFermiDiracFunc = 0, which causes sfac = 0 and has effects on
!  finding the Fermi energy properly.
!  ===================================================================
   if (Temperature > TEN2m6) then
!     ----------------------------------------------------------------
      sfac = sfac*getFermiDiracFunc(adjustEnergy(is,e),chempot,       &
                                    Boltzmann*Temperature)
!     ----------------------------------------------------------------
   endif
!
   n = getLocalNumSpecies(id)
   allocate(wks_loc(4*n+(4*n+1)*NumPEsInEGroup))
   dos => wks_loc(1:n)
   dos_mt => wks_loc(n+1:2*n)
   dos_out => wks_loc(2*n+1:3*n)
   tps => wks_loc(3*n+1:4*n)
   p1 => wks_loc(4*n+1:4*n+(4*n+1)*NumPEsInEGroup)
   msgbuf => aliasArray2_r(p1,4*n+1,NumPEsInEGroup)
!
   if (present(rfac)) then
      rmul = rfac
   else
      rmul = ONE
   endif
!
   if (present(redundant)) then
      red = redundant
   else
      red = .false.
   endif
!
   energy = adjustEnergy(is,e)
!
   InResonancePeak = .false.
!
   ib = info(6)
   if (ib > 0) then
      if (atom < 1) then
         call ErrorHandler('getSScatteringDOS','Ill condition: info(6) > 0 while atom < 1', &
                           info(6),atom)
      endif
   endif
!
!  ===================================================================
!  For now, I haven't found it useful to calculate the density differently 
!  for those energy points close to the resonance energy, so I will
!  replace the density for those energies with the density calculated for
!  the nearby energy which is outside of the resonance region.
!  -Y.W. Aug. 12, 2021
!  ===================================================================
!  InResonancePeak = .false.; ib = 0
!
   if (ib == 0) then ! In case the density is calculated with the normal scheme
!     write(6,'(a,i3,a,f12.8)')'Not in the resonance peak, print_dos = ',print_dos,', e = ',real(energy,kind=RealKind)
!     if (info(6) > 0) then
!        do while (isEnergyInResonanceRange(e=real(energy,kind=RealKind),site=id,atom=atom,spin=is,res_id=ib))
!           if (real(energy,kind=RealKind) < getResonanceStateEnergy(id,atom,is,ib)) then
!              energy = energy - HALF*resonance_width
!           else
!              energy = energy + HALF*resonance_width
!           endif
!        enddo
!     endif
!     ----------------------------------------------------------------
      call solveSingleScattering(spin=is,site=id,atom=atom,e=energy,vshift=CZERO)
!     ----------------------------------------------------------------
!
      if (getLocalAtomicNumber(id,1) == 0) then
!        call computeDOS(add_highl_fec=.true.)
         call computeDOS(spin=is,site=id,atom=atom)
      else
         call computeDOS(spin=is,site=id,atom=atom)
      endif
!
      Grid => getGrid(id)
!
      n = 0
      do ia = 1, getLocalNumSpecies(id)
         if (atom < 0 .or. ia == atom) then
            dos_r_jl => getDOS(spin=is,site=id,atom=ia)
            if (rad_derivative) then
               der_dos_r_jl => getDOSDerivative(spin=is,site=id,atom=ia)
            endif
!           ==========================================================
!           Becareful that iend may not be equal to NumRs.............
!           ==========================================================
            iend = size(dos_r_jl,1); jmax_dos = size(dos_r_jl,2)
            if (NumRs < iend) then
!              -------------------------------------------------------
               call WarningHandler('getSScatteringDOS','iend > NumRs',iend,NumRs)
!              -------------------------------------------------------
               iend = NumRs
            endif
            if (isASAPotential()) then
               dos(ia) = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
                                                   jmax_dos, 2, dos_r_jl, dos_mt(ia), truncated=.false. )
               dos_mt(ia) = dos(ia)
!              =======================================================
!              Another way to calculate the ASA DOS is using getRadialIntegration
!              -------------------------------------------------------
!              dos_mt(ia) = sfac*getRadialIntegration(id, Grid%jmt, dos_r_jl(:,1))/Y0
!              dos(ia) = dos_mt(ia)
!              =======================================================
               dos_out(ia) = ZERO
            else
!              -------------------------------------------------------
               dos(ia) = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
                                                   jmax_dos, 2, dos_r_jl, dos_mt(ia) )
               dos_mt(ia) = sfac*dos_mt(ia)
               dos_out(ia) = sfac*getOutsideDOS(spin=is,site=id,atom=ia)
!              -------------------------------------------------------
            endif
!
            do jl = 1, jmax_dos
               do ir = 1, iend
                  aux(n+ir) = rmul*sfac*dos_r_jl(ir,jl)
               enddo
!              -------------------------------------------------------
!              call zcopy(iend,dos_r_jl(1,jl),1,aux(n+1),1)
!              -------------------------------------------------------
               n = n + NumRs
            enddo
            if (rad_derivative) then
               do jl = 1, jmax_dos
                  do ir = 1, iend
                     aux(n+ir) = rmul*sfac*der_dos_r_jl(ir,jl)
                  enddo
                  n = n + NumRs
               enddo
            endif
            aux(n+1) = rmul*dos(ia)
            aux(n+2) = rmul*dos_mt(ia)
            aux(n+3) = rmul*dos(ia)*energy
            aux(n+4) = rmul*dos_out(ia)
            n = n + 4
         endif
      enddo
   else ! In this case, InResonancePeak = true, and Lorentzian contribution is excluded
!     ================================================================
!     This routine is placed in a separate module which will be used by
!     SMatrixPolesModule and GFMethodModule.
!     I set an error condition here to make sure ib = 0.
!     ----------------------------------------------------------------
      call ErrorHandler('getSScatteringDOS','Ill condition: info(6) > 0',info(6))
!     ----------------------------------------------------------------
!
!     ================================================================
!     I have commented out, with "!!!", the following lines of the code,
!     implemented previously inside GFMethodModule.
!     02/25/2022 -Yang
!     ================================================================
!     In case of ib = info(6) > 0:
!     info(6) = the resonance state index
!             ! the number of resonance states with the
!             ! width of the resonance peak less than MaxWidth,
!             ! which is specified when calling findSMatrixPoles
!  !  if NumResonanceStates > 0 and e is in the resonance range, for which
!  !            InResonancePeak = .true.
!     it needs to separate the Lorentzian near the resonance energy from the DOS.
!     The Lorentzian contribution to the DOS integration is treated separately.
!     ================================================================
!!!   if (isEnergyInResonanceRange(e=real(energy,kind=RealKind),site=id,atom=atom,spin=is,res_id=ib0)) then
!!!      InResonancePeak = .true.
!!!      if (ib0 /= ib) then
!!!         call ErrorHandler('getSScatteringDOS','resonance state ib0 <> ib',ib0,ib)
!!!      endif
!!!   else
!!!      call ErrorHandler('getSScatteringDOS','ib > 0, while e is not at a resonance energy',e)
!!!   endif
!     if (abs(energy-getResonanceStateEnergy(id,atom,is,1)) < 0.005) then
!        write(6,'(a,2f12.8,a,i3)')'e, re = ',real(energy,kind=RealKind),getResonanceStateEnergy(id,atom,is,1), &
!                                  ', ib = ',ib
!     endif
!!!   write(6,'(a,i3,a,f12.8)')'At the resonance peak, ib = ',ib,', e = ',real(energy,kind=RealKind)
!!!   n = 0
!!!   do ia = 1, getLocalNumSpecies(id)
!!!      if (atom < 0 .or. ia == atom) then
!!!         if (rad_derivative) then
!!!            dos_r_jl => getResidualResStateDensity(site=id,atom=ia,spin=is,rstate=ib, &
!!!                                                   e=real(energy,kind=RealKind),      &
!!!                                                   NumRs=iend,jmax_rho=jmax_dos,      &
!!!                                                   derivative=der_dos_r_jl)
!!!         else
!!!            dos_r_jl => getResidualResStateDensity(site=id,atom=ia,spin=is,rstate=ib, &
!!!                                                   e=real(energy,kind=RealKind),      &
!!!                                                   NumRs=iend,jmax_rho=jmax_dos)
!!!         endif
!
!!!         if (isASAPotential()) then
!!!            dos(ia) = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
!!!                                                 jmax_dos, 2, dos_r_jl, dos_mt(ia), truncated=.false. )
!!!            dos_mt(ia) = dos(ia)
!!!            dos_out(ia) = ZERO
!!!         else
!!!            dos(ia) = sfac*getVolumeIntegration( id, iend, Grid%r_mesh, jmax_dos, 2, dos_r_jl, dos_mt(ia) )
!!!            dos_mt(ia) = sfac*dos_mt(ia)
!!!            dos_out(ia) = sfac*getOutsideDOS(spin=is,site=id,atom=ia)
!!!         endif
!
!!!         do jl = 1, jmax_dos
!!!            do ir = 1, iend
!!!               aux(n+ir) = rmul*sfac*dos_r_jl(ir,jl)
!!!            enddo
!              -------------------------------------------------------
!              call zcopy(NumRs,dos_r_jl(1,jl),1,aux(n+1),1)
!              -------------------------------------------------------
!!!            n = n + NumRs
!!!         enddo
!!!         if (rad_derivative) then
!!!            do jl = 1, jmax_dos
!!!               do ir = 1, iend
!!!                  aux(n+ir) = rmul*sfac*der_dos_r_jl(ir,jl)
!!!               enddo
!!!               n = n + NumRs
!!!            enddo
!!!         endif
!!!         aux(n+1) = rmul*dos(ia)
!!!         aux(n+2) = rmul*dos_mt(ia)
!!!         aux(n+3) = rmul*dos(ia)*energy
!!!         aux(n+4) = rmul*dos_out(ia)
!!!         n = n + 4
!!!      endif
!!!   enddo
   endif
!
!  -------------------------------------------------------------------
   call computePhaseShift(spin=is,site=id,atom=atom)
!  -------------------------------------------------------------------
   kmax_phi = (lmax_phi+1)**2
   tps = ZERO
   do ia = 1, getLocalNumSpecies(id)
      if (atom < 0 .or. ia == atom) then
         ps => getPhaseShift(spin=is,site=id,atom=ia)
         do kl = 1, kmax_phi
            tps(ia) = tps(ia) + ps(kl)
         enddo
      endif
   enddo
!
   if (print_dos > 0) then
      if (.not.red) then
         msgbuf = ZERO
         msgbuf(1,MyPEinEGroup+1) = real(energy,kind=RealKind)
         n = 1
         do ia = 1, getLocalNumSpecies(id)
            if (atom < 0 .or. ia == atom) then
               msgbuf(n+1,MyPEinEGroup+1) = dos(ia)
               msgbuf(n+2,MyPEinEGroup+1) = dos_mt(ia)
               msgbuf(n+3,MyPEinEGroup+1) = dos_out(ia)
               msgbuf(n+4,MyPEinEGroup+1) = tps(ia)
               n = n + 4
            endif
         enddo
!        -------------------------------------------------------------
         call GlobalSumInGroup(eGID,msgbuf,4*getLocalNumSpecies(id)+1,&
                                  NumPEsInEGroup)
!        -------------------------------------------------------------
         if ( node_print_level >= 0) then
            if (getLocalNumSpecies(id) == 1) then
               do i = 1, NumPEsInEGroup
!                 write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')real(energy,kind=RealKind),dos,dos_mt,dos_out,tps
                  write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')msgbuf(1:5,i)
               enddo
            else
               do i = 1, NumPEsInEGroup
!                 write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')real(energy,kind=RealKind),dos,dos_mt,dos_out,tps
                  n = 1
                  do ia = 1, getLocalNumSpecies(id)
                     if (atom < 0 .or. ia == atom) then
                        write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8),4x,a)')msgbuf(1,i),msgbuf(n+1,i), &
                              msgbuf(n+2,i), msgbuf(n+3,i), msgbuf(n+4,i), getLocalAtomName(id,ia)
                        n = n + 4
                     endif
                  enddo
               enddo
            endif
         endif
      else if ( node_print_level >= 0) then
         if (getLocalNumSpecies(id) == 1) then
            write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')real(energy,kind=RealKind),dos(1),dos_mt(1), &
                  dos_out(1),tps(1)
         else
            do ia = 1, getLocalNumSpecies(id)
               if (atom < 0 .or. ia == atom) then
                  write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8),4x,a)')real(energy,kind=RealKind),dos(ia),dos_mt(ia), &
                        dos_out(ia),tps(ia), getLocalAtomName(id,ia)
               endif
            enddo
         endif
      endif
   endif
!
   if (atom < 1) then
      ssdos = dos(1)
   else
      ssdos = dos(atom)
   endif
!
   nullify(msgbuf, dos, dos_mt, dos_out, tps)
   deallocate(wks_loc)
!
   end function getSScatteringDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSScatteringDOSofCmplxE(info,e,aux,cfac) result(dos)
!  ===================================================================
!
!  This function returns single site DOS for a given energy in the complex 
!  plane. The returned dos value includes a factor cfac if it is present.
!
!  The function also returns aux, a data set that can be integrated
!  on an energy grid
!
!  In the no-spin-polarized case, the returned dos is per spin.
!  ===================================================================
   use PotentialTypeModule, only : isASAPotential, isFullPotential
!
   use RadialGridModule, only : getGrid, getRadialIntegration
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use SSSolverModule, only : solveSingleScattering
   use SSSolverModule, only : computeGreenFunction, getGreenFunction
   use SSSolverModule, only : getGreenFunctionDerivative
!
   use AtomModule, only : getLocalNumSpecies
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), intent(in), optional :: cfac
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: is, id, atom, ia, lmax_phi
   integer (kind=IntKind) :: kmax, jmax, iend, n, ir, jl, l, m, kl, klc
!
   real (kind=RealKind) :: dos, sfac, t0
!
   complex (kind=CmplxKind), intent(out), target :: aux(:)
   complex (kind=CmplxKind) :: energy, cmul, greenint, greenint_mt, ede
   complex (kind=CmplxKind), pointer :: green(:,:)
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: der_green(:,:), der_dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: pcv_x(:), pcv_y(:), pcv_z(:), p1(:)
!
   type (GridStruct), pointer :: Grid
!
!
   is = info(1); id = info(2); atom = info(3); lmax_phi = info(5)
!
   if (id /= local_index) then
      call ErrorHandler('getSScatteringDOSofCmplxE','Inconsistent local index',id,local_index)
   endif
!
   sfac= ONE ! TWO/real(n_spin_pola,kind=RealKind)
!
   if (present(cfac)) then
      cmul = cfac
   else
      cmul = CONE
   endif
!
!  ===================================================================
!  For energy in the upper complex plane, use spherical solver
!  together with irregular solution
!  ===================================================================
   energy = adjustEnergy(is,e)
!
!! -------------------------------------------------------------------
 ! call solveSingleScattering(spin=is,site=id,atom=atom,e=energy,vshift=CZERO,    &
 !                            isSphSolver=.true.,useIrrSol='h')
!! -------------------------------------------------------------------
   if (isFullPotential()) then
!     ================================================================
!     Note: In the full-potential case, isSphSolver is disabled so that
!           the complete green function is not calculated. We need to be 
!           careful with the pole of the sine matrix.
!     ----------------------------------------------------------------
      call solveSingleScattering(spin=is,site=id,atom=atom,e=energy,vshift=CZERO, &
                                 isSphSolver=.false.)
!     ----------------------------------------------------------------
      call computeGreenFunction(spin=is,site=id,atom=atom,noIrrTerm=.true.)
!     ----------------------------------------------------------------
   else
!     ================================================================
!     Note: If isSphSolver is enabled, unlike the full-potential case,
!           useIrrSol can be set so that complete green function is
!           calculated.
!     ----------------------------------------------------------------
      call solveSingleScattering(spin=is,site=id,atom=atom,e=energy,vshift=CZERO, &
                                 isSphSolver=.true.,useIrrSol='h')
      call computeGreenFunction(spin=is,site=id,atom=atom)
!     ----------------------------------------------------------------
   endif
!
!  -------------------------------------------------------------------
!! call computeGreenFunction(spin=is,site=id,atom=atom)
!  -------------------------------------------------------------------
!
   Grid => getGrid(id)
   n = 0
   do ia = 1, getLocalNumSpecies(id)
      if (atom < 0 .or. ia == atom) then
!        -------------------------------------------------------------
         green=>getGreenFunction(spin=is,site=id,atom=ia)
!        -------------------------------------------------------------
!
!        =============================================================
!        Becareful that iend may not be equal to NumRs................
!        =============================================================
         iend = size(green,1); kmax = size(green,2); jmax = jofk(kmax)
         if (jmax /= jmax_green) then
!           ----------------------------------------------------------
            call ErrorHandler('getSScatteringDOSofCmplxE',            &
                              'Inconsistent jmax value',jofk(kmax),jmax_green)
!           ----------------------------------------------------------
         else if (NumRs < iend) then
!           ----------------------------------------------------------
            call WarningHandler('getSScatteringDOSofCmplxE','iend > NumRs',iend,NumRs)
!           ----------------------------------------------------------
            iend = NumRs
         endif
         if (isASAPotential()) then
            greenint_mt = sfac*getRadialIntegration(id, Grid%jmt, green(:,1))/Y0
            greenint = greenint_mt
!           Alternatively, you can call getVolumeIntegration with truncated=.false.
!           to run the integration in ASA volume.
         else
!           ----------------------------------------------------------
            greenint = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,    &
                                                  kmax, 2, green, greenint_mt )
!           ----------------------------------------------------------
            greenint_mt = sfac*greenint_mt
         endif
!
         if (atom > 0 .or. ia == 1) then
            dos = real(SQRTm1*cmul*greenint/PI,kind=RealKind)
         endif
!
         if (iharris <= 1) then
            ede = energy
         else
            ede = (energy-chempot)
         endif
!
         p1 => aux(n+1:n+NumRs*jmax)
         dos_r_jl => aliasArray2_c(p1, NumRs, jmax)
         do jl = 1, jmax
            l = lofj(jl); m = mofj(jl)
            kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
            pcv_x => green(1:iend,kl)
            pcv_y => green(1:iend,klc)
            pcv_z => dos_r_jl(1:iend,jl)
            pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
         enddo
         n = n + NumRs*jmax
         if (rad_derivative) then
            der_green=>getGreenFunctionDerivative(spin=is,site=id,atom=ia)
            p1 => aux(n+1:n+NumRs*jmax)
            der_dos_r_jl => aliasArray2_c(p1, NumRs, jmax)
            do jl = 1, jmax
               l = lofj(jl); m = mofj(jl)
               kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
               pcv_x => der_green(1:iend,kl)
               pcv_y => der_green(1:iend,klc)
               pcv_z => der_dos_r_jl(1:iend,jl)
               pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
            enddo
            n = n + NumRs*jmax
         endif
         aux(n+1) = SQRTm1*cmul*greenint/PI
         aux(n+2) = SQRTm1*cmul*greenint_mt/PI
         aux(n+3) = SQRTm1*cmul*ede*greenint/PI
!        -------------------------------------------------------------
         greenint = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,       &
                                               kmax, 2, green, greenint_mt, &
                                               truncated=.false. ) - greenint
!        -------------------------------------------------------------
         aux(n+4) = SQRTm1*greenint/PI
         if (node_print_level >= 1) then
            write(6,'(f12.8,3x,d15.8,2x,2(4x,d15.8))')real(energy,kind=RealKind), &
                                                      real(aux(n+1),RealKind),    &
                                                      real(aux(n+2),RealKind),    &
                                                      real(aux(n+4),RealKind)
         endif
         n = n + 4
      endif
   enddo
!
   end function getSScatteringDOSofCmplxE
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRelSScatteringDOS(info,e,aux,rfac,redundant) result(dos)
!  ===================================================================
!  The following are subroutines by xianglin for relativistic calculation
!
!  This function returns single site DOS for a given energy on the real 
!  energy axis.
!
!  The function also returns aux, a data set that can be integrated
!  on an energy grid
!  ===================================================================
   use PhysParamModule, only : Boltzmann
   use PotentialTypeModule, only : isASAPotential
   use GroupCommModule, only : GlobalSumInGroup
   use RadialGridModule, only : getGrid
   use StepFunctionModule, only : getVolumeIntegration
   use RelSSSolverModule, only : getRelSSDOS, SingleDiracScattering
!
   implicit none
!
   real (kind=RealKind), intent(in) :: e
   real (kind=RealKind), intent(in), optional :: rfac
!
   logical, intent(in), optional :: redundant
   logical :: red
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: jmax_dos, kmax_phi, iend, n, kl, jl, ir, n0
   integer (kind=IntKind) :: is, id, print_dos, atom, lmax_phi
!
   real (kind=RealKind) :: dos, sfac, dos_mt, dos_out, tps, rmul, dos1, t0
   real (kind=RealKind), pointer :: ps(:)
   real (kind=RealKind) :: msgbuf(5,NumPEsInEGroup)
!
   complex (kind=CmplxKind), intent(out), target :: aux(:)
   complex (kind=CmplxKind) :: energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
!
   type (GridStruct), pointer :: Grid
!
   interface 
      function getFermiDiracFunc(z,mu,kBT) result(fd)
         use KindParamModule, only : IntKind, RealKind
         implicit none
!
         real (kind=RealKind), intent(in) :: mu
         real (kind=RealKind), intent(in) :: kBT
         real (kind=RealKind), intent(in) :: z
         real (kind=RealKind) :: fd
      end function getFermiDiracFunc
   end interface
!
   is = info(1); id = info(2); atom = info(3); print_dos = info(4)
!
   if (id /= local_index) then
      call ErrorHandler('getRelSScatteringDOS','Inconsistent local index',id,local_index)
   endif
!
   lmax_phi = info(5)
!
   sfac= ONE !TWO/real(n_spin_pola,kind=RealKind)
   sfac = sfac*getFermiDiracFunc(e,chempot,Boltzmann*Temperature)
!
   if (present(rfac)) then
      rmul = rfac
   else
      rmul = ONE
   endif
!
   if (present(redundant)) then
      red = redundant
   else
      red = .false.
   endif
!
!   energy = adjustEnergy(is,e)
   energy = e
!
!  t0 = getTime()
!  -------------------------------------------------------------------
   call SingleDiracScattering(id,energy)
!  -------------------------------------------------------------------
!  Timing_SS = Timing_SS + (getTime() - t0)
!  NumCalls_SS = NumCalls_SS + 1
!
!   if (AtomicNumber(id) == 0) then
!!     call computeDOS(add_highl_fec=.true.)
!      call computeDOS()
!   else
!      call computeDOS()
!   endif
!  -------------------------------------------------------------------
!   dos_out = sfac*getOutsideDOS()
!  commented since no phase shift is calculated relativistically right now
!  -------------------------------------------------------------------
!   call computePhaseShift()
!   ps => getPhaseShift()
!   kmax_phi = (lmax_phi+1)**2
!   tps = ZERO
!   do kl = 1, kmax_phi
!      tps = tps + ps(kl)
!   enddo
!!
!   if (print_dos > 0) then
!      msgbuf = ZERO
!      msgbuf(1,MyPEinEGroup+1) = real(energy,kind=RealKind)
!      msgbuf(2,MyPEinEGroup+1) = dos
!      msgbuf(3,MyPEinEGroup+1) = dos_mt
!      msgbuf(4,MyPEinEGroup+1) = dos_out
!      msgbuf(5,MyPEinEGroup+1) = tps
!!     ----------------------------------------------------------------
!      call GlobalSumInGroup(eGID,msgbuf,5,NumPEsInEGroup)
!!     ----------------------------------------------------------------
!      if ( node_print_level >= 0) then
!         do n = 1, NumPEsInEGroup
!!           write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')real(energy,kind=RealKind),dos,dos_mt,dos_out,tps
!            write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')msgbuf(1:5,n)
!         enddo
!      endif
!   endif
!
   dos_r_jl => getRelSSDOS(1,id) !is 1:4
   iend = size(dos_r_jl,1); jmax_dos = size(dos_r_jl,2)
   Grid => getGrid(id)
   do is = 1, n_spin_cant*n_spin_cant
      n0=(is-1)*(NumRs*jmax_dos+4)
      dos_r_jl => getRelSSDOS(is,id) !is 1:4
      do jl = 1, jmax_dos
         n = (jl-1)*NumRs
         do ir = 1, NumRs
            aux(n0+n+ir) = rmul*sfac*dos_r_jl(ir,jl)
         enddo
!        ----------------------------------------------------------------
!        call zcopy(NumRs,dos_r_jl(1,jl),1,aux(n+1),1)
!        ----------------------------------------------------------------
      enddo
      n = NumRs*jmax_dos
      if (isASAPotential()) then
         dos = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
                                          jmax_dos, 2, dos_r_jl, dos_mt, truncated=.false. )
      else
         dos = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
                                          jmax_dos, 2, dos_r_jl, dos_mt )
      endif
      dos_mt = sfac*dos_mt
!      print*,"single-site dos at real axis", "  is=", is
!      print*,"energy=",real(energy,kind=RealKind), "dos=", dos, "dos_mt=", dos_mt
      aux(n0+n+1) = rmul*dos
      aux(n0+n+2) = rmul*dos_mt
      aux(n0+n+3) = rmul*dos*energy
      aux(n0+n+4) = CZERO!rmul*dos_out
      if (is==1) then
         dos1=dos
      endif
   enddo
!
   dos=dos1
   end function getRelSScatteringDOS
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRelSScatteringDOSatPole(info,e,aux,rfac) result(dos)
!  ===================================================================
!  The following are subroutines by xianglin for relativistic calculation
!  This function is not used because the DOS calculated is unstable
!
!  This function returns single site DOS for a given energy on the real 
!  energy axis.
!
!  The function also returns aux, a data set that can be integrated
!  on an energy grid
!  ===================================================================
   use GroupCommModule, only : GlobalSumInGroup
   use PublicTypeDefinitionsModule, only : GridStruct
   use PotentialTypeModule, only : isASAPotential
   use RadialGridModule, only : getGrid
   use StepFunctionModule, only : getVolumeIntegration
   use RelSSSolverModule, only :getRelSSDOS,getRelSSGreen, getRelSSGreenOut,&
                SingleDiracScattering
!
   implicit none
!
   real (kind=RealKind), intent(in) :: e
   real (kind=RealKind), intent(in), optional :: rfac
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: jmax_dos, kmax_phi, iend, n, kl, jl, ir, n0
   integer (kind=IntKind) :: is, id, print_dos, atom, lmax_phi
!
   real (kind=RealKind) :: dos, sfac, dos_mt, dos_out, tps, rmul, dos1, dos_end, t0
   real (kind=RealKind), pointer :: ps(:)
   real (kind=RealKind) :: msgbuf(5,NumPEsInEGroup)
!
   complex (kind=CmplxKind), intent(out), target :: aux(:)
   complex (kind=CmplxKind) :: energy
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
!
   type (GridStruct), pointer :: Grid
!
   is = info(1); id = info(2); atom = info(3); print_dos = info(4)
!
   if (id /= local_index) then
      call ErrorHandler('getRelSScatteringDOSatPole','Inconsistent local index',id,local_index)
   endif
!
   lmax_phi = info(5)
   sfac= ONE !TWO/real(n_spin_pola,kind=RealKind)
!
!   if (present(rfac)) then
!      rmul = rfac
!   else
!      rmul = ONE
!   endif
!
!   energy = adjustEnergy(is,e)
!   energy = cmplx(e,0.d0)!cmplx(e,TEN2m6,kind=CmplxKind)
    energy = cmplx(e,TEN2m6,kind=CmplxKind)
!
!  t0 = getTime()
!  -------------------------------------------------------------------
   call SingleDiracScattering(id,energy)
!  -------------------------------------------------------------------
!  Timing_SS = Timing_SS + (getTime() - t0)
!  NumCalls_SS = NumCalls_SS + 1
!
!   if (AtomicNumber(id) == 0) then
!!     call computeDOS(add_highl_fec=.true.)
!      call computeDOS()
!   else
!      call computeDOS()
!   endif
!  -------------------------------------------------------------------
!   dos_out = sfac*getOutsideDOS()
!  commented since no phase shift is calculated relativistically right now
!  -------------------------------------------------------------------
!   call computePhaseShift()
!   ps => getPhaseShift()
!   kmax_phi = (lmax_phi+1)**2
!   tps = ZERO
!   do kl = 1, kmax_phi
!      tps = tps + ps(kl)
!   enddo
!!
!   if (print_dos > 0) then
!      msgbuf = ZERO
!      msgbuf(1,MyPEinEGroup+1) = real(energy,kind=RealKind)
!      msgbuf(2,MyPEinEGroup+1) = dos
!      msgbuf(3,MyPEinEGroup+1) = dos_mt
!      msgbuf(4,MyPEinEGroup+1) = dos_out
!      msgbuf(5,MyPEinEGroup+1) = tps
!!     ----------------------------------------------------------------
!      call GlobalSumInGroup(eGID,msgbuf,5,NumPEsInEGroup)
!!     ----------------------------------------------------------------
!      if ( node_print_level >= 0) then
!         do n = 1, NumPEsInEGroup
!!           write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')real(energy,kind=RealKind),dos,dos_mt,dos_out,tps
!            write(6,'(f12.8,3x,d15.8,2x,3(4x,d15.8))')msgbuf(1:5,n)
!         enddo
!      endif
!   endif
!
   dos_r_jl => getRelSSDOS(1,id) !is 1:4
   iend = size(dos_r_jl,1); jmax_dos = size(dos_r_jl,2)
   Grid => getGrid(id)
!-------------------------------------------------
!  find the correct rmul
!-------------------------------------------------
   dos = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
         jmax_dos, 2, dos_r_jl, dos_mt, truncated=.false. )
!
!  dos_out = int_0^infty A*r^2*exp(-2*k*x)*Y_00; dos(iend,1)=A*r^2 not used since it's wrong
   dos_end = real(dos_r_jl(iend,1),RealKind)
   if (e > 0) then
      call ErrorHandler("getRelSScatteringDOSatPole","unimplemented for poles on positive axis")
   endif
   dos_out = -AIMAG(getRelSSGreenOut(id))
!  SQRT(4.d0*PI)*dos_end/( 4.d0*(SQRT(-e))**3*(Grid%r_mesh(iend))**2 )
   rmul = 1.d0/(dos + dos_out)
   if (MyPE == 0) then
      print*, "dos=", dos, " dos_out=",dos_out, getRelSSGreenOut(id), " dos_total", dos+dos_out
      print*,"e=", e," number of electrons in this resonance: ", dos/(dos+dos_out)
   endif
!
   do is = 1, n_spin_cant*n_spin_cant
      n0=(is-1)*(NumRs*jmax_dos+4)
      dos_r_jl => getRelSSDOS(is,id) !is 1:4
      do jl = 1, jmax_dos
         n = (jl-1)*NumRs
         do ir = 1, NumRs
            aux(n0+n+ir) = rmul*sfac*dos_r_jl(ir,jl)
         enddo
!     ----------------------------------------------------------------
!     call zcopy(NumRs,dos_r_jl(1,jl),1,aux(n+1),1)
!     ----------------------------------------------------------------
      enddo
      n = NumRs*jmax_dos
      if (isASAPotential()) then
         dos = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
                                          jmax_dos, 2, dos_r_jl, dos_mt, truncated = .false. )
      else
         dos = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,   &
                                          jmax_dos, 2, dos_r_jl, dos_mt )
      endif
      dos_mt = sfac*dos_mt
      aux(n0+n+1) = rmul*dos
      aux(n0+n+2) = rmul*dos_mt
      aux(n0+n+3) = rmul*dos*real(energy,RealKind)
      aux(n0+n+4) = CZERO!rmul*dos_out
      if (is==1) then
         dos1=dos
         if (MyPE == 0) then
            print*,"aux(n0+n+1)=",aux(n0+n+1)
         endif
      endif
   enddo
!
   dos=dos1
!   if (MyPE == 0 .and. e<-.4 ) then
!      open (104,file="dos_r_jl1",action="write")
!      do jl=1,jmax_dos
!         do ir=1,iend
!            write(104,*) Grid%r_mesh(ir),real(dos_r_jl(ir,jl)*rmul,kind=RealKind)
!         enddo
!      enddo
!      close (104)
!   else
!      open (105,file="dos_r_jl2",action="write")
!      do jl=1,jmax_dos
!         do ir=1,iend
!            write(105,*) Grid%r_mesh(ir),real(dos_r_jl(ir,jl)*rmul,kind=RealKind)
!         enddo
!      enddo
!      close (105)
!   endif
   end function getRelSScatteringDOSatPole
!
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRelSScatteringDOSofCmplxE(info,e,aux,cfac) result(dos)
!  ===================================================================
!  The following are subroutines by xianglin for relativistic calculation
!
!  This function returns single site DOS for a given energy in the complex 
!  plane. The returned dos value includes a factor cfac if it is present.
!
!  The function also returns aux, a data set that can be integrated
!  on an energy grid
!  ===================================================================
   use PotentialTypeModule, only : isASAPotential
!
   use RadialGridModule, only : getGrid, getRadialIntegration
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use RelSSSolverModule, only : getRelSSGreen, SingleDiracScattering
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), intent(in), optional :: cfac
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: is, id, ks, atom
   integer (kind=IntKind) :: kmax, jmax, iend, n, ir, jl, l, m, kl, klc, p0, n0
!
   real (kind=RealKind) :: dos, sfac, t0
!
   complex (kind=CmplxKind), intent(out), target :: aux(:)
   complex (kind=CmplxKind) :: energy, cmul, greenint, greenint_mt, ede
   complex (kind=CmplxKind), pointer :: green(:,:)
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: pcv_x(:), pcv_y(:), pcv_z(:), p1(:)
!
   type (GridStruct), pointer :: Grid
!
   is = info(1); id = info(2); atom = info(3)
!
   if (id /= local_index) then
      call ErrorHandler('getRelSScatteringDOSofCmplxE','Inconsistent local index',id,local_index)
   endif
!
   sfac= ONE !TWO/real(n_spin_pola,kind=RealKind)
!
   if (present(cfac)) then
      cmul = cfac
   else
      cmul = CONE
   endif
!
!  ===================================================================
!  For energy in the upper complex plane, use spherical solver
!  together with irregular solution
!  ===================================================================
   energy = adjustEnergy(is,e)
!
!  t0 = getTime()
!  -------------------------------------------------------------------
   call SingleDiracScattering(id,energy)
!  -------------------------------------------------------------------
!  Timing_SS = Timing_SS + (getTime() - t0)
!  NumCalls_SS = NumCalls_SS + 1
!
!  -------------------------------------------------------------------
   green=>getRelSSGreen(1,id)
!  -------------------------------------------------------------------
!
   iend = size(green,1); kmax = size(green,2); jmax = jofk(kmax)
   Grid => getGrid(id)
   if (isASAPotential()) then
      greenint_mt = sfac*getRadialIntegration(id, Grid%jmt, green(:,1))/Y0
      greenint = greenint_mt
   else
!     ----------------------------------------------------------------
      greenint = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,            &
                                            kmax, 2, green, greenint_mt )
!     ----------------------------------------------------------------
      greenint_mt = sfac*greenint_mt
   endif
!
   dos = real(SQRTm1*cmul*greenint/PI,kind=RealKind)
!
   if (iharris <= 1) then
      ede = energy
   else
      ede = (energy-chempot)
   endif
!
   n = iend*jmax
   do ks = 1, n_spin_cant*n_spin_cant
      n0=(ks-1)*(iend*jmax+4)
      green => getRelSSGreen(ks,id)
      p1 => aux(n0+1:n0+n)
      dos_r_jl => aliasArray2_c(p1,iend,jmax)
      do jl = 1, jmax
         l = lofj(jl); m = mofj(jl)
         kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
         pcv_x => green(1:iend,kl)
         pcv_y => green(1:iend,klc)
         pcv_z => dos_r_jl(1:iend,jl)
         pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
      enddo
   enddo
   aux(n0+n+1) = SQRTm1*cmul*greenint/PI
   aux(n0+n+2) = SQRTm1*cmul*greenint_mt/PI
   aux(n0+n+3) = SQRTm1*cmul*ede*greenint/PI
!  -------------------------------------------------------------------
   greenint = sfac*getVolumeIntegration( id, iend, Grid%r_mesh,       &
                                         kmax, 2, green, greenint_mt, &
                                         truncated=.false. ) - greenint
!  -------------------------------------------------------------------
   aux(n+4) = SQRTm1*cmul*greenint/PI
!  print*,"aux(n0+n+1)=",aux(n0+n+1)
!  print*,"aux(n0+n+2)=",aux(n0+n+2)
!  print*,"aux(n0+n+3)=",aux(n0+n+3)
!  print*,"aux(n0+n+4)=",aux(n0+n+4)
!
   end function getRelSScatteringDOSofCmplxE
!  ===================================================================
end module SingleScatteringDOSModule
