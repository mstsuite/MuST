Module ConvergenceCheckModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, CZERO, PI4, THIRD, HALF, Ten2m14
   use PublicTypeDefinitionsModule, only : GridStruct
!
public :: initConvergenceCheck,     &
          endConvergenceCheck,      &
          checkConvergence
!
private
   logical :: Initialized = .false.
   logical :: isFullPot = .false.
!
   integer (kind=IntKind) :: iteration
   integer (kind=IntKind) :: LocalNumAtoms
   integer (kind=IntKind) :: n_spin_pola
   integer (kind=IntKind) :: nr_max, jmax, lmax
   integer (kind=IntKind), allocatable :: Print_Level(:)
!
   real (kind=RealKind) :: Old_FermiEnergy
   real (kind=RealKind) :: Old_TotalEnergy
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initConvergenceCheck(nlocal,pola,iprint)
!  ===================================================================
   use Atom2ProcModule, only : getGlobalIndex
   use SystemModule, only : getLmaxRho, getLmaxPot
   use PotentialTypeModule, only : isFullPotential
!
   use RadialGridModule, only : getNumRmesh, getGrid
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: nlocal, pola
   integer (kind=IntKind), intent(in) :: iprint(nlocal)
!
   integer (kind=IntKind) :: id, gid
!
   LocalNumAtoms = nlocal
   n_spin_pola = pola
   allocate(Print_Level(nlocal))
   Print_Level(1:nlocal) = iprint(1:nlocal)
!
   Old_FermiEnergy = ZERO
   Old_TotalEnergy = ZERO
!
   iteration = 0
   isFullPot = isFullPotential()
!
   if ( .not.isFullPot ) then
      nr_max = 0
      do id = 1, LocalNumAtoms
         nr_max = max(nr_max, getNumRmesh(id))
      enddo
   else
      lmax = 0
      nr_max = 0
      do id = 1,LocalNumAtoms
         gid = getGlobalIndex(id)
         nr_max = max(nr_max, getNumRmesh(id))
         lmax = max(lmax,getLmaxRho(gid))
         lmax = max(lmax,getLmaxPot(gid))
      enddo
      jmax = (lmax+1)*(lmax+2)/2
   endif
!
   Initialized = .true.
!
   end subroutine initConvergenceCheck
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endConvergenceCheck()
!  ===================================================================
   implicit none
!
   iteration = 0
   Initialized = .false.
!
   end subroutine endConvergenceCheck
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine checkConvergence( rho_rms, pot_rms, ef_diff, et_diff,   &
                                evec_rms, bcon_rms, efermi, etot )
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
!
   use Atom2ProcModule, only : getGlobalIndex
!
   use DataServiceCenterModule, only : getDataStorage,        &
                                       isDataStorageExisting, &
                                       RealType, ComplexType, &
                                       RealMark, ComplexMark
!
   use RadialGridModule, only : getNumRmesh, getGrid
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use ChargeDensityModule, only : getRhoLmax
   use ChargeDensityModule, only : getSphChargeDensity
   use ChargeDensityModule, only : getChargeDensity
   use ChargeDensityModule, only : getSphMomentDensity
   use ChargeDensityModule, only : getMomentDensity
   use ChargeDensityModule, only : getChargeComponentFlag
!
   use AtomModule, only : getLocalNumSpecies, getLocalSpeciesContent
!
   use PotentialModule, only : getPotLmax
   use PotentialModule, only : getOldSphPotr => getSphPotr
   use PotentialModule, only : getOldPotential => getPotential
!
   use PotentialGenerationModule, only : getNewSphPotr => getSphPotr
   use PotentialGenerationModule, only : getNewPotential => getPotential
   use PotentialGenerationModule, only : getPotComponentFlag
!
   use PolyhedraModule, only : getVolume
!
   implicit none
!
   integer (kind=IntKind) :: is, ir, jmt, id, jl, jend, ia, n
   integer (kind=IntKind) :: lmax_l, jmax_l, kmax_l
   integer (kind=IntKind) :: lmax_prod, jmax_prod, kmax_prod
   integer (kind=IntKind), pointer :: flag_jl(:)
!
   real (kind=RealKind), intent(in) :: efermi, etot
   real (kind=RealKind), intent(out) :: rho_rms(:,:)
   real (kind=RealKind), intent(out) :: pot_rms(:,:)
   real (kind=RealKind), intent(out) :: evec_rms(LocalNumAtoms)
   real (kind=RealKind), intent(out) :: bcon_rms(LocalNumAtoms)
   real (kind=RealKind), intent(out) :: ef_diff, et_diff
!
   real (kind=RealKind), pointer :: vec1(:)
   real (kind=RealKind), pointer :: vec2(:)
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), allocatable :: wk(:)
   real (kind=RealKind) :: mt_volume, vol_int, VP_volume, tol, vint_mt
   real (kind=RealKind) :: cfac, rho_rms_av(n_spin_pola), pot_rms_av(n_spin_pola)
!
   complex (kind=CmplxKind), pointer :: p_den1(:,:)
   complex (kind=CmplxKind), pointer :: p_den2(:,:)
   complex (kind=CmplxKind), pointer :: p_wk(:,:), prod(:,:)
!
   complex (kind=CmplxKind), allocatable, target :: wk_c(:), wk_prod(:)
!
   type (GridStruct), pointer :: Grid
!
   interface
      subroutine computeProdExpan(n,lf,f,lg,g,lh,h)
        use KindParamModule, only : IntKind, CmplxKind
        integer (kind=IntKind), intent(in) :: n, lf, lg, lh
        complex (kind=CmplxKind), intent(in) :: f(:,:), g(:,:)
        complex (kind=CmplxKind), intent(out) :: h(:,:)
      end subroutine computeProdExpan 
   end interface
!
   if (.not.Initialized) then
      call ErrorHandler('checkConvergence','Module is not initialized')
   endif
!
   if (iteration == 0) then
      ef_diff = ZERO
      et_diff = ZERO
      evec_rms = ZERO
      bcon_rms = ZERO
   else
      ef_diff = efermi - Old_FermiEnergy;    ef_diff = abs(ef_diff)
      et_diff = etot - Old_TotalEnergy;  et_diff = abs(et_diff)
   endif
!
   Old_FermiEnergy = efermi
   Old_TotalEnergy = etot
!
   iteration = iteration + 1
!
   if ( .not.isFullPot ) then
      allocate( wk(nr_max) )
      n = 0
      rho_rms = ZERO
      pot_rms = ZERO
      do id = 1, LocalNumAtoms
         Grid => getGrid(id)
         jmt = Grid%jmt
         r_mesh => Grid%r_mesh(1:jmt)
         mt_volume = PI4*(Grid%rmt)**3*THIRD
!
!        =============================================================
!        calculate RMS of density.
!        =============================================================
         rho_rms_av = ZERO
         pot_rms_av = ZERO
         do ia = 1, getLocalNumSpecies(id)
            n = n + 1
            cfac = getLocalSpeciesContent(id,ia)
            if (n_spin_pola == 1) then
               vec1 => getSphChargeDensity("TotalOld",id,ia)
               vec2 => getSphChargeDensity("TotalNew",id,ia)
               do ir=1,jmt
                  wk(ir)=(vec1(ir)-vec2(ir))**2
               enddo
               vol_int = getVolumeIntegration( id, jmt, r_mesh, 0, wk )
               rho_rms(1,n) = sqrt(vol_int)/mt_volume
               rho_rms_av(1) = rho_rms_av(1) + cfac*rho_rms(1,n)
            else
               vec1 => getSphChargeDensity("TotalOld",id,ia)
               vec2 => getSphChargeDensity("TotalNew",id,ia)
               do ir=1,jmt
                  wk(ir)=(vec1(ir)-vec2(ir))**2
               enddo
               vol_int = getVolumeIntegration( id, jmt, r_mesh, 0, wk )
               rho_rms(1,n) = sqrt(vol_int)/mt_volume
               rho_rms_av(1) = rho_rms_av(1) + cfac*rho_rms(1,n)
               vec1 => getSphMomentDensity("TotalOld",id,ia)
               vec2 => getSphMomentDensity("TotalNew",id,ia)
               do ir=1,jmt
                  wk(ir)=(vec1(ir)-vec2(ir))**2
               enddo
               vol_int = getVolumeIntegration( id, jmt, r_mesh, 0, wk )
               rho_rms(2,n) = sqrt(vol_int)/mt_volume
               rho_rms_av(2) = rho_rms_av(2) + cfac*rho_rms(2,n)
            endif
!
!           ==========================================================
!           calculate RMS of potential.
!           ==========================================================
            do is = 1, n_spin_pola
               vec1 => getOldSphPotr(id,ia,is)
               vec2 => getNewSphPotr(id,ia,is)
               do ir=1,jmt
                  wk(ir)=( vec1(ir)-vec2(ir) )**2
               enddo
               vol_int = getVolumeIntegration( id, jmt, r_mesh, 0, wk )
               pot_rms(is,n) = sqrt(vol_int)/mt_volume
               pot_rms_av(is) = pot_rms_av(is) + cfac*pot_rms(is,n)
            enddo
         enddo
!
         if (Print_Level(id) >= 0) then
            write(6,'(/,80(''=''))')
            write(6,'(1x,a)')                                         &
                    'AtomId    Iter    RMS of Rho   RMS of Pot      Diff Ef    Diff Etot'
            write(6,'(80(''-''))')
            write(6,'(i7,1x,i7,1x,4(1x,d12.5))')getGlobalIndex(id),   &
                   iteration,max(rho_rms_av(1),rho_rms_av(n_spin_pola)), &
                   max(pot_rms_av(1),pot_rms_av(n_spin_pola)),ef_diff,et_diff
            write(6,'(80(''=''),/)')
         endif
      enddo
!
      deallocate( wk )
   else
      tol = Ten2m14
      allocate( wk_c(nr_max*jmax), wk_prod(nr_max*(2*lmax+1)*(lmax+1)))
      n = 0
      rho_rms = ZERO
      pot_rms = ZERO
      do id = 1, LocalNumAtoms
         Grid => getGrid(id)
         jend = Grid%jend
         r_mesh => Grid%r_mesh(1:jend)
         VP_volume = getVolume(id)
!
         rho_rms_av = ZERO
         pot_rms_av = ZERO
         do ia = 1, getLocalNumSpecies(id)
            n = n + 1
            cfac = getLocalSpeciesContent(id,ia)
!           ==========================================================
!           calculate RMS of density.
!           ==========================================================
            lmax_l = getRhoLmax(id)
            jmax_l = (lmax_l+1)*(lmax_l+2)/2
            kmax_l = (lmax_l+1)*(lmax_l+1)
            p_wk => aliasArray2_c(wk_c,jend,jmax_l)
            flag_jl => getChargeComponentFlag(id)
            lmax_prod = 2*lmax_l
            jmax_prod = (lmax_prod+1)*(lmax_prod+2)/2
            kmax_prod = (lmax_prod+1)**2
            prod => aliasArray2_c(wk_prod,jend,jmax_prod)
            if (n_spin_pola == 1) then
               p_den1 => getChargeDensity("TotalOld",id,ia)
               p_den2 => getChargeDensity("TotalNew",id,ia)
               p_wk = CZERO
               do jl = 1,jmax_l
                  if ( flag_jl(jl) /= 0 ) then
                     do ir=1,jend
                        p_wk(ir,jl) = (p_den1(ir,jl)-p_den2(ir,jl))!*     &
                                   !   conjg(p_den1(ir,jl)-p_den2(ir,jl))
                     enddo
                  endif
               enddo
!              -------------------------------------------------------
               call computeProdExpan(jend,lmax_l,p_wk,lmax_l,p_wk,lmax_prod,prod)
               vol_int = getVolumeIntegration( id, jend, r_mesh,         &
                                               kmax_prod, jmax_prod, 0, prod, vint_mt, tol_in=tol)
!              -------------------------------------------------------
               rho_rms(1,n) = sqrt(abs(vol_int))/VP_volume
               rho_rms_av(1) = rho_rms_av(1) + cfac*rho_rms(1,n)
            else
               p_den1 => getChargeDensity("TotalOld",id,ia)
               p_den2 => getChargeDensity("TotalNew",id,ia)
               p_wk = CZERO
               do jl = 1,jmax_l
                  if ( flag_jl(jl) /= 0 ) then
                     do ir=1,jend
                        p_wk(ir,jl) = (p_den1(ir,jl)-p_den2(ir,jl)) !*     &
                                   !conjg(p_den1(ir,jl)-p_den2(ir,jl))
                     enddo
                  endif
               enddo
!              -------------------------------------------------------
               call computeProdExpan(jend,lmax_l,p_wk,lmax_l,p_wk,lmax_prod,prod)
               vol_int = getVolumeIntegration( id, jend, r_mesh,         &
                                               kmax_prod, jmax_prod, 0, prod, vint_mt, tol_in=tol)
!              -------------------------------------------------------
               rho_rms(1,n) = sqrt(abs(vol_int))/VP_volume
               rho_rms_av(1) = rho_rms_av(1) + cfac*rho_rms(1,n)
               p_den1 => getMomentDensity("TotalOld",id,ia)
               p_den2 => getMomentDensity("TotalNew",id,ia)
               p_wk = CZERO
               do jl = 1,jmax_l
                  if ( flag_jl(jl) /= 0 ) then
                     do ir=1,jend
                        p_wk(ir,jl) = (p_den1(ir,jl)-p_den2(ir,jl)) !*     &
                                     ! conjg(p_den1(ir,jl)-p_den2(ir,jl))
                     enddo
                  endif
               enddo
!              -------------------------------------------------------
               call computeProdExpan(jend,lmax_l,p_wk,lmax_l,p_wk,lmax_prod,prod)
               vol_int = getVolumeIntegration( id, jend, r_mesh,         &
                                               kmax_prod, jmax_prod, 0, prod, vint_mt, tol_in=tol)
!              -------------------------------------------------------
               rho_rms(2,n) = sqrt(abs(vol_int))/VP_volume
               rho_rms_av(2) = rho_rms_av(2) + cfac*rho_rms(2,n)
            endif
!
!           ==========================================================
!           calculate RMS of potential.
!           ==========================================================
            lmax_l = getPotLmax(id)
            jmax_l = (lmax_l+1)*(lmax_l+2)/2
            kmax_l = (lmax_l+1)*(lmax_l+1)
            p_wk => aliasArray2_c(wk_c,jend,jmax_l)
            lmax_prod = 2*lmax_l
            jmax_prod = (lmax_prod+1)*(lmax_prod+2)/2
            kmax_prod = (lmax_prod+1)**2
            prod => aliasArray2_c(wk_prod,jend,jmax_prod)
            do is = 1, n_spin_pola
               flag_jl => getPotComponentFlag(id)
               p_den1 => getOldPotential(id,ia,is)
               p_den2 => getNewPotential("Total",id,ia,is)
               p_wk = CZERO
               do jl = 1,jmax_l
                  if ( flag_jl(jl) /= 0 ) then
                     do ir=1,jend
                        p_wk(ir,jl) = (p_den1(ir,jl)-p_den2(ir,jl)) ! *     &
                                     ! conjg(p_den1(ir,jl)-p_den2(ir,jl))
                     enddo
                  endif
               enddo
!              -------------------------------------------------------
               call computeProdExpan(jend,lmax_l,p_wk,lmax_l,p_wk,lmax_prod,prod)
               vol_int = getVolumeIntegration( id, jend, r_mesh,         &
                                               kmax_prod, jmax_prod, 0, prod, vint_mt, tol_in=tol)
!              -------------------------------------------------------
               pot_rms(is,n) = sqrt(abs(vol_int))/VP_volume
               pot_rms_av(is) = pot_rms_av(is) + cfac*pot_rms(is,n)
            enddo
         enddo
!
         if (Print_Level(id) >= 0) then
            write(6,'(/,80(''=''))')
            write(6,'(1x,a)')                                            &
                    'AtomId    Iter    RMS of Rho   RMS of Pot      Diff Ef    Diff Etot'
            write(6,'(80(''-''))')
            write(6,'(i7,1x,i7,4(1x,d12.5))')getGlobalIndex(id),   &
                   iteration, maxval(rho_rms_av(1:n_spin_pola)),         &
                   maxval(pot_rms_av(1:n_spin_pola)), ef_diff, et_diff
            write(6,'(80(''=''),/)')
         endif
      enddo
      nullify( p_wk, p_den1, p_den2, prod )
      deallocate( wk_c, wk_prod )
   endif
!
   end subroutine checkConvergence
!  ===================================================================
!
end module ConvergenceCheckModule
