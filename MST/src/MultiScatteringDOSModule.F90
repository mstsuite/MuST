module MultiScatteringDOSModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, ONE, TWO, TEN2m6, Y0, PI, PI2, CZERO, SQRTm1, CONE
!
   use PublicTypeDefinitionsModule, only : GridStruct
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use TimerModule, only : getTime
!
   use IntegerFactorsModule, only : lofk, lofj, mofj, m1m, mofk, jofk
!
public :: setMScatteringDOSParam, &
          getMScatteringDOS,        &
          getRelMScatteringDOS
!
   interface setMScatteringDOSParam
      module procedure setMSP_generic, setMSP_ef, setMSP_specific
   end interface setMScatteringDOSParam
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
   logical :: isDensityMatrixNeeded = .false.
   logical :: isZtauZ = .false.
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
   subroutine setMSP_generic(ns,nc,dev,harris,iprint,dm)
!  ===================================================================
   use GroupCommModule, only : getGroupID, getNumPEsInGroup, getMyPEinGroup
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: ns, nc, harris, iprint
!
   logical, intent(in) :: dev
   logical, intent(in), optional :: dm
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
   if (present(dm)) then
      isDensityMatrixNeeded = dm
   else
      isDensityMatrixNeeded = .false.
   endif
!
   end subroutine setMSP_generic
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setMSP_specific(id,nr,jmax_in,ztauz)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, nr, jmax_in
!
   logical, intent(in), optional :: ztauz
!
   local_index = id
   NumRs = nr
   jmax_green = jmax_in
!
   if (present(ztauz)) then
      isZtauZ = ztauz
   else 
      isZtauZ = .false.
   endif
!
   end subroutine setMSP_specific
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setMSP_ef(ef,temp)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: ef, temp
!
   chempot = ef
   Temperature = temp
!
   end subroutine setMSP_ef
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMScatteringDOS(info,e,dos_array,cfac) result(dos)
!  ===================================================================
   use MathParamModule, only : CONE, ZERO, HALF, PI, Y0
!
   use MSSolverModule, only : getMSGreenFunction, getMSGreenMatrix
   use MSSolverModule, only : getMSGreenFunctionDerivative
!
   use SSSolverModule, only : getFreeElectronDOS, getTMatrix
!
   use PotentialTypeModule, only : isASAPotential
!
   use RadialGridModule, only : getGrid, getRadialIntegration
!
   use AtomModule, only : getLocalEvec, getLocalNumSpecies
!
   use StepFunctionModule, only : getVolumeIntegration
!
   use SpinRotationModule, only : transformDensityMatrix
!
   use WriteMatrixModule,  only : writeMatrix
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), intent(in), optional :: cfac
!
   complex (kind=CmplxKind), intent(inout), target :: dos_array(:)
!
   integer (kind=IntKind), intent(in) :: info(:)
   integer (kind=IntKind) :: is, id, atom, ia, ir, nr
   integer (kind=IntKind) :: ks, jid, j, l, m, jl, kl, klc, p0, n, i
   integer (kind=IntKind) :: jmax, kmax, iend, ns_sqr, gform
!
   real (kind=RealKind) :: r0j(3), evec_new(3), sfac, dos(2)
   real (kind=RealKind), pointer :: r_mesh(:), pra_y(:,:)
!
   complex (kind=CmplxKind) :: dosmt_pola(2)
   complex (kind=CmplxKind) :: dosws_pola(2), ede, energy, cmul
!
   complex (kind=CmplxKind) :: greenint(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: greenint_mt(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: gm(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: gm_mt(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind), pointer :: green_matrix(:,:,:,:)
   complex (kind=CmplxKind), pointer :: green(:,:,:,:)
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:), pca_x(:,:), pca_y(:,:)
   complex (kind=CmplxKind), pointer :: der_green(:,:,:,:), der_dos_r_jl(:,:)
   complex (kind=CmplxKind), pointer :: pcv_x(:), pcv_y(:), pcv_z(:), p1(:)
!
   real (kind=RealKind), pointer :: fedos(:)
   real (kind=RealKind), allocatable :: der_fedos(:)
!
   type (GridStruct), pointer :: Grid
!
   logical :: isPositive
!
   is = info(1); id = info(2); atom = info(3)
!
   if (id /= local_index) then
      call ErrorHandler('getMScatteringDOS','Inconsistent local index',id,local_index)
   endif
!
   energy = adjustEnergy(is,e)
!
   if (present(cfac)) then
      cmul = cfac
   else
      cmul = CONE
   endif
!
   if (rad_derivative) then
      allocate(der_fedos(NumRs))
   endif
!
!  =================================================================
!  Note: green is the multiple scattering term of the Green function
!        multiplied by r^2 and is expanded in spherical harmonics.
!
!        For nonspin-polarized case, green, greenint, and
!        greenint_mt also contain a factor 2.0.
!
!        green and greenint are in local spin reference frame
!  -----------------------------------------------------------------
   green => getMSGreenFunction(id,gform)
!  -----------------------------------------------------------------
!  write(6,'(a,2d15.8)')'green(900) = ',green(900,1,1,1)
   if (rad_derivative) then
!     --------------------------------------------------------------
      der_green => getMSGreenFunctionDerivative(id)
!     --------------------------------------------------------------
   endif
!
   sfac= TWO/real(n_spin_pola,kind=RealKind)
   ns_sqr = n_spin_cant*n_spin_cant
   iend = size(green,1); kmax = size(green,2); jmax = jofk(kmax)
!
   if (iend /= NumRs) then
      call ErrorHandler('getMScatteringDOS','Inconsistent r-mesh size of green and dos_r_jl', &
                        iend,NumRs)
   else if (jmax /= jmax_green) then
      call ErrorHandler('getMScatteringDOS','Inconsistent jmax size of green and dos_r_jl', &
                        jmax,jmax_green)
   endif
!
   Grid => getGrid(id)
!
   p0 = 0
   do ia = 1, getLocalNumSpecies(id)
      if (atom < 0 .or. ia == atom) then
         do ks = 1, ns_sqr
            if (isASAPotential()) then
               greenint_mt(ks) = cmul*sfac*getRadialIntegration(id, Grid%jmt, green(:,1,ks,ia))/Y0
               greenint(ks) = greenint_mt(ks)
!              Alternatively, you can call getVolumeIntegration with truncated=.false.
!              to run the integration in ASA volume.
            else
!              -----------------------------------------------------
               greenint(ks) = cmul*sfac*getVolumeIntegration( id, iend, Grid%r_mesh, &
                                                              kmax, 2, green(:,:,ks,ia), greenint_mt(ks) )
!              -----------------------------------------------------
               greenint_mt(ks) = cmul*sfac*greenint_mt(ks)
            endif
         enddo
!
         if (isASAPotential()) then
            gm_mt = greenint_mt
            if (n_spin_cant == 2) then
!              ------------------------------------------------------
               call transformDensityMatrix(id,gm_mt)
!              ------------------------------------------------------
            endif
            gm = gm_mt
         else
            gm = greenint
            gm_mt = greenint_mt
            if (n_spin_cant == 2) then
!              ------------------------------------------------------
               call transformDensityMatrix(id,gm)
               call transformDensityMatrix(id,gm_mt)
!              ------------------------------------------------------
            endif
         endif
!
!        ===========================================================
!        Note: dosmt, dosws, and dos_array(p0+1:p0+3) are complex type, and
!              their real part are MT-DOS, WS-DOS, and e*DOS, respectively.
!        ===========================================================
         dosmt_pola = CZERO; dosws_pola = CZERO
         if(n_spin_cant == 2) then
!           ---------------------------------------------------------
            evec_new(1:3) = getLocalEvec(id,'new')
            dosmt_pola(1) = SQRTm1*HALF*(gm_mt(1)+gm_mt(2)*evec_new(1)+gm_mt(3)*evec_new(2)+gm_mt(4)*evec_new(3))/PI
            dosmt_pola(2) = SQRTm1*HALF*(gm_mt(1)-gm_mt(2)*evec_new(1)-gm_mt(3)*evec_new(2)-gm_mt(4)*evec_new(3))/PI
            dosws_pola(1) = SQRTm1*HALF*(gm(1)+gm(2)*evec_new(1)+gm(3)*evec_new(2)+gm(4)*evec_new(3))/PI
            dosws_pola(2) = SQRTm1*HALF*(gm(1)-gm(2)*evec_new(1)-gm(3)*evec_new(2)-gm(4)*evec_new(3))/PI
         else
            dosmt_pola(1) = SQRTm1*greenint_mt(1)/PI
            dosws_pola(1) = SQRTm1*greenint(1)/PI
         endif
!
         if (node_print_level >= 0) then
            write(6,'(/,''getMScatteringDOS: energy ='',2f18.12,'', id ='',i4,'', ia ='',i4)')energy,id,ia
            if (gform == 0) then
               write(6,'(  ''                       Int [Z*(Tau-t)*Z] on MT         ='',2f18.12)') &
                     dosmt_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Int [Z*(Tau-t)*Z] on VP         ='',2f18.12)') &
                     dosws_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Im{-Int [Z*(Tau-t)*Z]/pi on MT} ='',f18.12)')  &
                     real(dosmt_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
               write(6,'(  ''                       Im{-Int [Z*(Tau-t)*Z]/pi on VP} ='',f18.12)')  &
                     real(dosws_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
            else if (gform == 1) then
               write(6,'(  ''                       Int [Z*Tau*Z-Z*J] on MT         ='',2f18.12)') &
                     dosmt_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Int [Z*Tau*Z-Z*J] on VP         ='',2f18.12)') &
                     dosws_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Im{-Int [Z*Tau*Z-Z*J]/pi on MT} ='',f18.12)')  &
                     real(dosmt_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
               write(6,'(  ''                       Im{-Int [Z*Tau*Z-Z*J]/pi on VP} ='',f18.12)')  &
                     real(dosws_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
            else
               write(6,'(  ''                       Int [Z*Tau*Z] on MT         ='',2f18.12)')     &
                     dosmt_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Int [Z*Tau*Z] on VP         ='',2f18.12)')     &
                     dosws_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
               write(6,'(  ''                       Im{-Int [Z*Tau*Z]/pi on MT} ='',f18.12)')      &
                     real(dosmt_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
               write(6,'(  ''                       Im{-Int [Z*Tau*Z]/pi on VP} ='',f18.12)')      &
                     real(dosws_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
            endif
         endif ! print_level
!
         if (atom > 0 .or. ia == 1) then
            dos = real(dosws_pola,kind=RealKind)
         endif
!
         if (isDensityMatrixNeeded) then
!           ---------------------------------------------------------
            green_matrix => getMSGreenMatrix(id)
!           ---------------------------------------------------------
         endif
!
         if (iharris <= 1) then
            ede = energy
         else
            ede = (energy-chempot)
         endif
!
         n = iend*jmax
         do ks = 1, ns_sqr
            p1 => dos_array(p0+1:p0+n)
            dos_r_jl => aliasArray2_c(p1, iend, jmax)
            do jl = 1, jmax
               l = lofj(jl); m = mofj(jl)
               kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
               pcv_x => green(:,kl,ks,ia)
               pcv_y => green(:,klc,ks,ia)
               pcv_z => dos_r_jl(:,jl)
               pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
            enddo
            p0 = p0 + n
            if (rad_derivative) then
               p1 => dos_array(p0+1:p0+n)
               der_dos_r_jl => aliasArray2_c(p1, iend, jmax)
               do jl = 1, jmax
                  l = lofj(jl); m = mofj(jl)
                  kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
                  pcv_x => der_green(:,kl,ks,ia)
                  pcv_y => der_green(:,klc,ks,ia)
                  pcv_z => der_dos_r_jl(:,jl)
                  pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
               enddo
               p0 = p0 + n
            endif
            dos_array(p0+1)  = SQRTm1*greenint(ks)/PI
            dos_array(p0+2) = SQRTm1*greenint_mt(ks)/PI
            dos_array(p0+3) = SQRTm1*ede*greenint(ks)/PI
            dos_array(p0+4) = CZERO
            p0 = p0 + 4
            if (isDensityMatrixNeeded) then
               kmax = size(green_matrix,1)
               pca_x => green_matrix(:,:,ks,ia)
               p1 => dos_array(p0+1:p0+kmax*kmax)
               pca_y => aliasArray2_c(p1,kmax,kmax)
               pca_y = cmul*SQRTm1*pca_x
               p0 = p0 + kmax*kmax
            endif
! =============================================================================
! Note: It is possible that due to numerical round off error, the imaginary part
!       of Green function becomes negative. If this happens, in non-full-potential 
!       case, DOS is replaced with the free-electron DOS. An investigation of this
!       issue is necessary in the future.  06/14/2020....
! =============================================================================
            if (jmax == 1 .and. isZtauZ) then
               isPositive = .false.
               LOOP_ir: do ir = Grid%jend, 1, -1
                  if (aimag(green(ir,1,ks,ia)) > ZERO) then
!                    write(6,'(a,2i5,2d16.8)')'id,ir,r,-Im[green(ir)]/PI = ',     &
!                          id,ir,Grid%r_mesh(ir),-aimag(green(ir,1,ks,ia))/Grid%r_mesh(ir)**2/PI
                     nr = ir
                     isPositive = .true.
                     exit LOOP_ir
                  endif
               enddo LOOP_ir
               if (isPositive) then
                  if (rad_derivative) then
                     fedos => getFreeElectronDOS(site=id,gfac=cmul,derivative=der_fedos)
                  else
                     fedos => getFreeElectronDOS(site=id,gfac=cmul)
                  endif
                  if (rad_derivative) then
                     do ir = 1, nr
                        if (aimag(green(ir,1,ks,ia)) > ZERO) then
                           dos_r_jl(ir,1) = sfac*fedos(ir)
                           der_dos_r_jl(ir,1) = sfac*der_fedos(ir)
                        endif
                     enddo
                  else
                     do ir = 1, nr
                        if (aimag(green(ir,1,ks,ia)) > ZERO) then
                           dos_r_jl(ir,1) = sfac*fedos(ir)
                        endif
                     enddo
                  endif
               endif
            endif
! =============================================================================
         enddo
      endif
   enddo
   if (rad_derivative) then
      deallocate(der_fedos)
   endif
   nullify(fedos)
!
   end function getMScatteringDOS
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRelMScatteringDOS(info,e,dos_array,cfac) result(dos)
!  ===================================================================
   use MathParamModule, only : CONE, ZERO, HALF, PI
!
   use RelMSSolverModule, only : getRelMSGreenFunction!, getMSGreenMatrix
!
   use PotentialTypeModule, only : isASAPotential
!
   use RadialGridModule, only : getGrid, getRadialIntegration
!
   use AtomModule, only : getLocalEvec
!
   use StepFunctionModule, only : getVolumeIntegration
!
!   use SpinRotationModule, only : transformDensityMatrix
!
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: e
   complex (kind=CmplxKind), intent(in), optional :: cfac
!
   complex (kind=CmplxKind), intent(inout), target :: dos_array(*)
!
   integer (kind=IntKind), intent(in) :: info(*)
   integer (kind=IntKind) :: is, id, atom
   integer (kind=IntKind) :: ks, jid, j, l, m, jl, kl, klc, p0, n, i
   integer (kind=IntKind) :: jmax, kmax, iend, ns_sqr, gform
!
   real (kind=RealKind) :: r0j(3), evec_new(3), sfac, dos(2)
   real (kind=RealKind), pointer :: r_mesh(:), pra_y(:,:)
!
   complex (kind=CmplxKind) :: dosmt_pola(2)
   complex (kind=CmplxKind) :: dosws_pola(2), ede, energy, cmul
!
   complex (kind=CmplxKind) :: greenint(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: greenint_mt(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: gm(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind) :: gm_mt(n_spin_cant*n_spin_cant)
   complex (kind=CmplxKind), pointer :: green_matrix(:,:,:)
   complex (kind=CmplxKind), pointer :: green(:,:,:)
   complex (kind=CmplxKind), pointer :: dos_r_jl(:,:), pca_x(:,:), pca_y(:,:)
   complex (kind=CmplxKind), pointer :: pcv_x(:), pcv_y(:), pcv_z(:), p1(:)
!
   type (GridStruct), pointer :: Grid
!
   is = info(1); id = info(2); atom = info(3)
!
   if (id /= local_index) then
      call ErrorHandler('getRelMScatteringDOS','Inconsistent local index',id,local_index)
   endif
!
   energy = adjustEnergy(is,e)
!
   if (present(cfac)) then
      cmul = cfac
   else
      cmul = CONE
   endif
!
!  =================================================================
!  Note: green is the multiple scattering term of the Green function
!        multiplied by r^2 and is expanded in spherical harmonics.
!
!        For nonspin-polarized case, green, greenint, and
!        greenint_mt also contain a factor 2.0.
!
!        green and greenint are in local spin reference frame
!  -----------------------------------------------------------------
   green => getRelMSGreenFunction(id)
!  -----------------------------------------------------------------
!
   sfac= ONE !TWO/real(n_spin_pola,kind=RealKind)
   ns_sqr = n_spin_cant*n_spin_cant
   iend = size(green,1); kmax = size(green,2); jmax = jofk(kmax)
!
   if (iend /= NumRs) then
      call ErrorHandler('getRelMScatteringDOS','Inconsistent r-mesh size of green and dos_r_jl', &
                        iend,NumRs)
   else if (jmax /= jmax_green) then
      call ErrorHandler('getRelMScatteringDOS','Inconsistent jmax size of green and dos_r_jl', &
                        jmax,jmax_green)
   endif
!
   Grid => getGrid(id)
!
   do ks = 1, ns_sqr
      if (isASAPotential()) then
         greenint_mt(ks) = cmul*sfac*getRadialIntegration(id, Grid%jmt, green(:,1,ks))/Y0
         greenint(ks) = greenint_mt(ks)
!        Alternatively, you can call getVolumeIntegration with truncated=.false.
!        to run the integration in ASA volume.
      else
!        -----------------------------------------------------------
         greenint(ks) = cmul*sfac*getVolumeIntegration( id, iend, Grid%r_mesh, &
                                                        kmax, 2, green(:,:,ks), greenint_mt(ks) )
!        -----------------------------------------------------------
         greenint_mt(ks) = cmul*sfac*greenint_mt(ks)
      endif
   enddo
!
!   if (isASAPotential()) then
!      gm_mt = greenint_mt
!      if (n_spin_cant == 2) then
!        ------------------------------------------------------------
!         call transformDensityMatrix(id,gm_mt)
!        ------------------------------------------------------------
!      endif
!      greenint = greenint_mt
!      gm = gm_mt
!   else
      gm = greenint
      gm_mt = greenint_mt
!      if (n_spin_cant == 2) then
!        ------------------------------------------------------------
!         call transformDensityMatrix(id,gm)
!         call transformDensityMatrix(id,gm_mt)
!        ------------------------------------------------------------
!      endif
!   endif
!
   dosmt_pola = CZERO; dosws_pola = CZERO
   if(n_spin_cant == 2) then
!     ---------------------------------------------------------------
      evec_new(1:3) = getLocalEvec(id,'new')
      dosmt_pola(1) = SQRTm1*HALF*(gm_mt(1)+gm_mt(2)*evec_new(1)+gm_mt(3)*evec_new(2)+gm_mt(4)*evec_new(3))/PI
      dosmt_pola(2) = SQRTm1*HALF*(gm_mt(1)-gm_mt(2)*evec_new(1)-gm_mt(3)*evec_new(2)-gm_mt(4)*evec_new(3))/PI
      dosws_pola(1) = SQRTm1*HALF*(gm(1)+gm(2)*evec_new(1)+gm(3)*evec_new(2)+gm(4)*evec_new(3))/PI
      dosws_pola(2) = SQRTm1*HALF*(gm(1)-gm(2)*evec_new(1)-gm(3)*evec_new(2)-gm(4)*evec_new(3))/PI
   else
      dosmt_pola(1) = SQRTm1*greenint_mt(1)/PI
      dosws_pola(1) = SQRTm1*greenint(1)/PI
   endif

   if (node_print_level >= 0) then
      write(6,'(/,''getRelMScatteringDOS: energy ='',2f18.12,'', id ='',i4)')energy,id
      write(6,'(  ''                       Int [Z*(Tau-t)*Z] on MT         ='',2f18.12)') &
            dosmt_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
      write(6,'(  ''                       Int [Z*(Tau-t)*Z] on VP         ='',2f18.12)') &
            dosws_pola(1:n_spin_cant)*PI/(SQRTm1*cmul*sfac)
      write(6,'(  ''                       Im{-Int [Z*(Tau-t)*Z]/pi on MT} ='',f18.12)')  &
            real(dosmt_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
      write(6,'(  ''                       Im{-Int [Z*(Tau-t)*Z]/pi on VP} ='',f18.12)')  &
            real(dosws_pola(1:n_spin_cant)/(cmul*sfac),RealKind)
   endif ! print_level

   dos = real(dosws_pola,kind=RealKind)

   if (iharris /= 0) then
      call ErrorHandler('getRelMScatteringDOS','Relativistic Harris calculation not implemented')
   endif
!
   if (iharris <= 1) then !harris energy calculation,
      ede = energy
   else
      ede = (energy-chempot)
   endif
!
   n = iend*jmax
   p0 = 0
   do ks = 1, ns_sqr
      p1 => dos_array(p0+1:p0+n)
      dos_r_jl => aliasArray2_c(p1, iend, jmax)
      do jl = 1, jmax
         l = lofj(jl); m = mofj(jl)
         kl = (l+1)**2-l+m; klc = (l+1)**2-l-m
         pcv_x => green(:,kl,ks)
         pcv_y => green(:,klc,ks)
         pcv_z => dos_r_jl(:,jl)
         pcv_z = sfac*SQRTm1*(cmul*pcv_x-m1m(m)*conjg(cmul*pcv_y))/PI2
      enddo
      p0 = p0 + n
      dos_array(p0+1)  = SQRTm1*greenint(ks)/PI
      dos_array(p0+2) = SQRTm1*greenint_mt(ks)/PI
      dos_array(p0+3) = SQRTm1*ede*greenint(ks)/PI
      dos_array(p0+4) = CZERO
      p0 = p0 + 4
!      if (isDensityMatrixNeeded) then
!         kmax = size(green_matrix,1)
!         pca_x => green_matrix(:,:,ks)
!         p1 => dos_array(p0+1:p0+kmax*kmax)
!         pca_y => aliasArray2_c(p1,kmax,kmax)
!         pca_y = cmul*SQRTm1*pca_x
!         p0 = p0 + kmax*kmax
!      endif
   enddo
   end function getRelMScatteringDOS
!  ===================================================================
end module MultiScatteringDOSModule
