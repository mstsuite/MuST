module SCPolesModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   implicit none
!
public :: initSCPoles,      &
          endSCPoles,       &
          computeSCPoles,   &
          getNumSCPoles,    &
          getSCPole,        &
          getSCPoleResidule   
!
private
   integer (kind=IntKind) :: LocalNumAtoms, n_spin_cant
   integer (kind=IntKind), allocatable :: lmax_kkr(:)
   integer (kind=IntKind), allocatable :: lmax_phi(:)
!
!
   type PoleStruct
      integer (kind=IntKind) :: NumPoles
      real (kind=RealKind), pointer :: Poles(:)
      complex (kind=CmplxKind), pointer :: Residules(:,:)
   end type PoleStruct
!
   type (PoleStruct), allocatable :: SCPoles(:,:)
   real (kind=RealKind), allocatable :: Poles(:)
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initSCPoles(lna, lkkr, lphi, nsc)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: lna, nsc
   integer (kind=IntKind), intent(in) :: lkkr(lna), lphi(lna)
   integer (kind=IntKind) :: id, lmax_kkr_max
!
   LocalNumAtoms = lna
   n_spin_cant = nsc
!
!  ===================================================================
!  setup the Lmax values.             
!  ===================================================================
   allocate(lmax_kkr(LocalNumAtoms), lmax_phi(LocalNumAtoms))
   lmax_kkr_max = 0
   do id=1,LocalNumAtoms
      lmax_kkr(id) = lkkr(id)
      lmax_phi(id) = lphi(id)
      lmax_kkr_max = max(lmax_kkr_max, lmax_kkr(id))
   enddo
!
   allocate( SCPoles(n_spin_cant,LocalNumAtoms) )
   allocate( Poles((lmax_kkr_max+1)**2) )
!
   end subroutine initSCPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endSCPoles()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: is, id
!
   deallocate(lmax_kkr, lmax_phi)
!
   do id = 1, LocalNumAtoms
      do is = 1, n_spin_cant
         if (SCPoles(is,id)%NumPoles > 0) then
            deallocate( SCPoles(is,id)%Poles )
            if ( associated(SCPoles(is,id)%Residules) ) then
               deallocate( SCPoles(is,id)%Residules )
            endif
            nullify( SCPoles(is,id)%Poles, SCPoles(is,id)%Residules )
         endif
      enddo
   enddo
   deallocate( SCPoles, Poles )
!
   end subroutine endSCPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeSCPoles(eb,et)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: is, id, ip
!
   real (kind=RealKind), intent(in) :: eb, et
!
!  -------------------------------------------------------------------
   call calQuadraticPoles(eb,et)
!  -------------------------------------------------------------------
!
   do id = 1,LocalNumAtoms
      do is = 1,n_spin_cant
         write(6,'(a,i3,a,i3)')'Spin index: ',is,',  Atom index: ',id
         write(6,'(a,f10.5,a,f10.5,a,i5)')'The Number of poles found within (',Eb,', ', &
                                           Et,'): ',SCPoles(is,id)%NumPoles
         do ip = 1, SCPoles(is,id)%NumPoles
            write(6,'(f20.12)')SCPoles(is,id)%Poles(ip)
         enddo
      enddo
   enddo
!
   end subroutine computeSCPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumSCPoles(is,id) result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id, is
   integer (kind=IntKind) :: n
!
   n = SCPoles(is,id)%NumPoles
!
   end function getNumSCPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSCPole(is,id,ip) result(e)
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, ip
   real (kind=RealKind) :: e
!
   if (ip < 1 .or. ip > SCPoles(is,id)%NumPoles) then
      call ErrorHandler('getSCPole','The pole index is out of range',ip)
   endif
!
   e = SCPoles(is,id)%Poles(ip)
!
   end function getSCPole
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getSCPoleResidule(is,id,ip) result(res)
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
   implicit none
!
   integer (kind=IntKind), intent(in) :: id, is, ip
   integer (kind=IntKind) :: nsize
!
   complex (kind=CmplxKind), pointer :: res(:)
!
   if (ip < 1 .or. ip > SCPoles(is,id)%NumPoles) then
      call ErrorHandler('getSCPole','The pole index is out of range',ip)
   endif
!
   nsize = (lmax_kkr(id)+1)**4
!
   res => SCPoles(is,id)%Residules(1:nsize,ip)
!
   end function getSCPoleResidule   
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calQuadraticPoles(eb,et)
!  ===================================================================
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CZERO, CONE, Ten2m8, SQRTm1
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use SSSolverModule, only : getSineMatrix, getCosineMatrix, solveSingleScattering, &
                              getJostMatrix
!
   use QuadraticMatrixModule, only : initQuadraticMatrix, endQuadraticMatrix, &
                                     solveQuadraticEquation, getEigenValue
!
   use MatrixDeterminantModule, only : MtxDet
!
   implicit none
!
   integer (kind=IntKind) :: id, is, ie, iw, j, n, kmax_kkr, kmax_phi, info
   integer (kind=IntKind) :: lmax_kkr_max, lmax_phi_max
   integer (kind=IntKind) :: kmax_kkr_max, kmax_phi_max, kl, klp, m, mp, nv
!
   real (kind=RealKind), intent(in) :: eb, et
   real (kind=RealKind) :: WindowWidth = 0.005d0  ! Ryd unit
   real (kind=RealKind) :: e0, de, de2, dede2, pe, w0, w1
!
   complex (kind=CmplxKind) :: e, det
   complex (kind=CmplxKind), pointer :: sin_mat(:,:), jost_mat(:,:)
   complex (kind=CmplxKind), pointer :: pv(:)
   complex (kind=CmplxKind), pointer :: s0(:,:), s1(:,:), s2(:,:), st(:,:)
   complex (kind=CmplxKind), pointer :: jost0(:,:), jost1(:,:), jost2(:,:)
   complex (kind=CmplxKind), pointer :: sj0(:,:), sj1(:,:), sj2(:,:)
   complex (kind=CmplxKind), pointer :: am(:,:)
   complex (kind=CmplxKind), pointer :: sin_t(:,:), p1c(:)
   complex (kind=CmplxKind), allocatable, target :: wks_sc(:,:)
!
   if (abs(eb-et) < TEN2m6) then
      call ErrorHandler('calQuadraticPoles','et - eb < 0.000001',eb,et)
   else
      WindowWidth = min(abs(et-eb),WindowWidth)
   endif
!
   lmax_kkr_max = 0
   lmax_phi_max = 0
   do id = 1,LocalNumAtoms
      lmax_kkr_max = max(lmax_kkr_max,lmax_kkr(id))
      lmax_phi_max = max(lmax_phi_max,lmax_phi(id))
   enddo 
   kmax_kkr_max = (lmax_kkr(1)+1)**2
   kmax_phi_max = (lmax_phi(1)+1)**2
   allocate( wks_sc(kmax_phi_max*kmax_kkr_max,14) )
!
   kmax_kkr = kmax_kkr_max
!  -------------------------------------------------------------------
   call initQuadraticMatrix(kmax_kkr)
!  -------------------------------------------------------------------
   do id = 1,LocalNumAtoms
      if ((lmax_kkr(id)+1)**2 /= kmax_kkr) then
!        -------------------------------------------------------------
         call endQuadraticMatrix()
         call initQuadraticMatrix(kmax_kkr)
!        -------------------------------------------------------------
      endif
!
      kmax_kkr = (lmax_kkr(id)+1)**2
      kmax_phi = (lmax_phi(id)+1)**2
      p1c => wks_sc(1:kmax_phi*kmax_kkr,1)
      s0 => aliasArray2_c(p1c,kmax_phi,kmax_kkr)
      p1c => wks_sc(1:kmax_phi*kmax_kkr,2)
      s1 => aliasArray2_c(p1c,kmax_phi,kmax_kkr)
      p1c => wks_sc(1:kmax_phi*kmax_kkr,3)
      s2 => aliasArray2_c(p1c,kmax_phi,kmax_kkr)
      p1c => wks_sc(1:kmax_phi*kmax_kkr,4)
      st => aliasArray2_c(p1c,kmax_phi,kmax_kkr)
      p1c => wks_sc(1:kmax_phi*kmax_kkr,5)
      jost0 => aliasArray2_c(p1c,kmax_phi,kmax_kkr)
      p1c => wks_sc(1:kmax_phi*kmax_kkr,6)
      jost1 => aliasArray2_c(p1c,kmax_phi,kmax_kkr)
      p1c => wks_sc(1:kmax_phi*kmax_kkr,7)
      jost2 => aliasArray2_c(p1c,kmax_phi,kmax_kkr)
      p1c => wks_sc(1:kmax_kkr*kmax_kkr,8)
      sj0 => aliasArray2_c(p1c,kmax_kkr,kmax_kkr)
      p1c => wks_sc(1:kmax_kkr*kmax_kkr,9)
      sj1 => aliasArray2_c(p1c,kmax_kkr,kmax_kkr)
      p1c => wks_sc(1:kmax_kkr*kmax_kkr,10)
      sj2 => aliasArray2_c(p1c,kmax_kkr,kmax_kkr)
      p1c => wks_sc(1:kmax_kkr*kmax_kkr,11)
      am => aliasArray2_c(p1c,kmax_kkr,kmax_kkr)
      p1c => wks_sc(1:kmax_phi*kmax_kkr,12)
      sin_t => aliasArray2_c(p1c,kmax_phi,kmax_kkr)
!
      do is = 1,n_spin_cant
         n = 0
         w1 = min(eb,et)
         DO_WHILE: do while (w1 < max(eb,et))
            if (w1 < ZERO .and. w1 + WindowWidth > -0.001) then  ! Need to avoid the panel including e = 0
               if (w1 < -0.01d0) then     ! In this case, keep the panel on the left of e = 0
                  w0 = w1; w1 = -0.005d0
               else                       ! In this case, shift the panel to the right of e = 0
                  w0 = 0.005d0
                  w1 = w0 + WindowWidth
               endif
            else                          ! In this case, panel will not cover e = 0
               w0 = w1
               w1 = w0 + WindowWidth
            endif
            if (w1 > et) then
               w1 = et
            endif
            if (w1 < w0) then
               exit DO_WHILE
            else if (w0 <= 0.005) then  ! ignore the panels too close to the origin or in the negative energy region.
               cycle DO_WHILE
            endif
!
            e0 = (w0+w1)*HALF
            de = (w1-w0)/3.0d0; de2 = de*TWO; dede2 = de*de*TWO
!
            e = cmplx(e0,ZERO,kind=CmplxKind)
!           ----------------------------------------------------------
            call solveSingleScattering(is,id,e,CZERO)
!           ----------------------------------------------------------
            sin_mat => getSineMatrix()
            call constructSinT(kmax_phi,kmax_kkr,sin_mat,s0)
            jost_mat => getJostMatrix()
            jost0 = jost_mat
!           ----------------------------------------------------------
!
            e = cmplx(e0+de,ZERO,kind=CmplxKind)
!           ----------------------------------------------------------
            call solveSingleScattering(is,id,e,CZERO)
!           ----------------------------------------------------------
            sin_mat => getSineMatrix()
            call constructSinT(kmax_phi,kmax_kkr,sin_mat,s2)
            jost_mat => getJostMatrix()
            jost2 = jost_mat
!           ----------------------------------------------------------
!
            e = cmplx(e0-de,ZERO,kind=CmplxKind)
!           ----------------------------------------------------------
            call solveSingleScattering(is,id,e,CZERO)
!           ----------------------------------------------------------
            sin_mat => getSineMatrix()
            call constructSinT(kmax_phi,kmax_kkr,sin_mat,st)
            jost_mat => getJostMatrix()
!           ----------------------------------------------------------
!
            s1 = (s2 - st)/de2
            s2 = (s2 + st - TWO*s0)/dede2
!
            jost1 = (jost2 - jost_mat)/de2
            jost2 = (jost2 + jost_mat - TWO*jost0)/dede2
!
!           ----------------------------------------------------------
            call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_phi, CONE,          &
                        s0, kmax_phi, jost0, kmax_phi, CZERO, sj0, kmax_kkr)
!           ----------------------------------------------------------
!
!           ----------------------------------------------------------
            call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_phi, CONE,          &
                        s0, kmax_phi, jost1, kmax_phi, CZERO, sj1, kmax_kkr)
!           ----------------------------------------------------------
            call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_phi, CONE,          &
                        s1, kmax_phi, jost0, kmax_phi, CONE, sj1, kmax_kkr)
!           ----------------------------------------------------------
!
!           ----------------------------------------------------------
            call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_phi, CONE,          &
                        s0, kmax_phi, jost2, kmax_phi, CZERO, sj2, kmax_kkr)
!           ----------------------------------------------------------
            call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_phi, CONE,          &
                        s1, kmax_phi, jost1, kmax_phi, CONE, sj2, kmax_kkr)
!           ----------------------------------------------------------
            call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_phi, CONE,          &
                        s2, kmax_phi, jost0, kmax_phi, CONE, sj2, kmax_kkr)
!           ----------------------------------------------------------
!
!           ----------------------------------------------------------
            call solveQuadraticEquation(sj0,sj1,sj2,info)
!           ----------------------------------------------------------
!
            if (info /= 0) then
!              -------------------------------------------------------
               call ErrorHandler('calQuadraticPoles','Encountered ill conditioned sc0, sc1, or sc2 matrix')
!              -------------------------------------------------------
            endif
!
!           ----------------------------------------------------------
            pv => getEigenValue(nv)
!           ----------------------------------------------------------
            LOOP_ie: do ie = 1, kmax_kkr*2
               if (abs(aimag(pv(ie))) < Ten2m8) then
                  pe = real(pv(ie),kind=RealKind) + e0
                  if (pe >= w0 .and. pe <= w1) then
                     do j = 1, n
                        if (abs(Poles(j) - pe) < TEN2m8) then
                           cycle LOOP_ie
                        endif
                     enddo
                     n = n + 1
                     Poles(n) = pe
                  endif
               endif
            enddo LOOP_ie
         enddo DO_WHILE
!
         SCPoles(is,id)%NumPoles = n
         if (n > 0) then
            allocate( SCPoles(is,id)%Poles(n) )
         endif
         do n = 1, SCPoles(is,id)%NumPoles
            SCPoles(is,id)%Poles(n) = Poles(n)
            if (.false.) then
               do j = -1, 1
                  e = cmplx(Poles(n)+j*de,ZERO,kind=CmplxKind)
!                 ----------------------------------------------------
                  call solveSingleScattering(is,id,e,CZERO)
!                 ----------------------------------------------------
                  sin_mat => getSineMatrix()
                  jost_mat => getJostMatrix()
                  call constructSinT(kmax_phi,kmax_kkr,sin_mat,sin_t)
!                 ----------------------------------------------------
                  call zgemm( 't', 'n', kmax_kkr, kmax_kkr, kmax_phi, CONE,          &
                              sin_t, kmax_phi, jost_mat, kmax_phi, CZERO, am, kmax_kkr)
!                 ----------------------------------------------------
                  call MtxDet(kmax_kkr,am,det)
                  write(6,'(a,f12.6,2x,2d15.8)')'energy, determ = ',real(e,kind=RealKind),det
!                 ----------------------------------------------------
               enddo
            endif
         enddo
!
      enddo
   enddo
!
!  -------------------------------------------------------------------
   call endQuadraticMatrix()
!  -------------------------------------------------------------------
!
   nullify( s0, s1, s2, st, jost0, jost1, jost2, sj0, sj1, sj2, sin_t, am, sin_mat, jost_mat, pv )
   deallocate( wks_sc )
!
   end subroutine calQuadraticPoles
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine constructSinT(kmax_phi,kmax_kkr,sin_mat,sin_t)
!  ===================================================================
   use IntegerFactorsModule, only : mofk, m1m
   implicit none
!
   integer (kind=IntKind), intent(in) :: kmax_phi, kmax_kkr
!
   complex (kind=CmplxKind), intent(in) :: sin_mat(kmax_phi,kmax_kkr)
   complex (kind=CmplxKind), intent(out) :: sin_t(kmax_phi,kmax_kkr)
!
   integer (kind=IntKind) :: kl, m, klp, mp, klpp
!
   do kl = 1,kmax_kkr
      m = mofk(kl)
      do klp = 1,kmax_phi
         mp = mofk(klp)
         sin_t(klp,kl) = m1m(m+mp)*sin_mat(klp-2*mp,kl-2*m)
      enddo
   enddo
!
   end subroutine constructSinT
!  ===================================================================
end module SCPolesModule
