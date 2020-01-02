!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calScatteringPoles(LocalNumAtoms,lkkr,lphi,             &
                                 n_spin_pola,n_spin_cant,e0,e1,ei,nstep)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : ZERO, TEN2m6, HALF, ONE, TWO, CONE, Ten2m8, CZERO
!
   use MatrixInverseModule, only : MtxInv_GE
   use LUdcmpModule, only : LUdcmp
   use SSSolverModule, only : getSineMatrix,                    &
                                     solveSingleScattering
!
   implicit none
!
   character(len=9)  :: char_lm
!
   logical :: isLUdcmp = .false.
!
   integer (kind=IntKind), parameter :: offset_lm = 100
!
   integer (kind=IntKind), intent(in) :: LocalNumAtoms
   integer (kind=IntKind), intent(in) :: lkkr(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: lphi(1:LocalNumAtoms)
   integer (kind=IntKind), intent(in) :: n_spin_pola,n_spin_cant,nstep
!
   integer (kind=IntKind) :: id, ie, is, kmax_kkr, num_e, kl, klp, info, zoomDone
   integer (kind=IntKind) :: lkkr0, l, m, lp, mp, sgnS, sgnS_prev, sgnS_save
   integer (kind=IntKind), allocatable :: indx(:), num_ep(:,:)
!
   real (kind=RealKind), intent(in) :: e0, e1, ei
   real (kind=RealKind) :: step_e, step_eh, sgn
!
   complex (kind=CmplxKind) :: e, e_tmp, diff_e, e_head, kappa, det
   complex (kind=CmplxKind) :: a2l, a2lp, d
   complex (kind=CmplxKind), pointer :: pshift(:,:)
   complex (kind=CmplxKind), pointer :: sin_mat(:,:)
   complex (kind=CmplxKind), allocatable :: SinU(:,:)
!
   type PolesStruct
      integer (kind=IntKind) :: sgnS
      real (kind=RealKind) :: step_e
      complex (kind=CmplxKind) :: energy
      type (PolesStruct), pointer :: next
   end type PolesStruct
!
   type (PolesStruct), allocatable, target :: Poles(:,:)
   type (PolesStruct), pointer :: pPoles, pP_prev
!
   step_e = (e1-e0)/real(nstep,kind=RealKind)
   num_e = nstep + 1
!
   lkkr0 = lkkr(1)
   kmax_kkr = (lkkr0+1)**2
   allocate( SinU(1:kmax_kkr,1:kmax_kkr) )
   allocate( indx(1:kmax_kkr) )
   allocate( Poles(LocalNumAtoms,n_spin_pola-n_spin_cant+1), num_EP(LocalNumAtoms,n_spin_pola-n_spin_cant+1) )
!
   open(unit=11,file='SineMatrixDiag.dat',status='unknown',form='formatted')
!
   write(11,'(a,$)') " #Poles   ReE    ImE     Re[Det]      Im[Det] "
   do l = 1,kmax_kkr
      char_lm(1:3) = "ReL"
      write(char_lm(4:6),'(i3)') offset_lm + l
      char_lm(4:4) = '_'
      write(11,'(9x,a6,$)') char_lm
      char_lm(1:3) = "ImL"
      write(char_lm(4:6),'(i3)') offset_lm + l
      char_lm(4:4) = '_'
      write(11,'(9x,a6,$)') char_lm
   enddo
   write(11,'(a)') "  "
!
   num_Ep(:,:) = 0
!
   do is = 1,n_spin_pola-n_spin_cant+1
!     do id = 1,LocalNumAtoms
      do id = 1,1               ! we run only for the 1st atom
         num_ep(id,is) = 0
         pPoles => Poles(id,is)
         nullify( pPoles%next )
         if (lkkr(id) /= lkkr0) then
            kmax_kkr = (lkkr(id)+1)**2
            lkkr0 = lkkr(id)
            deallocate( SinU, indx )
            allocate( SinU(1:kmax_kkr,1:kmax_kkr) )
            allocate( indx(1:kmax_kkr) )
         endif
         do ie = 1, num_e
            e = cmplx(e0+(ie-1)*step_e,ei,kind=CmplxKind)
            if (abs(e) < TEN2m6) then
               e = e + HALF*step_e
            endif
            kappa = sqrt(e)
!
!           ----------------------------------------------------------
            call solveSingleScattering(is,id,e,CZERO)
!           ----------------------------------------------------------
!
            sin_mat => getSineMatrix(is,id)
            do kl = 1, kmax_kkr
               SinU(1:kmax_kkr,kl) = sin_mat(1:kmax_kkr,kl)
            enddo
!
!           ==========================================================
!           rescaling sine matrix
!           ==========================================================
            kl = 0
            a2l = CONE/kappa
            do l = 0, lkkr(id)
               do m = -l, l
                  kl = kl + 1
                  klp = 0
                  a2lp = CONE
                  do lp = 0, lkkr(id)
                     do mp = -lp, lp
                        klp = klp + 1
                        SinU(klp,kl) = SinU(klp,kl)*a2l*a2lp
                     enddo
                     a2lp = a2lp*(2*lp+ONE)/kappa
                  enddo
               enddo
               a2l = a2l*(2*l+ONE)/kappa
            enddo
            if ( .not.isLUdcmp ) then
!              -------------------------------------------------------
               call GaussianElim(SinU,kmax_kkr)
!              -------------------------------------------------------
               det = CONE
               do kl = 1,kmax_kkr
                  det = det*SinU(kl,kl)
               enddo
            else
!              -------------------------------------------------------
               call LUdcmp(SinU,kmax_kkr,kmax_kkr,indx,d)
!              -------------------------------------------------------
               det = d
               do kl = 1, kmax_kkr
                  det = det*SinU(kl,kl)
               enddo
            endif
            sgn = ONE
            sgnS = sign(sgn,real(det,kind=RealKind))
            if ( ie==1 ) then
               sgnS_prev = sgnS
            else if ( sgnS/=sgnS_prev ) then
               num_ep(id,is) = num_ep(id,is)+1
               pPoles%energy = e
               pPoles%sgnS = sgnS
               allocate(pPoles%next)
               pPoles => pPoles%next
               nullify( pPoles%next )
               sgnS_prev = sgnS
            endif
!
            kl = 0
!            a2l = kappa
            a2l = CONE
            do l = 0, lkkr(id)
               do m = -l, l
                  kl = kl + 1
                  klp = 0
                  a2lp = CONE
                  do lp = 0, lkkr(id)
                     do mp = -lp, lp
                        klp = klp + 1
                        SinU(klp,kl) = SinU(klp,kl)*a2l*a2lp
                     enddo
                     a2lp = a2lp*kappa/cmplx((2*lp+ONE),ZERO,kind=CmplxKind)
                  enddo
               enddo
               a2l = a2l*kappa/cmplx((2*l+ONE),ZERO,kind=CmplxKind)
            enddo
            write(11,'(i4,2(2x,2d17.8),$)') num_ep(id,is), e, det
            do kl = 1,kmax_kkr
                write(11,'(2x,2d17.8,$)') SinU(kl,kl)
            enddo
            write(11,'(" ")')
!
         enddo
!
         if ( num_ep(id,is)/=0) then
            pPoles => Poles(id,is)
            step_e = step_e*TWO
            do ie = 1,num_EP(id,is)
               step_eh = step_e/TWO
               e = pPoles%energy
               e_tmp = e - step_eh
               zoomDone = 0
               nullify(sin_mat)
               do while ( zoomDone==0 )
!                 ----------------------------------------------------
                  call solveSingleScattering(is,id,e_tmp,CZERO)
!                 ----------------------------------------------------
                  sin_mat => getSineMatrix(is,id)
                  do kl = 1, kmax_kkr
                     SinU(1:kmax_kkr,kl) = sin_mat(1:kmax_kkr,kl)
                  enddo
!                 ====================================================
!                 rescaling sine matrix
!                 ====================================================
                  kl = 0
                  a2l = CONE/kappa
                  do l = 0, lkkr(id)
                     do m = -l, l
                        kl = kl + 1
                        klp = 0
                        a2lp = CONE
                        do lp = 0, lkkr(id)
                           do mp = -lp, lp
                              klp = klp + 1
                              SinU(klp,kl) = SinU(klp,kl)*a2l*a2lp
                           enddo
                           a2lp = a2lp*(2*lp+ONE)/kappa
                        enddo
                     enddo
                     a2l = a2l*(2*l+ONE)/kappa
                  enddo
                  if ( .not.isLUdcmp ) then
!                    -------------------------------------------------
                     call GaussianElim(SinU,kmax_kkr)
!                    -------------------------------------------------
                     det = CONE
                     do kl = 1,kmax_kkr
                        det = det*SinU(kl,kl)
                     enddo
                  else
!                    -------------------------------------------------
                     call LUdcmp(SinU,kmax_kkr,kmax_kkr,indx,d)
!                    -------------------------------------------------
                     det = d
                     do kl = 1, kmax_kkr
                        det = det*SinU(kl,kl)
                     enddo
                  endif
                  sgn = ONE
                  sgnS = sign(sgn,real(det,kind=RealKind))
                  step_eh = step_eh/TWO
                  if ( sgnS == pPoles%sgnS) then
                     e_tmp = e_tmp - step_eh
!
                     write(6,'(a,3(2x,2d17.8))') " Refine Root  - :",e_tmp, step_eh
                  else
                     e_tmp = e_tmp + step_eh
!
                     write(6,'(a,3(2x,2d17.8))') " Refine Root  + :",e_tmp, step_eh
                  endif
                  if ( step_eh < Ten2m8 ) then
                     write(6,'(a,i4,3(2x,2d17.8))') " Root   :",num_ep(id,is), e_tmp, step_eh
                     zoomDone = 1
                  endif
               enddo
               pPoles%energy = e_tmp
               if ( sgnS == pPoles%sgnS ) then
                  pPoles%step_e = step_eh
               else
                  pPoles%step_e = -step_eh
               endif
               pPoles%sgnS = sgnS
               if ( ie/=num_ep(id,is) ) then
                  pPoles => pPoles%next
               else
                  pP_prev => pPoles
                  pPoles => pPoles%next
                  nullify( pP_prev%next )
                  deallocate(pPoles)
               endif
            enddo
         endif
      enddo
   enddo
!
   do is = 1 ,n_spin_pola-n_spin_cant+1
      do id = 1,LocalNumAtoms
         pPoles => Poles(id,is)%next
         if ( num_Ep(id,is) /=0 ) then
            do ie = 2,num_EP(id,is)
               pP_prev => pPoles%next
               nullify(pPoles%next)
               deallocate(pPoles)
               pPoles => pP_prev
            enddo
         endif
         nullify( Poles(id,is)%next )
      enddo
   enddo
!
   close(unit=11)
   deallocate( Poles, num_EP )
   deallocate( SinU, indx )
!
   end subroutine calScatteringPoles
!  ===================================================================
