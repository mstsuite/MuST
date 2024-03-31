module AccelerateCPAModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : ZERO, ONE, CZERO, CONE, SQRTm1, TEN2m8
   use ErrorHandlerModule, only : ErrorHandler, WarningHandler
   use PublicParamDefinitionsModule, only : SimpleMixing, AndersonMixing, &
                                            BroydenMixing, AndersonMixingOld
!
public :: initAccelerateCPA,      &
          setAccelerationParam,   &
          initializeAcceleration, &
          getAccelerationType,    &
          accelerateCPA,          &
          endAccelerateCPA
!
private
   integer (kind=IntKind) :: MaxIteration = 30
   integer (kind=IntKind) :: Kmax_KKR_Save
!
   real (kind=RealKind) :: alpha = 0.15d0
   real (kind=RealKind) :: CPA_tol = TEN2m8
!
!  ===========
   integer (kind=IntKind), parameter :: cdim=30
   integer (kind=IntKind), parameter :: iscf_cpa=1
   real (kind=RealKind), parameter :: cw0=5.0d-03
!
   integer (kind=IntKind) :: cpaiter_m, cpaiter_k
   integer (kind=IntKind) :: AccelerationType
!
   real (kind=RealKind) :: cpaiter_pq, cpaiter_p, cpaiter_ppq
   real (kind=RealKind) :: cpaiter_tp, cpaiter_u
   real (kind=RealKind), allocatable, target :: RWORK(:)
   real (kind=RealKind), pointer :: cpaiter_x(:)
   real (kind=RealKind), allocatable :: cpaiter_t(:)
   real (kind=RealKind), allocatable :: cpaiter_nm(:), cpaiter_nml(:)
   real (kind=RealKind), allocatable :: cpaiter_fm(:), cpaiter_fml(:)
   real (kind=RealKind), allocatable :: cpaiter_delta(:,:,:)
   real (kind=RealKind), allocatable :: cpaiter_bkni(:,:)
!
!  ===========
   integer (kind=IntKind), parameter :: ipits = 4
!
   complex (kind=CmplxKind), allocatable, target :: CWORK1(:)
   complex (kind=CmplxKind), allocatable, target :: CWORK2(:)
   complex (kind=CmplxKind), pointer :: tcin(:,:)
   complex (kind=CmplxKind), pointer :: tcout(:,:)
!
!  ===========
   complex (kind=CmplxKind), pointer :: tc_save(:)
!
contains
!
   include '../lib/arrayTools.F90'
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initAccelerateCPA(acc_type,max_iter,acc_mix,ctol,kmax_kkr)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: acc_type, kmax_kkr
   integer (kind=IntKind), optional, intent(in) :: max_iter
   integer (kind=IntKind) :: i, n
!
   real (kind=RealKind), optional, intent(in) :: acc_mix, ctol
!
   AccelerationType = acc_type
   Kmax_KKR_Save = kmax_kkr
!
   if (present(acc_mix)) then
      alpha = acc_mix
   endif
!
   if (present(ctol)) then
      CPA_tol = ctol
   endif
!
   if (present(max_iter)) then
      MaxIteration = max_iter
   endif
!
   if (AccelerationType == BroydenMixing) then
      n = 2*kmax_kkr**2
      if (iscf_cpa.eq.1) then
         allocate (cpaiter_nm(n),cpaiter_fm(n),cpaiter_nml(n),cpaiter_fml(n))
         if(cdim.gt.1) allocate (cpaiter_delta(n,cdim,2),cpaiter_bkni(cdim,cdim))
      else if(iscf_cpa.eq.2) then
         i=2*n*(1+cdim)+(1+cdim)**2+2*cdim**2
         allocate (cpaiter_t(i))
      endif
      allocate (RWORK(n))
      cpaiter_x => RWORK
   else if (AccelerationType == AndersonMixing .or. AccelerationType == AndersonMixingOld) then
      n = kmax_kkr**2*ipits
      allocate (CWORK1(n))
      allocate (CWORK2(n))
      tcin => aliasArray2_c(CWORK1,kmax_kkr**2,ipits)
      tcout => aliasArray2_c(CWORK2,kmax_kkr**2,ipits)
   else
      n = kmax_kkr**2
      allocate (CWORK1(n))
      tc_save => CWORK1
   endif
!
   end subroutine initAccelerateCPA
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endAccelerateCPA()
!  ===================================================================
   implicit none
!
   if (AccelerationType == BroydenMixing) then
      if (iscf_cpa == 1) then
         deallocate(cpaiter_nm)
         deallocate(cpaiter_nml)
         deallocate(cpaiter_fm)
         deallocate(cpaiter_fml)
         if (cdim > 1) then
            deallocate(cpaiter_delta)
            deallocate(cpaiter_bkni)
         endif
      else if (iscf_cpa == 2) then
         deallocate(cpaiter_t)
      endif
      deallocate(RWORK)
      nullify(cpaiter_x)
   else if (AccelerationType == AndersonMixing .or. AccelerationType == AndersonMixingOld) then
      nullify(tcin, tcout)
      deallocate (CWORK1)
      deallocate (CWORK2)
   else
      nullify(tc_save)
      deallocate(CWORK1)
   endif
!
   end subroutine endAccelerateCPA
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setAccelerationParam(acc_mix,acc_type,max_iter)
!  ===================================================================
   implicit none
!
   real (kind=RealKind), intent(in) :: acc_mix
   integer (kind=IntKind), intent(in), optional :: acc_type
   integer (kind=IntKind), intent(in), optional :: max_iter
   integer (kind=IntKind) :: n, kmax_kkr
!
   if (present(acc_type)) then
      if (AccelerationType /= acc_type) then
         call endAccelerateCPA()
         kmax_kkr = Kmax_KKR_Save
         call initAccelerateCPA(acc_type=acc_type,acc_mix=acc_mix,kmax_kkr=kmax_kkr)
      endif
   endif
!
   if (present(max_iter)) then
      MaxIteration = max_iter
   endif
!
   alpha = acc_mix
!
   end subroutine setAccelerationParam
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getAccelerationType() result(a)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: a
!
   a = AccelerationType
!  
   end function getAccelerationType
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initializeAcceleration(tcpa,ndim,itcpa)
!  ===================================================================
   implicit none 
!
   integer (kind=IntKind), intent(in) :: itcpa, ndim
   integer (kind=IntKind) :: i, n
!
   complex (kind=CmplxKind), intent(in) :: tcpa(ndim)
!
   if (AccelerationType == BroydenMixing .and. itcpa == 1) then
      n=2*ndim ! calculate length of the vector.......................
      cpaiter_x(1:ndim)  = real(tcpa(1:ndim))
      cpaiter_x(1+ndim:n)= aimag(tcpa(1:ndim))
      if (iscf_cpa == 1) then
         cpaiter_nm(1:n)=cpaiter_x(1:n)
         cpaiter_fm(1:n)=ZERO
         cpaiter_delta=ZERO
         cpaiter_bkni=ZERO
      else if(iscf_cpa == 2) then
         cpaiter_t=ZERO
         cpaiter_x(1:ndim)  = real(tcpa(1:ndim))
         cpaiter_x(1+ndim:n)= aimag(tcpa(1:ndim))
         do i=1,n
            cpaiter_t(i)=cpaiter_x(i)
         enddo
         cpaiter_ppq=max(cpaiter_ppq,cw0)
      endif
      cpaiter_k=0
      cpaiter_m=0
      cpaiter_pq=ONE
      cpaiter_p=alpha
      cpaiter_ppq=alpha
      cpaiter_tp=ZERO
      cpaiter_u=ZERO
   else if (AccelerationType == AndersonMixing .or. AccelerationType == AndersonMixingOld) then
      i=mod(itcpa-1,ipits)+1
      call zcopy(ndim,tcpa,1,tcin(1,i),1)
   else if (AccelerationType == SimpleMixing) then
      call zcopy(ndim,tcpa,1,tc_save,1)
   endif
!
   end subroutine initializeAcceleration
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine accelerateCPA(tcpa,ndim,itcpa)
!  ===================================================================
   implicit none 
!
   integer (kind=IntKind), intent(in) :: itcpa, ndim
   integer (kind=IntKind) :: i
!
   complex (kind=CmplxKind), intent(inout) :: tcpa(ndim)
   complex (kind=CmplxKind) :: cfac
!
   if (AccelerationType == BroydenMixing) then
!     ----------------------------------------------------------------
      call cpaiter(tcpa,ndim)
!     ----------------------------------------------------------------
   else if (AccelerationType == AndersonMixing) then
      i=mod(itcpa-1,ipits)+1
!     ----------------------------------------------------------------
      call zcopy(ndim,tcpa,1,tcout(1,i),1)
!     ----------------------------------------------------------------
!     
!     ================================================================
!     The following lines of code are modified on 5/22/2020.
!     The new code will start calling acctc at itcpa =1, while the old
!     code calls acctc starting itcpa > ipits
!     ----------------------------------------------------------------
      call acctc(tcpa,ndim,itcpa)
!     ----------------------------------------------------------------
   else if (AccelerationType == AndersonMixingOld) then
!     ================================================================
!     The code before 5/22/2020 is as follows. I now call it AndersonMixingOld
!     ================================================================
      i=mod(itcpa-1,ipits)+1
!     ----------------------------------------------------------------
      call zcopy(ndim,tcpa,1,tcout(1,i),1)
!     ----------------------------------------------------------------
!
      if ( itcpa >= ipits ) then
!        =============================================================
!        speed up the convergence by loading the accelarator.
!        -------------------------------------------------------------
         call acctc(tcpa,ndim,itcpa)
!        -------------------------------------------------------------
      else
         cfac = alpha
         tcpa = cfac*tcpa
         cfac = CONE-cfac
!        -------------------------------------------------------------
         call zaxpy(ndim,cfac,tcin(1,i),1,tcpa,1)
!        -------------------------------------------------------------
      endif
   else
      cfac = alpha
      tcpa = cfac*tcpa
      cfac = CONE-cfac
!     ----------------------------------------------------------------
      call zaxpy(ndim,cfac,tc_save,1,tcpa,1)
!     ----------------------------------------------------------------
   endif
!
   end subroutine accelerateCPA
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine cpaiter(matrix,ndim)
!  ===================================================================
   implicit none 
!
   integer (kind=IntKind), intent(in) :: ndim
   integer (kind=IntKind) :: i,n
!
   complex (kind=CmplxKind), intent(inout) :: matrix(ndim)
!
   interface
      subroutine giterat(niter,mdim,pmix,x,n,t,u,m,k,p,ppq,pq,w0,tp,tolerance)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind) :: niter,mdim,n,m,k
         real (kind=RealKind) :: x(:),t(:),pmix,u,p,ppq,pq,w0,tp,tolerance
      end subroutine giterat
   end interface
!
   interface
      subroutine giterbr(pmix,w0,mdim,x,nx,m,nm,nml,fm,fml,delta,bkni,tp)
         use KindParamModule, only : IntKind, RealKind
         implicit none
         integer (kind=IntKind) :: nx,mdim,m
         real (kind=RealKind) :: x(nx),pmix,w0,nm(nx),nml(nx),fm(nx),fml(nx),&
       &                         delta(nx,mdim,2),bkni(mdim,mdim),tp
      end subroutine giterbr
   end interface
!
   n=2*ndim
   cpaiter_x(1:ndim)   = real(matrix(1:ndim))
   cpaiter_x(1+ndim:n) = aimag(matrix(1:ndim))
   if(iscf_cpa == 1) then
      call giterbr(alpha,cw0,cdim,cpaiter_x,n,cpaiter_m,cpaiter_nm,cpaiter_nml, &
                   cpaiter_fm,cpaiter_fml,cpaiter_delta,cpaiter_bkni,cpaiter_tp)
   else if(iscf_cpa == 2) then
      call giterat(MaxIteration,cdim,alpha,cpaiter_x,n,cpaiter_t,cpaiter_u, &
                   cpaiter_m,cpaiter_k,cpaiter_p,cpaiter_ppq,cpaiter_pq,cw0,&
                   cpaiter_tp,CPA_tol)
   endif
   matrix(1:ndim)=cmplx(cpaiter_x(1:ndim),cpaiter_x(ndim+1:n),kind=CmplxKind)
!
   end subroutine cpaiter
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine acctc(tc,ndim,nt)
!  ===================================================================
   implicit   none
!
   integer (kind=IntKind), intent(in) :: ndim
   integer (kind=IntKind), intent(in) :: nt
   integer (kind=IntKind) :: nit
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: info
!
   real (kind=RealKind) :: dptc
   real (kind=RealKind) :: acctc_a(ipits+1,ipits+1)
   real (kind=RealKind) :: acctc_b(ipits+1)
!
   complex (kind=CmplxKind), intent(out) :: tc(ndim)
   complex (kind=CmplxKind) :: sum, sum_old
!
!  interface
!     function dptc(tcini,tcouti,tcinj,tcoutj,ndim) result(d)
!        use KindParamModule, only : IntKind, RealKind
!        integer (kind=IntKind), intent(in) :: ndim
!        real (kind=RealKind), intent(in) :: tcini(ndim)
!        real (kind=RealKind), intent(in) :: tcouti(ndim)
!        real (kind=RealKind), intent(in) :: tcinj(ndim)
!        real (kind=RealKind), intent(in) :: tcoutj(ndim)
!        real (kind=RealKind) :: d
!     end function dptc
!  end interface
!
   interface
      subroutine dgaleq(a,y,n,ipits,info)
         use KindParamModule, only : IntKind, RealKind
         integer (kind=IntKind), intent(in) :: n
         integer (kind=IntKind), intent(in) :: ipits
         integer (kind=IntKind), intent(out) :: info
         real (kind=RealKind), intent(inout) :: a(ipits+1,ipits+1)
         real (kind=RealKind), intent(inout) :: y(ipits+1)
      end subroutine dgaleq
   end interface
!
   acctc_a = ZERO
   acctc_b = ZERO
   nit=min(ipits,nt)
   n=nit+1
   do j=1,nit
      do i=1,j
         acctc_a(i,j)=dptc(tcin(1,i),tcout(1,i),tcin(1,j),tcout(1,j),ndim*2)
         if(i.ne.j) acctc_a(j,i)=acctc_a(i,j)
      enddo
      acctc_b(j)=ZERO
      acctc_a(j,n)=ONE
      acctc_a(n,j)=ONE
   enddo
   acctc_a(n,n)=ZERO
   acctc_b(n)=ONE
!
!  -------------------------------------------------------------------
   call dgaleq(acctc_a,acctc_b,n,ipits,info)
!  -------------------------------------------------------------------
   if (info == 0) then
      if (alpha > 0.999 .or. iscf_cpa == 2) then
         do k=1,ndim
            sum=CZERO
            do i=1,nit
               sum=sum+acctc_b(i)*tcout(k,i)
            enddo
            tc(k)=sum
         enddo
      else
         do k=1,ndim
            sum_old=CZERO
            sum=CZERO
            do i=1,nit
               sum_old=sum_old+acctc_b(i)*tcin(k,i)
               sum=sum+acctc_b(i)*tcout(k,i)
            enddo
            tc(k)=sum_old+alpha*(sum-sum_old)
         enddo
      endif
   else ! Use simple mixing to get new tc
!     ----------------------------------------------------------------
!     call WarningHandler('acctc','Ill condition appeared in DGA mixing',force_to_print=.true.)
!     ----------------------------------------------------------------
      n = mod(nt-1,ipits)+1
      do k=1,ndim
         tc(k)=tcin(k,n)+alpha*(tcout(k,n)-tcin(k,n))
      enddo
   endif
!
   end subroutine acctc
!  ===================================================================
end module AccelerateCPAModule
