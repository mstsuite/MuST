!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine computeTraceComplex(isExchangeParamNeeded,                &
                                  kkrsz,tmat_g,tau00_g,                 &
                                  tr_pxtau,tr_pmtau,tr_pdott,tr_ptpt00, &
                                  tr_ptptab)
!  ===================================================================
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use MathParamModule, only : CZERO, SQRTm1, FOURTH
!
   implicit none
!
   logical, intent(in) :: isExchangeParamNeeded
!
   integer (kind=IntKind), intent(in) :: kkrsz
!
   complex (kind=CmplxKind), intent(in) :: tmat_g(:,:)
   complex (kind=CmplxKind), intent(in) :: tau00_g(:,:)
!
   complex (kind=CmplxKind), intent(out) :: tr_pxtau(3)
   complex (kind=CmplxKind), intent(out) :: tr_pmtau(3)
   complex (kind=CmplxKind), intent(out) :: tr_pdott
   complex (kind=CmplxKind), intent(out) :: tr_ptpt00
   complex (kind=CmplxKind), intent(out) :: tr_ptptab(9)
!
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: info
   integer (kind=IntKind), allocatable :: ipvt(:)
   integer (kind=IntKind) :: ka(9)
   integer (kind=IntKind) :: kb(9)
!
   data (ka(i),i=1,9)/1,2,3,1,2,3,1,2,3/
   data (kb(i),i=1,9)/1,1,1,2,2,2,3,3,3/
!
   complex (kind=CmplxKind), allocatable :: pmat(:,:)
   complex (kind=CmplxKind), allocatable :: pmat_m(:,:)
   complex (kind=CmplxKind), allocatable :: p_vec(:,:,:)
   complex (kind=CmplxKind), allocatable :: tau_vec(:,:,:)
   complex (kind=CmplxKind), allocatable :: tau_0(:,:)
   complex (kind=CmplxKind), allocatable :: wspace0(:)
   complex (kind=CmplxKind), allocatable :: wspace1(:,:)
   complex (kind=CmplxKind), allocatable :: wspace2(:,:)
   complex (kind=CmplxKind) :: term
!
   allocate( ipvt(2*kkrsz) )
   allocate( pmat(2*kkrsz,2*kkrsz) )
   allocate( pmat_m(kkrsz,kkrsz) )
   allocate( p_vec(kkrsz,kkrsz,3) )
   allocate( tau_vec(kkrsz,kkrsz,3) )
   allocate( tau_0(kkrsz,kkrsz) )
   allocate( wspace0(4*kkrsz*kkrsz) )
   allocate( wspace1(kkrsz,kkrsz) )
   allocate( wspace2(kkrsz,kkrsz) )
!
!  ===================================================================
!  call zgetrf and zgetri to invert tmat_g........................
!  -------------------------------------------------------------------
   call zcopy(4*kkrsz*kkrsz,tmat_g,1,pmat,1)
   call zgetrf(kkrsz*2,kkrsz*2,pmat,kkrsz*2,ipvt,info)
   call zgetri(kkrsz*2,pmat,kkrsz*2,ipvt,wspace0,4*kkrsz*kkrsz,info)
!  -------------------------------------------------------------------
!
!  ===================================================================
!  calculate tr_pxtau, tr_pmtau, tr_pdott,tr_ptpt00, tr_ptptab,
!  and p_vec, tau_vec
!     
!
!                          ^  00
!     input:    tau00_g  = Tau                    in the global frame
!
!                          ^   -1   ^
!               pmat     = Tmat   = Pmat          in the global frame         
!
!                               ^      ->        -> 
!               pmat_m   = Tr [ Pmat * Sigma ] * e
!                            s
!
!                               ->  ->
!     output:   tr_pxtau = Tr [ p X Tau ] 
!                            L                             
!
!                                   ->
!               tr_pmtau = Tr [ p * Tau ] 
!                            L   
!
!                               ->  ->
!               tr_pdott = Tr [ p * Tau ]
!                            L
!
!
!               tr_ptpt00= Tr [ p * Tau * p * Tau ]
!                            L         0         0
!
!
!               tr_ptptab= Tr [ p * Tau * p * Tau  ]
!                            L         a         b
!
!                               ^      ->        ->
!     where     p        = Tr [ Pmat * Sigma ] * e / 2 = pmat_m / 2
!                            s
!
!               ->              ^      ->
!               p        = Tr [ Pmat * Sigma ] / 2     = p_vec / 2
!                            s
!
!               ->              ^  00  ->
!               Tau      = Tr [ Tau  * Sigma ] / 2     = tau_vec / 2
!                            s
!
!                        = (Tau , Tau , Tau )
!                              1     2     3
!
!                               ^  00 
!               Tau      = Tr [ Tau   ] / 2  = tau_0 / 2
!                  0         s
!
!               a, b     = 1, 2, 3, or x, y, z
!
!               d        = 1, if a =  b
!                ab
!                        = 0, if a <> b
!  ===================================================================
!
!  -------------------------------------------------------------------
   tr_pxtau(1:3) = CZERO
   tr_pmtau(1:3) = CZERO
   tr_ptptab(1:9) = CZERO
   tr_pdott = CZERO
   tr_ptpt00 = CZERO
!  -------------------------------------------------------------------
!
      do j=1,kkrsz
	 do i=1,kkrsz
	    p_vec(i,j,1)=pmat(kkrsz+i,j)
         enddo
	 do i=1,kkrsz
	    p_vec(i,j,1)=p_vec(i,j,1)+pmat(i,kkrsz+j)
         enddo
      enddo
      do j=1,kkrsz
	 do i=1,kkrsz
	    tau_vec(i,j,1)=tau00_g(kkrsz+i,j)
         enddo
	 do i=1,kkrsz
	    tau_vec(i,j,1)=tau_vec(i,j,1)+tau00_g(i,kkrsz+j)
         enddo
      enddo
!
      do j=1,kkrsz
	 do i=1,kkrsz
	    p_vec(i,j,2)=-pmat(kkrsz+i,j)
         enddo
	 do i=1,kkrsz
	    p_vec(i,j,2)=sqrtm1*(p_vec(i,j,2)+pmat(i,kkrsz+j))
         enddo
      enddo
      do j=1,kkrsz
	 do i=1,kkrsz
	    tau_vec(i,j,2)=-tau00_g(kkrsz+i,j)
         enddo
	 do i=1,kkrsz
	    tau_vec(i,j,2)=sqrtm1*(tau_vec(i,j,2)+tau00_g(i,kkrsz+j))
         enddo
      enddo
!
      do j=1,kkrsz
	 do i=1,kkrsz
	    p_vec(i,j,3)=pmat(i,j)
         enddo
	 do i=1,kkrsz
	    p_vec(i,j,3)=p_vec(i,j,3)-pmat(kkrsz+i,kkrsz+j)
         enddo
      enddo
      do j=1,kkrsz
	 do i=1,kkrsz
	    tau_vec(i,j,3)=tau00_g(i,j)
         enddo
	 do i=1,kkrsz
	    tau_vec(i,j,3)=tau_vec(i,j,3)-tau00_g(kkrsz+i,kkrsz+j)
         enddo
      enddo
!
      do j=1,kkrsz
	 do i=1,kkrsz
            tr_pxtau(1)=tr_pxtau(1)+p_vec(i,j,2)*tau_vec(j,i,3)-      &
                                    p_vec(i,j,3)*tau_vec(j,i,2)
         enddo
      enddo
      do j=1,kkrsz
	 do i=1,kkrsz
            tr_pxtau(2)=tr_pxtau(2)+p_vec(i,j,3)*tau_vec(j,i,1)-      &
                                    p_vec(i,j,1)*tau_vec(j,i,3)
         enddo
      enddo
      do j=1,kkrsz
	 do i=1,kkrsz
            tr_pxtau(3)=tr_pxtau(3)+p_vec(i,j,1)*tau_vec(j,i,2)-      &
                                    p_vec(i,j,2)*tau_vec(j,i,1)
         enddo
      enddo
      tr_pxtau(1)=fourth*tr_pxtau(1)
      tr_pxtau(2)=fourth*tr_pxtau(2)
      tr_pxtau(3)=fourth*tr_pxtau(3)
!
      if(.not.isExchangeParamNeeded) then
         return     ! only torque needs to be calculated
      endif
!
      do j=1,kkrsz
	 do i=1,kkrsz
	    tau_0(i,j)=tau00_g(i,j)
         enddo
	 do i=1,kkrsz
	    tau_0(i,j)=tau_0(i,j)+tau00_g(kkrsz+i,kkrsz+j)
         enddo
      enddo
!
      do j=1,kkrsz
	 do i=1,kkrsz
            tr_pmtau(1)=tr_pmtau(1)+pmat_m(i,j)*tau_vec(j,i,1)
         enddo
      enddo
      do j=1,kkrsz
	 do i=1,kkrsz
            tr_pmtau(2)=tr_pmtau(2)+pmat_m(i,j)*tau_vec(j,i,2)
         enddo
      enddo
      do j=1,kkrsz
	 do i=1,kkrsz
            tr_pmtau(3)=tr_pmtau(3)+pmat_m(i,j)*tau_vec(j,i,3)
         enddo
      enddo
      tr_pmtau(1)=fourth*tr_pmtau(1)
      tr_pmtau(2)=fourth*tr_pmtau(2)
      tr_pmtau(3)=fourth*tr_pmtau(3)
!
      term=czero
      do k=1,3
         do j=1,kkrsz
            do i=1,kkrsz
               term=term+p_vec(i,j,k)*tau_vec(j,i,k)
            enddo
         enddo
      enddo
      tr_pdott=fourth*term
!
      do j=1,kkrsz
         do i=1,kkrsz
            wspace1(i,j)=CZERO
         enddo
         do n=1,kkrsz
            do i=1,kkrsz
               wspace1(i,j)=wspace1(i,j)+pmat_m(i,n)*tau_0(n,j)
            enddo
         enddo
      enddo
      term=czero
      do j=1,kkrsz
         do i=1,kkrsz
            term=term+wspace1(i,j)*wspace1(j,i)
         enddo
      enddo
      tr_ptpt00=fourth*fourth*term
!
      do k=1,9
         do j=1,kkrsz
            do i=1,kkrsz
               wspace1(i,j)=CZERO
            enddo
            do i=1,kkrsz
               wspace2(i,j)=CZERO
            enddo
            do n=1,kkrsz
               do i=1,kkrsz
                  wspace1(i,j)=wspace1(i,j)+pmat_m(i,n)*tau_vec(n,j,ka(k))
               enddo
               do i=1,kkrsz
                  wspace2(i,j)=wspace2(i,j)+pmat_m(i,n)*tau_vec(n,j,kb(k))
               enddo
            enddo
         enddo
         term=czero
         do j=1,kkrsz
            do i=1,kkrsz
               term=term+wspace1(i,j)*wspace2(j,i)
            enddo
         enddo
         tr_ptptab(k)=fourth*fourth*term
      enddo
!
   deallocate( ipvt,pmat,pmat_m,p_vec,tau_vec,tau_0,wspace0,wspace1,wspace2 )
!
   end subroutine computeTraceComplex
