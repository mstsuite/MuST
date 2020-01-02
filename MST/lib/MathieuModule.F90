Module MathieuModule
   use KindParamModule, only : IntKind, RealKind
   use ErrorHandlerModule, only : ErrorHandler
!
public :: initMathieu,   &
          endMathieu,    &
          getMathieuBand
!
private
   integer (kind=IntKind), parameter :: MaxBands = 100
!
   integer (kind=IntKind) :: num_degen(MaxBands)
   integer (kind=IntKind) :: num_bands
!
   real (kind=RealKind) :: mathieu_band(MaxBands)
!
   logical :: Initialized = .false.
   logical :: LattModInit = .false.
!
   integer (kind=IntKind) :: print_instruction
!
   character (len=50) :: stop_routine
!
   real (kind=RealKind) :: V1,V2
   real (kind=RealKind) :: pi2oa
   real (kind=RealKind) :: kc
!
contains
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initMathieu(U1,U2,bravais,istop,iprint)
!  ====================================================================
   use MathParamModule, only : pi2, TEN2m6
   use LatticeModule, only : isInitialized
!
   implicit   none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: iprint
!
   real (kind=RealKind), intent(in) :: U1,U2
   real (kind=RealKind), intent(in) :: bravais(3,3)
!
   real (kind=RealKind) :: alat
!
   V1=U1
   V2=U2
!
   alat = bravais(1,1)
   if (alat < TEN2m6) then
      call ErrorHandler('initMathieu','Invalid Bravais lattice vector')
   else if (.not.isInitialized()) then
      call ErrorHandler('initMathieu','LatticeModule needs to be initialized')
   endif
!
   pi2oa=pi2/alat
   Initialized = .true.
   stop_routine = istop
   print_instruction = iprint
!
   end subroutine initMathieu
!  ====================================================================
!
!  ********************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endMathieu
!  ====================================================================
   implicit   none
!
   Initialized = .false.
   end subroutine endMathieu
!  ====================================================================
!
!  ********************************************************************
!
!  cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine getMathieuBand(kvec,emin,emax)
!  ====================================================================
   use MathParamModule, only : ZERO, TWO, TEN2m8, TEN2m14
   use SortModule, only : QuickSort
   use LatticeModule, only : getKSpaceBravais, getKnVectors
!
   implicit   none
!
   integer (kind=IntKind), allocatable :: info(:)
   integer (kind=IntKind) :: i
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: k
   integer (kind=IntKind) :: n
   integer (kind=IntKind) :: m
   integer (kind=IntKind) :: nb(3)
!
   integer (kind=IntKind) :: NumKnVecs
   integer (kind=IntKind) :: NumEpts
!
   real (kind=RealKind), intent(in) :: kvec(3)
   real (kind=RealKind), intent(in) :: emin
   real (kind=RealKind), intent(in) :: emax
!
   real (kind=RealKind), allocatable :: e(:)
   real (kind=RealKind) :: kbas(3,3), kcut
   real (kind=RealKind) :: e1
   real (kind=RealKind) :: e2
   real (kind=RealKind) ::  f1,f2,df
   real (kind=RealKind) :: ec(MaxBands,3)
   real (kind=RealKind) :: e_n(MaxBands)
   real (kind=RealKind) :: zriddr
   real (kind=RealKind), pointer :: kn(:,:)
!
   real (kind=RealKind), parameter :: de=0.001
!
   if (.not.Initialized) then
      call ErrorHandler('getMathieuBand','call initMathieu first!')
   endif
!
   if (abs(V2).lt.ten2m8) then
      NumEpts=1
   else
      NumEpts=int((emax-emin)/de+0.5)
   endif
!
   allocate(info(1:NumEpts),e(1:NumEpts))
!
   kbas(1:3,1:3) = getKSpaceBravais()
!
   kcut = max(kbas(1,1),kbas(2,2),kbas(3,3))
   kcut = abs(kcut)
!
   kn => getKnVectors(kcut)
   NumKnVecs = size(kn,2)
   do k=1,3
      n=0
      do j=1,NumKnVecs
         if (abs(V2).lt.ten2m8) then
            kc=kvec(k)+kn(k,j)
            e(1)=kc*kc
            info(1)=0
         else
!           print *,'j = ',j,'  kn = ',kn(k,j)
!           ----------------------------------------------------------
            call calfg(emin,f1,df)
!           ----------------------------------------------------------
            if (abs(f1).lt.ten2m8) then
               e(1)=emin
               info(1)=0
            endif
            do i=1,NumEpts
               e1=emin+de*(i-1)
               e2=e1+de
               kc=kvec(k)+kn(k,j)
               call calfg(e2,f2,df)
               if (abs(f2).lt.ten2m8) then
                  e(i)=e2
                  info(i)=0
               else if (abs(f1).gt.ten2m8 .and. f1*f2.lt.zero) then
!                 ----------------------------------------------------
                  e(i)=zriddr(calfg,e1,e2,ten2m14,info(i))
!                 ----------------------------------------------------
!                 write(6,'(3f10.5,1i5)')e1,e2,e(i),info(i)
               else
                  e(i)=zero
                  info(i)=1
               endif
               f1=f2
            enddo
         endif
         do i=1,NumEpts
            if (info(i).eq.0) then
               m=1
               do while(m.le.n)
                  if(abs(e(i)-ec(m,k)).lt.ten2m8) then
                     m=0
                     exit
                  endif
                  m=m+1
               enddo
               if(m.gt.0) then
                  if (n.eq.MaxBands) then
                     call ErrorHandler('getMathieuBand','n = MaxBands')
                  else
                     n=n+1
                     ec(n,k)=e(i)
!                    print *,n,'  ',e(i)
                  endif
               endif
            endif
         enddo
      enddo
      nb(k)=n
   enddo
   deallocate(info,e)
!
   n=nb(1)*nb(2)*nb(3)
   if(n .gt. MaxBands) then
      call ErrorHandler('getMathieuBand','nb(1)*nb(2)*nb(3) > MaxBands',&
                        n,MaxBands)
   else if(n .eq. 0) then
!     call ErrorHandler('getMathieuBand','no band exists')
      nb(1)=1; nb(2)=1; nb(3)=1
      ec(1,1)=zero; ec(1,2)=zero; ec(1,3)=zero
   endif
!
   n=0
   do i=1,nb(1)
      do j=1,nb(2)
         do k=1,nb(3)
            n=n+1
            e_n(n)=V1+ec(i,1)+ec(j,2)+ec(k,3)
         enddo
      enddo
   enddo
!  -------------------------------------------------------------------
   call QuickSort(n,e_n)
!  -------------------------------------------------------------------
   m=0
   do i=1,n
      if(e_n(i).ge.emin .and. e_n(i).le.emax) then
         m=m+1
         mathieu_band(m)=e_n(i)
      endif
   enddo
   n=m
   if(n.eq.0) then
      call ErrorHandler('getMathieuBand','no band exists between',emin,emax)
   endif
   do i=1,n
      e_n(i)=mathieu_band(i)
   enddo
!
   num_bands=1
   mathieu_band(1)=e_n(1)
   num_degen(1)=1
   do i=2,n
      if(abs(e_n(i)-mathieu_band(num_bands)).lt.ten2m8) then
         num_degen(num_bands)=num_degen(num_bands)+1
      else
         num_bands=num_bands+1
         mathieu_band(num_bands)=e_n(i)
         num_degen(num_bands)=1
      endif
   enddo
   write(6,'('' K-vector = '',3f10.5)')kvec(1:3)
   write(6,'('' =================================================='')')
   write(6,'(''  Band Index        Band Energy         Degeneracy '')')
   write(6,'('' --------------------------------------------------'')')
   do j=1,num_bands
      write(6,'(3x,1i5,8x,f15.8,14x,1i1)')j,mathieu_band(j),num_degen(j)
   enddo
   write(6,'(''===================================================='',/)')
   end subroutine getMathieuBand
!  ===================================================================
!
!  *******************************************************************
!
!  ===================================================================
   subroutine calfg(E,FGm1,dFG)
   use MathParamModule, only : zero, one, two
!
   implicit   none
!
   integer (kind=IntKind), parameter :: nmax=50
   integer (kind=IntKind) :: n
!
   real (kind=RealKind) :: E
   real (kind=RealKind) :: FGm1
   real (kind=RealKind) :: dFG
!
   real (kind=RealKind) :: ToV2
   real (kind=RealKind) :: Gn
   real (kind=RealKind) :: Gnp
   real (kind=RealKind) :: Fn
   real (kind=RealKind) :: Fnp
   real (kind=RealKind) :: Sn
   real (kind=RealKind) :: Tn
   real (kind=RealKind) :: dE
!
!  *******************************************************************
!  Solve the one dimensional Mathieu equation  
!
!     y'' + [ E - V(x)]y=0   where V(x)= V2 * cos(2*pi*x/a0)
!                                                    
!  input, kc   : k-vector component
!         V2   : constant potential
!         pi2oa: pi*2/a0
!         E    : initial guess of eigenvalue
!    
!  output,E    : final value of eigenvalue
!  *******************************************************************
!
!  ===================================================================
!  Let Sn = -2*(E-Kn*Kn)/V2 and Tn = -2*(E-K(-n)*K(-n))/V2,
!  evaluate the following continued fraction to obtain Gn = Cn/Cn-1, 
!  and Fn = C-n-1/C-n:
!
!      Sn * Cn  + Cn+1  + Cn-1  = 0,  or,  Sn + Gn+1 + 1/Gn   = 0
!      Tn * C-n + C-n-1 + C-n+1 = 0,  or,  Tn + Fn   + 1/Fn-1 = 0
!
!  A valid solution (k and E consistent) must have G0F0=1.
!
!  Gnp and Fnp are energy derivatives of Gn and Fn respectively 
!
!      Sn' + Gn+1' -  Gn' / Gn**2  = 0,    Sn' = -2/V2
!      Tn' + Fn'   - Fn-1'/Fn-1**2 = 0,    Tn' = -2/V2
!
!  They are used for Newton-Rapheson improvement.
!  ===================================================================
   ToV2=two/V2
   Gn= zero
   Gnp=zero
   Fn= zero
   Fnp=zero
   do n=nmax,0,-1
      Sn = -ToV2*( E - (kc +  n  *pi2oa)**2 )
      Tn = -ToV2*( E - (kc -(n+1)*pi2oa)**2 )
      Gn = -one/(Sn+Gn)
      Fn = -one/(Tn+Fn)
      Gnp= Gn*Gn*(Gnp-ToV2)
      Fnp= Fn*Fn*(Fnp-ToV2)
   enddo
   FGm1=Gn*Fn-one
   dFG=Fn*Gnp+Fnp*Gn
!
   end subroutine calfg
!  ===================================================================
end module MathieuModule
