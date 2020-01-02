! ********************************************************************
! *  Based on LloydModule                                            *
! ********************************************************************
! stores quantities used for the Krein's fomula
! ********************************************************************
!
module KreinModule
  use KindParamModule, only : IntKind, RealKind, CmplxKind
  use MathParamModule, only : CZERO, TEN, TEN2m8, TEN2m7, TEN2m10
!  
  public :: initKrein,       &
            endKrein,        &
            calJost_E0,      &
            setLloydStatus,  &
            setLloydQ,       &
            isLloydOn,       &
            isLloydOnly,     &
            getLloydQ,       &
            aimag_tol,       &
            print_Krein_Info
!
  real (kind=RealKind), parameter :: aimag_tol = TEN2m7 ! CZERO ! TEN2m8
!
  complex (kind=CmplxKind) :: Krein_E0
!
  type jost_E0_Struct
     integer (kind=IntKind)                 :: kmax_kkr
     complex (kind=CmplxKind), allocatable  :: jost_at_E0(:,:,:)
     complex (kind=CmplxKind), allocatable  :: log_j_E0(:,:)
  end type jost_E0_Struct
  type (jost_E0_Struct), allocatable,public :: jost_E0(:)
!
  interface setLloydQ
     module procedure setLloydQ1, setLloydQ2, setLloydQ3
  end interface
!
private
!
  logical :: Initialized = .false.
  logical :: isJostE0_done = .false.
  logical :: lloydOn = .false.
  logical :: lloydOnly = .false.
!
  integer (kind=IntKind) :: LocalNumAtoms
  integer (kind=IntKind) :: nSpinCant, nSpinPola
  integer (kind=IntKind) :: NumKs
  integer (kind=IntKind) :: lmax_kkr_max, kmax_kkr_max
!
  complex (kind=CmplxKind) :: q_lloyd(2)
!
contains
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initKrein(nla, cant, pola, lmax_kkr, nk, istop)
!  ===================================================================
   use ScfDataModule,only: ErBottom
   use ScfDataModule              , only : isLloyd
   implicit none

   character (len=*), intent(in) :: istop 
!
   integer (kind=IntKind), intent(in) :: nla, cant, pola
   integer (kind=IntKind), intent(in) :: lmax_kkr(nla)
   integer (kind=IntKind), intent(in) :: nk
!
   integer(kind=IntKind) :: is,j,kmax,ns
!
   if ( .not.isLloyd() ) return
!
   if ( Initialized ) call endKrein()
!
   LocalNumAtoms=nla
   nSpinPola=pola
   nSpinCant=cant
   NumKs=nk

   allocate (jost_E0(LocalNumAtoms))
   ns = nSpinPola
   do j = 1,LocalNumAtoms
      kmax =(lmax_kkr(j)+1)**2
      jost_E0(j)%kmax_kkr = kmax
      allocate(jost_E0(j)%jost_at_E0(kmax,kmax,ns))
      allocate(jost_E0(j)%log_j_E0(kmax,ns))
      jost_E0(j)%jost_at_E0 = CZERO
      jost_E0(j)%log_j_E0 = CZERO
   enddo
   q_lloyd = CZERO
!    
   Krein_E0=cmplx(ErBottom,aimag_tol,kind=CmplxKind)
!
   Initialized = .true.
   write (6,'(a)') "initializing LLoyd module done"
!
   call print_Krein_Info()
!
   end subroutine initKrein
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endKrein()
!  ===================================================================
   implicit none
   integer(kind=IntKind) :: j
!
   if (allocated(jost_E0)) then
      do j = 1,LocalNumAtoms
         if (allocated(jost_E0(j)%jost_at_E0)) then
            deallocate(jost_E0(j)%jost_at_E0)
         end if
         if (allocated(jost_E0(j)%log_j_E0)) then
            deallocate(jost_E0(j)%log_j_E0)
         end if
      enddo
      deallocate(jost_E0)
   endif
!
   isJostE0_done = .false.
   Initialized = .false.
!
   end subroutine endKrein
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setLloydStatus(y1,y2)
!  ===================================================================
   implicit none
!
   logical, intent(in) :: y1, y2
!
   LloydOn = y1
   LloydOnly = y2
!
   end subroutine setLloydStatus
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setLloydQ1(v)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind), intent(in) :: v
!
   q_lloyd = v
!
   end subroutine setLloydQ1
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setLloydQ2(is,v)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: is
!
   complex (kind=CmplxKind), intent(in) :: v
!
   q_lloyd(is) = v
!
   end subroutine setLloydQ2
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine setLloydQ3(v,n)
!  ===================================================================
   implicit none
!
   integer (kind=IntKind), intent(in) :: n
!
   complex (kind=CmplxKind), intent(in) :: v(n)
!
   q_lloyd(1:n) = v(1:n)
!
   end subroutine setLloydQ3
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isLloydOn() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   y = LloydOn
!
   end function isLloydOn
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function isLloydOnly() result(y)
!  ===================================================================
   implicit none
!
   logical :: y
!
   y = LloydOnly
!
   end function isLloydOnly
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getLloydQ() result (q)
!  ===================================================================
   implicit none
!
   complex (kind=CmplxKind) :: q(nSpinPola)
!
   q(1:nSpinPola) = q_lloyd(1:nSpinPola)
!
   end function getLloydQ
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine calJost_E0()
!  ===================================================================
   use ErrorHandlerModule, only : ErrorHandler
   use SSSolverModule, only : solveSingleScattering, getJostMatrix
   use MathParamModule, only : pi
!
   implicit none
!
   integer(kind=IntKind) :: is, js, ks, id, kmax, i, kl
   complex(kind=CmplxKind), pointer :: jost(:,:) 
!
   if (.not.Initialized) then
      call ErrorHandler('calJost_E0','KreinModule is not initialized')
   endif
!
   if ( isJostE0_done ) then
      return
   endif
!
   do id = 1, LocalNumAtoms
      kmax=jost_E0(id)%kmax_kkr
      do is = 1, nSpinPola-nSpinCant+1
!        ----------------------------------------
         call solveSingleScattering(is,id,Krein_E0,CZERO)
!        ----------------------------------------
         do js = 1,nSpinCant
            ks = max(is,js)
!           ----------------------------
            jost => getJostMatrix(id,ks)
!           ----------------------------
            do kl = 1,kmax
               jost_E0(id)%jost_at_E0(1:kmax,kl,ks) = jost(1:kmax,kl)
               jost_E0(id)%log_j_E0(kl,ks) = log(jost(kl,kl))
            enddo
         enddo
      enddo 
   enddo
!
   isJostE0_done = .true.
!    
   end subroutine calJost_E0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine print_Krein_Info()
!  ===================================================================
   use ScfDataModule              , only : isLloyd,getLloydMode
   implicit none
   integer :: myKreinMode
!
   myKreinMode = getLloydMode()
!
   if (isLloyd()) then
!
      write(6,'(/,80(''-''))')
      write(6,'(/,24x,a)')'********************************'
      write(6,'( 24x,a )')'*       Output from Krein      *'
      write(6,'(24x,a,/)')'********************************'
      write(6,'(80(''=''))')
!
      select case (myKreinMode)
!
      case(1)
         write (6,'(a)') "Krein is on, Krein mode = 1"
         write (6,'(a)') "This mode is used to run a self consistent calculation using Krein's formula."
         write (6,'(a)') "An additional point gets added to the contour with weight 0, the 2nd last one"
         write (6,'(a)') "to determine the correct charge at the Fermi level."
         write (6,'(a)') "The Green function, DOS, DOS at eF are scaled by a energy independent"
         write (6,'(a)') "correction factor, see output of mufind."

      case DEFAULT
         write (6,'(a)') 'Krein is on, Krein mode is unknown'
      end select
!       
      write(6,'(80(''=''))')      
   endif
!    
   end subroutine print_Krein_Info
! ====================================================================
end module KreinModule
