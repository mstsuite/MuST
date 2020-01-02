!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  module SurfJNYModule
! ====================================================================
  use KindParamModule, only : IntKind
  use KindParamModule, only : RealKind
  use KindParamModule, only : CmplxKind
!
  use DefinedTypeModule, only : AtomInfoStruc
!
  private :: IntKind, RealKind, CmplxKind, AtomInfoStruc
!
  integer (kind=1), pointer :: surface_index(:)
  integer (kind=IntKind) :: NumSurfPoints
  integer (kind=IntKind), pointer :: JHunt(:)
  integer (kind=IntKind) :: NumSurfR
  integer (kind=IntKind), pointer :: I2JSurf(:)
!
  real (kind=RealKind), pointer :: n_dot_r(:)
  real (kind=RealKind), pointer :: r_surf(:)
  real (kind=RealKind), pointer :: weight(:)
  real (kind=RealKind), pointer :: SurfR(:)
!
  complex (kind=CmplxKind), pointer :: ylm(:,:)
  complex (kind=CmplxKind), pointer :: nDotGradY(:,:)
  complex (kind=CmplxKind), pointer :: sjl(:,:)
  complex (kind=CmplxKind), pointer :: dsjl(:,:)
  complex (kind=CmplxKind), pointer :: snl(:,:)
  complex (kind=CmplxKind), pointer :: dsnl(:,:)
! complex (kind=CmplxKind), allocatable :: scfac(:)
!
  type (AtomInfoStruc), private, pointer :: ThisAtom
!
  private :: TransferPointer
!
  contains
!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine InitSurfJNY()
! ====================================================================
  use MathParamModule, only : ten2m6, ten2m8, czero
!
  use SortModule, only : HeapSort
!
  use ErrorHandlerModule, only : ErrorHandler
!
  use AtomInfoModule, only : CurrentAtom
!
  use SurfElemModule, only : num_surf_pts
  use SurfElemModule, only : first_surface_point
  use SurfElemModule, only : PointStruc
  use SurfElemModule, only : surf_norm
  use SurfElemModule, only : EndSurfElem
! 
  use SphericalHarmonicsModule, only : calYlm
! 
  use ArraySizeModule, only : lmax_kkr, lmax_phi, kmax_kkr, kmax_phi
! 
  implicit none
!
  integer (kind=IntKind) :: i,j,k,lmax,kmax
  integer (kind=IntKind), allocatable :: indx(:)
!
  real (kind=RealKind) :: r
  real (kind=RealKind) :: x,y,z
  real (kind=RealKind), allocatable :: rtmp(:)
!
  complex (kind=CmplxKind), allocatable :: grady(:,:)
  complex (kind=CmplxKind), pointer :: pylm(:)
!
  type (PointStruc), pointer :: current_point
!
  CurrentAtom%SurfJNY%NumSurfPoints=num_surf_pts
  NumSurfPoints=CurrentAtom%SurfJNY%NumSurfPoints
  ThisAtom=>CurrentAtom
  CurrentAtom%SurfJNY%kappa=czero
  lmax=max(lmax_kkr, lmax_phi)
  kmax=(lmax+1)*(lmax+1)
! --------------------------------------------------------------------
  allocate( CurrentAtom%SurfJNY%n_dot_r(1:NumSurfPoints),             &
            CurrentAtom%SurfJNY%r_surf(1:NumSurfPoints),              &
            CurrentAtom%SurfJNY%surface_index(1:NumSurfPoints),       &
            CurrentAtom%SurfJNY%weight(1:NumSurfPoints) )
  allocate( grady(1:kmax,3),                                          &
            CurrentAtom%SurfJNY%ylm(1:kmax,1:NumSurfPoints),          &
            CurrentAtom%SurfJNY%nDotGradY(1:kmax,1:NumSurfPoints) )
! allocate(scfac(0:lmax))
! --------------------------------------------------------------------
  surface_index => CurrentAtom%SurfJNY%surface_index
  nDotGradY => CurrentAtom%SurfJNY%nDotGradY
  n_dot_r => CurrentAtom%SurfJNY%n_dot_r
  r_surf => CurrentAtom%SurfJNY%r_surf
  weight => CurrentAtom%SurfJNY%weight
  ylm => CurrentAtom%SurfJNY%ylm
!
! ====================================================================
! loop over the points on the boundary planes that defines the 
! Voronoi Polyhedron.
! ====================================================================
  current_point => first_surface_point
  do i=1,NumSurfPoints
     x=current_point%position(1) 
     y=current_point%position(2) 
     z=current_point%position(3) 
     r=sqrt(x*x+y*y+z*z)
     if (r.lt.ten2m8) then
!       --------------------------------------------------------------
        call ErrorHandler('IntSurfJNY','r is too small',r)
!       --------------------------------------------------------------
     endif
     r_surf(i)=r
!
!    =================================================================
!    calculating Ylm => ylm, Grad{Ylm} => grady.
!    =================================================================
     pylm=>ylm(1:kmax,i)
!    -----------------------------------------------------------------
     call calYlm(x,y,z,lmax,pylm,grady)
!    -----------------------------------------------------------------
!
!    =================================================================
!                ->  ->
!    calculating n * Grad{Ylm} => nDotGradY
!    =================================================================
     weight(i)=current_point%weight
     surface_index(i)=current_point%surface_index
     k=surface_index(i)
     nDotGradY(1:kmax,i)=surf_norm(1,k)*grady(1:kmax,1)+              &
                         surf_norm(2,k)*grady(1:kmax,2)+              &
                         surf_norm(3,k)*grady(1:kmax,3)
     n_dot_r(i)=(surf_norm(1,k)*x+surf_norm(2,k)*y+surf_norm(3,k)*z)/r
!
!    =================================================================
!    update current_point and current_ylm
!    =================================================================
     current_point => current_point%next_point
  enddo
!
! --------------------------------------------------------------------
  call EndSurfElem()
! --------------------------------------------------------------------
  deallocate( grady )
! --------------------------------------------------------------------
!
! ====================================================================
! setup I2JSurf and SurfR. They are used for surface integration
! NumSurfR: number of distiguished r_surf
! SurfR:    stores distiguished r_surf
! I2JSurf:  mapping indeces for r_surf -> SurfR
! --------------------------------------------------------------------
  allocate(rtmp(1:NumSurfPoints), indx(1:NumSurfPoints))
! --------------------------------------------------------------------
  rtmp(1:NumSurfPoints)=r_surf(1:NumSurfPoints)
  call HeapSort(NumSurfPoints,rtmp,indx)
  CurrentAtom%SurfJNY%NumSurfR=1
  r=rtmp(1)
  do i=2,NumSurfPoints
     if (rtmp(i)-r .gt. ten2m6) then
         CurrentAtom%SurfJNY%NumSurfR=CurrentAtom%SurfJNY%NumSurfR+1
         r=rtmp(i)
     endif
  enddo
  NumSurfR=CurrentAtom%SurfJNY%NumSurfR
! --------------------------------------------------------------------
  allocate( CurrentAtom%SurfJNY%JHunt(1:NumSurfR),                    &
            CurrentAtom%SurfJNY%SurfR(1:NumSurfR),                    &
            CurrentAtom%SurfJNY%I2JSurf(1:NumSurfPoints) )
! --------------------------------------------------------------------
  CurrentAtom%SurfJNY%SurfR(1)=rtmp(1)
  k=indx(1)
  CurrentAtom%SurfJNY%I2JSurf(k)=1
  j=1
  do i=2,NumSurfPoints
     if (rtmp(i)-CurrentAtom%SurfJNY%SurfR(j) .gt. ten2m6) then
         j=j+1
         CurrentAtom%SurfJNY%SurfR(j)=rtmp(i)
     endif
     k=indx(i)
     CurrentAtom%SurfJNY%I2JSurf(k)=j
  enddo
  do j=1,NumSurfR
     CurrentAtom%SurfJNY%JHunt(j)=CurrentAtom%jmt
!    -----------------------------------------------------------------
     call hunt(CurrentAtom%r_mesh,CurrentAtom%jend_plus_n,            &
               CurrentAtom%SurfJNY%SurfR(j),CurrentAtom%SurfJNY%JHunt(j))
!    -----------------------------------------------------------------
  enddo
  JHunt=>CurrentAtom%SurfJNY%JHunt
  SurfR=>CurrentAtom%SurfJNY%SurfR
  I2JSurf=>CurrentAtom%SurfJNY%I2JSurf
! --------------------------------------------------------------------
  deallocate( rtmp, indx )
! --------------------------------------------------------------------
!
! --------------------------------------------------------------------
  allocate( CurrentAtom%SurfJNY%sjl(0:lmax_phi,1:NumSurfPoints),      &
            CurrentAtom%SurfJNY%snl(0:lmax_phi,1:NumSurfPoints),      &
            CurrentAtom%SurfJNY%dsjl(0:lmax_phi,1:NumSurfPoints),     &
            CurrentAtom%SurfJNY%dsnl(0:lmax_phi,1:NumSurfPoints) )
! --------------------------------------------------------------------
  sjl => CurrentAtom%SurfJNY%sjl
  dsjl => CurrentAtom%SurfJNY%dsjl
  snl => CurrentAtom%SurfJNY%snl
  dsnl => CurrentAtom%SurfJNY%dsnl
!
  end subroutine InitSurfJNY
!
!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine EndSurfJNY()
! ====================================================================
!
  use AtomInfoModule, only : CurrentAtom
!
  implicit none
!
  if (.not.associated(ThisAtom,CurrentAtom)) then
!    -----------------------------------------------------------------
     call TransferPointer(CurrentAtom%SurfJNY)
!    -----------------------------------------------------------------
     ThisAtom=>CurrentAtom
  endif
!
! --------------------------------------------------------------------
  deallocate( surface_index, weight, n_dot_r, r_surf, ylm, nDotGradY, &
              sjl, snl, dsjl, dsnl, JHunt, SurfR, I2JSurf )
! deallocate( scfac )
! --------------------------------------------------------------------
!
  end subroutine EndSurfJNY
!
!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine TransferPointer(SurfJNY)
! ====================================================================
!
  use DefinedTypeModule, only : SurfJNYStruc
!
  implicit none
!
  type (SurfJNYStruc) :: SurfJNY
!
  NumSurfPoints=SurfJNY%NumSurfPoints
  NumSurfR=SurfJNY%NumSurfR
!
  surface_index => SurfJNY%surface_index
!
  n_dot_r => SurfJNY%n_dot_r
  r_surf => SurfJNY%r_surf
  JHunt => SurfJNY%JHunt
  weight => SurfJNY%weight
  SurfR => SurfJNY%SurfR
  I2JSurf => SurfJNY%I2JSurf
!
  ylm => SurfJNY%ylm
  nDotGradY => SurfJNY%nDotGradY
  sjl => SurfJNY%sjl
  dsjl => SurfJNY%dsjl
  snl => SurfJNY%snl
  dsnl => SurfJNY%dsnl
!
  end subroutine TransferPointer
!
!
! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine UpdateSurfJNY(Energy)
! ====================================================================
!
! ********************************************************************
! Must call InitSurfJNY before this routine is called.
! It update spherical Bessel function on the atomic cell surface for a 
! new Energy, which is the square of kappa.
! ********************************************************************
!
  use MathParamModule, only : ten2m8
  use MathParamModule, only : cone
!
  use AtomInfoModule, only : CurrentAtom
  use AtomInfoModule, only : c2inv
!
  use BesselModule, only : SphericalBessel, SphericalNeumann
!
  use ArraySizeModule, only : lmax_phi, kmax_phi, lofk, mofk, m1m
!
  implicit none
!
  integer (kind=IntKind) :: i, l, m, kl, klc, mfac
!
  complex (kind=CmplxKind), intent(in) :: Energy
!
  complex (kind=CmplxKind) :: kr, kappa
  complex (kind=CmplxKind), pointer :: pb(:), pdb(:)
!
  kappa = sqrt(Energy+Energy*Energy*c2inv)
!
  if (.not.associated(ThisAtom,CurrentAtom)) then
!    -----------------------------------------------------------------
     call TransferPointer(CurrentAtom%SurfJNY)
!    -----------------------------------------------------------------
     ThisAtom=>CurrentAtom
  endif
  if (abs(CurrentAtom%SurfJNY%kappa-kappa).lt.ten2m8) then
     return
  else
     CurrentAtom%SurfJNY%kappa=kappa
     do i=1,NumSurfPoints
        kr=kappa*r_surf(i)
        pb => sjl(0:lmax_phi,i)
        pdb => dsjl(0:lmax_phi,i)
!       --------------------------------------------------------------
        call SphericalBessel(lmax_phi,kr,pb,pdb)
!       --------------------------------------------------------------
        dsjl(0:lmax_phi,i)=kappa*dsjl(0:lmax_phi,i)
        pb => snl(0:lmax_phi,i)
        pdb => dsnl(0:lmax_phi,i)
!       --------------------------------------------------------------
        call SphericalNeumann(lmax_phi,kr,pb,pdb)
!       --------------------------------------------------------------
        dsnl(0:lmax_phi,i)=kappa*dsnl(0:lmax_phi,i)
     enddo
!    do l=0,lmax_phi
!       scfac(l)=cone/dsnl(l,NumSurfR)
!    enddo
  endif
!
  end subroutine UpdateSurfJNY
!
!
  end module SurfJNYModule
