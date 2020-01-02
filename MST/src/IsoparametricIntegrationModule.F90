module IsoparametricIntegrationModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
   use MathParamModule, only : ZERO, THIRD, HALF, ONE, TWO, THREE
   use MathParamModule, only : PI, PI2, PI4
   use MathParamModule, only : TEN2m3, TEN2m4, TEN2m6, TEN2m8, TEN2m10, TEN2m12, TEN
   use MathParamModule, only : CZERO, CONE, SQRTm1
   use IntegerFactorsModule, only : lofj, mofj, jofk, lofk, mofk
!
public :: initIsoparametricIntegration,&
          endIsoparametricIntegration, &
          getVolumeIntegral
!
   interface getVolumeIntegral
      module procedure getVolumeIntegral_0
      module procedure getVolumeIntegral_jl, getVolumeIntegral_kl
      module procedure getVolumeIntegral_3kl
   end interface
!
private
   include "imp_inp.h"
!
   logical :: Initialized = .false.
!
   integer(kind=IntKind) :: LocalNumAtoms
   integer(kind=IntKind) :: Nqp
   integer(kind=IntKind) :: kmax_max, jmax_max, lmax_max
!
   integer (kind=IntKind), parameter :: nrint = 5 
!
   integer (kind=IntKind), allocatable:: nVertex1(:),nVertex2(:), &
                                         Ivdex2(:,:,:), nfaces(:),&
                                         Iface(:,:,:),fcount(:), iType(:)
   real (kind=RealKind) :: lattice(3), t(3,3)
   real (kind=RealKind), allocatable :: Rmt(:), weight(:), r(:), &
                                        rmag(:,:,:,:,:,:),         &
                                        vj(:,:,:,:,:),            &
                                        vertices(:,:,:)
   complex(kind=CmplxKind), allocatable, target :: Ylm_rmag(:,:,:,:,:,:)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initIsoparametricIntegration(lmax_in)
!  ===================================================================
   use IntegerFactorsModule, only : initIntegerFactors
   use Atom2ProcModule, only : getGlobalIndex, getLocalNumAtoms
!
   use PolyhedraModule, only : getNumPlanes, getCorner
   use PolyhedraModule, only : printPolyhedron, getNumCorners
   use PolyhedraModule, only : printPolyhedraTable, getNumPlaneCorners
   use PolyhedraModule, only : getOutscrSphRadius, getIndexCornerPlane, &
                               getInscrSphRadius, getWignerSeitzRadius
!
   use SphericalHarmonicsModule, only : calYlm
!
   implicit none
!
   integer(kind=IntKind) :: lmax_in
!
   integer(kind=IntKind) :: i, m, n, mn, ig, i0, j0, k0, iVP
   integer(kind=IntKind) :: jl,kl,l
!
   real(kind=RealKind) :: posi(3), rr
   real(kind=RealKind), allocatable :: VPInt(:)
!
   real(kind=RealKind), pointer :: p2r(:,:)
!
   complex(kind=CmplxKind), pointer :: ylm(:)
!
   if (Initialized) then
      return
   endif
!
   LocalNumAtoms = getLocalNumAtoms()
   lmax_max = lmax_in
   kmax_max = (lmax_in+1)**2
   jmax_max = (lmax_in+1)*(lmax_in+2)/2
!  -------------------------------------------------------------------
   call initIntegerFactors(lmax_max)
!  -------------------------------------------------------------------
!
   allocate( r(LocalNumAtoms),rmt(LocalNumAtoms),iType(LocalNumAtoms))
   allocate( nfaces(LocalNumAtoms))
   allocate( nVertex1(LocalNumAtoms), nVertex2(LocalNumAtoms) )
   allocate( Ivdex2(LocalNumAtoms,MNV,2), vertices(LocalNumAtoms,3,MNV),&
             Iface(LocalNumAtoms,MNF,MNV+1), fcount(LocalNumAtoms) )
!
   allocate( weight(MNqp), rmag(3,MNqp,MNqp,MNqp,MNF,LocalNumAtoms) )
   allocate( vj(MNqp,MNqp,MNqp,MNF,LocalNumAtoms) )
!
   allocate(VPInt(LocalNumAtoms))
!
   lattice = ONE
   t       = ONE
   Iface   = 0
   Ivdex2  = 0
   fcount  = 0
!
   Nqp = 8
   do i = 1,LocalNumAtoms
!
      r(i) = getInscrSphRadius(i) ! MT radius
      rmt(i) = r(i)! MT radius
      iType(i) = i
      nfaces(i)   = getNumPlanes(i)
      nVertex1(i) = getNumCorners(i)
      nVertex2(i) = nVertex1(i)
      p2r => getCorner(i)
!
      do n = 1,nVertex1(i)
         vertices(i,1,n) = p2r(1,n)
         vertices(i,2,n) = p2r(2,n)
         vertices(i,3,n) = p2r(3,n)
         Ivdex2(i,n,1)  = n
         Ivdex2(i,n,2)  = n
      enddo
      do m = 1,nfaces(i)
         mn = getNumPlaneCorners(i,m)
         Iface(i,m,1) = mn 
         do n = 1,mn
            Iface(i,m,1+n) = getIndexCornerPlane(i,m,n)
!           write(6,'(a,3i5,3(1x,f18.10))') "Iface::",m,n+1,Iface(i,m,1+n), p2r(1:3,Iface(i,m,1+n))
         enddo
      enddo
!
   enddo
!
   iVP = 0
   call VPI_main( iVP, Nqp, LocalNumAtoms, lattice, r, rmt,           &
                  LocalNumAtoms, iType, t, nVertex1, nVertex2,        &
                  Ivdex2, nfaces, Iface, vertices, fcount,            &
                  weight, rmag, vj )
!
   allocate(Ylm_rmag(kmax_max,MNqp,MNqp,MNqp,MNF,LocalNumAtoms))
!
   VPInt = 0.0D0
   do i =1,LocalNumAtoms
      do n = 1,fcount(i)
         do i0 = 1,Nqp
         do j0 = 1,Nqp
         do k0 = 1,Nqp 
            ylm => Ylm_rmag(1:kmax_max,i0,j0,k0,n,i)
            posi = rmag(1:3,i0,j0,k0,n,i)
            rr = sqrt(posi(1)**2+posi(2)**2+posi(3)**2)
            if ( rr<0.0000001d0) then
               ylm=CZERO
            else
               call calYlm(posi,lmax_max,ylm)
               VPInt(i) = VPInt(i) + weight(i0)*weight(j0)*weight(k0)* &
                                     vJ(i0,j0,k0,n,i)
            endif
         enddo
         enddo
         enddo
      enddo
      write(6,'(60(''=''))')
      write(6,'(3x,a)') "*  Isorarametric - Volumes  *"
      write(6,'(60(''-''))')
      write(6,'(a,1x,i4,2(1x,a,1x,f18.10))') "i=",i,"rmt =",rmt(i),"VPInt =", VPInt(i)
      write(6,'(60(''=''))')
   enddo
!
   deallocate( VPInt )
   deallocate( r, iType, nfaces, nVertex1, nVertex2 )
   deallocate( Ivdex2, vertices, Iface )
!
   Initialized = .true.
!
   end subroutine initIsoparametricIntegration
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endIsoparametricIntegration()
!  ===================================================================
   use IntegerFactorsModule, only : endIntegerFactors
   implicit none
!
   if (.not.Initialized) then
      return
   endif
!
   deallocate( weight, rmag, vj, Rmt, fcount, Ylm_rmag )
!  -------------------------------------------------------------------
   call endIntegerFactors()
!  -------------------------------------------------------------------
!
   Initialized = .false.
!
   end subroutine endIsoparametricIntegration
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVolumeIntegral_0(id,nr,rr,fr)       result(vint)
!  ===================================================================
   use InterpolationModule, only : PolyInterp
!
   implicit none
!
   complex(kind=CmplxKind) :: vint, v0
!
   integer(kind=IntKind), intent(in) :: id, nr
!
   real(kind=RealKind), target :: rr(1:nr)
!
   complex(kind=CmplxKind), target :: fr(1:nr)
!
   integer(kind=IntKind)   :: n, i0, j0, k0, l, kl, jl, ir, irp
   real(kind=RealKind)     :: r, err, posi(3)
   complex(kind=CmplxKind) :: fr_in
   complex(kind=CmplxKind), pointer :: ylm(:)
!
   vint = CZERO
!
   do n = 1,fcount(id)
      do i0 = 1,Nqp
      do j0 = 1,Nqp
      do k0 = 1,Nqp 
         posi = rmag(1:3,i0,j0,k0,n,id)
         r = sqrt(posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
!        -------------------------------------------------------------
         call hunt(nr,rr,r,ir)
!        -------------------------------------------------------------
         if (ir > nr-(nrint-1)/2) then
            irp=nr-nrint+1
         else if (2*ir+1 > nrint) then
            irp=ir-(nrint-1)/2
         else
            irp=1
         endif
!        -------------------------------------------------------------
         call PolyInterp( nrint, rr(irp:irp+nrint-1),                 &
                          fr(irp:irp+nrint-1), r, fr_in, err)
!        -------------------------------------------------------------
         if ( mod(i0,8)==0 .and. mod(j0,8)==0 .and. mod(k0,8)==0) then
            print*,"Iso:: r = ",r," fr_in =", fr_in
         endif
!
         vint = vint + weight(i0)*weight(j0)*weight(k0)*   &
                       vj(i0,j0,k0,n,id)*fr_in
      enddo
      enddo
      enddo
   enddo
!
   end function getVolumeIntegral_0
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVolumeIntegral_jl(id,nr,rr,jmax,frj)       result(vint)
!  ===================================================================
   use InterpolationModule, only : PolyInterp
!
   implicit none
!
   complex(kind=CmplxKind) :: vint, v0
!
   integer(kind=IntKind), intent(in) :: id, nr, jmax
!
   real(kind=RealKind), intent(in) :: rr(1:nr)
!
   complex(kind=CmplxKind), target :: frj(1:nr,1:jmax)
!
   integer(kind=IntKind)   :: n, i0, j0, k0, l, kl, jl, ir, irp, kmax
   real(kind=RealKind)     :: r, err, posi(3)
   complex(kind=CmplxKind) :: frj_in
   complex(kind=CmplxKind), pointer :: ylm(:)
!
   vint = CZERO
   kmax = (lofj(jmax)+1)**2
!
   do n = 1,fcount(id)
      do i0 = 1,Nqp
      do j0 = 1,Nqp
      do k0 = 1,Nqp 
         posi = rmag(1:3,i0,j0,k0,n,id)
         r = sqrt(posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
!        -------------------------------------------------------------
         call hunt(nr,rr,r,ir)
!        -------------------------------------------------------------
         if (ir > nr-(nrint-1)/2) then
            irp=nr-nrint+1
         else if (2*ir+1 > nrint) then
            irp=ir-(nrint-1)/2
         else
            irp=1
         endif
!
         ylm => Ylm_rmag(:,i0,j0,k0,n,id)
         v0 = CZERO
         do jl = 1,jmax
            l = lofj(jl)
!           ----------------------------------------------------------
            call PolyInterp( nrint, rr(irp:irp+nrint-1),              &
                             frj(irp:irp+nrint-1,jl), r, frj_in, err)
!           ----------------------------------------------------------
            kl = (l+1)*(l+1)-l+mofj(jl)
            if (mofj(jl) == 0) then
               v0 = v0 + frj_in*ylm(kl)
            else
               v0 = v0 + TWO*frj_in*ylm(kl)
            endif
         enddo
!
         vint = vint + weight(i0)*weight(j0)*weight(k0)*   &
                       vj(i0,j0,k0,n,id)*v0
      enddo
      enddo
      enddo
   enddo
!
   end function getVolumeIntegral_jl
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVolumeIntegral_kl(id,nr,rr,jmax,kmax,frk)  result(vint)
!  ===================================================================
   use InterpolationModule, only : PolyInterp
!
   implicit none
!
   complex(kind=CmplxKind) :: vint, v0
!
   integer(kind=IntKind), intent(in) :: id, nr, jmax, kmax
!
   real(kind=RealKind), intent(in) :: rr(1:nr)
!
   complex(kind=CmplxKind), target:: frk(1:nr,1:kmax)
!
   integer(kind=IntKind)   :: n, i0, j0, k0, l, kl, jl, ir, irp
   real(kind=RealKind)     :: r, err, posi(3)
   complex(kind=CmplxKind) :: frk_in
   complex(kind=CmplxKind), pointer :: ylm(:)
!
   vint = CZERO
!
   do n = 1,fcount(id)
      do i0 = 1,Nqp
      do j0 = 1,Nqp
      do k0 = 1,Nqp 
         posi = rmag(1:3,i0,j0,k0,n,id)
         r = sqrt(posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
!        -------------------------------------------------------------
         call hunt(nr,rr,r,ir)
!        -------------------------------------------------------------
         if (ir > nr-(nrint-1)/2) then
            irp=nr-nrint+1
         else if (2*ir+1 > nrint) then
            irp=ir-(nrint-1)/2
         else
            irp=1
         endif
!
         ylm => Ylm_rmag(:,i0,j0,k0,n,id)
         v0 = CZERO
         do kl = 1,kmax
            l = lofk(kl)
!           ----------------------------------------------------------
            call PolyInterp( nrint, rr(irp:irp+nrint-1),              &
                             frk(irp:irp+nrint-1,kl), r, frk_in, err)
!           ----------------------------------------------------------
            v0 = v0 + frk_in*ylm(kl)
         enddo
!
         vint = vint + weight(i0)*weight(j0)*weight(k0)*   &
                       vj(i0,j0,k0,n,id)*v0
      enddo
      enddo
      enddo
   enddo
!
   end function getVolumeIntegral_kl
!  ===================================================================
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getVolumeIntegral_3kl( id, nr, rr, k1, fr1, j2, k2, fr2,&
                                   ffr2, j3, k3, fr3, ffr3) result(vint)
!  ===================================================================
   use InterpolationModule, only : PolyInterp
   use GauntFactorsModule, only : getK3
   use GauntFactorsModule, only : getNumK3
!
   implicit none
!
   complex(kind=CmplxKind) :: vint, v02, v0
!
   integer(kind=IntKind), intent(in) :: id, nr, k1, j2, k2, j3, k3
!
   real(kind=RealKind), intent(in) :: rr(1:nr)
!
   integer(kind=IntKind) :: ffr2(1:k2)
   integer(kind=IntKind) :: ffr3(1:k3)
!
   complex(kind=CmplxKind) :: fr1(1:nr)
   complex(kind=CmplxKind) :: fr2(1:nr,1:k2)
   complex(kind=CmplxKind) :: fr3(1:nr,1:k3)
!
   integer(kind=IntKind)   :: n, i0, j0, k0, l, kl, jl, ir, irp, m1m
   integer(kind=IntKind)   :: kl2, kl3, i2, nj2_i2
   integer (kind=IntKind), pointer :: kj2_i2(:)
   integer (kind=IntKind), pointer :: nj3(:,:)
   integer (kind=IntKind), pointer :: kj3(:,:,:)
   real(kind=RealKind)     :: r, err, posi(3)
   complex(kind=CmplxKind) :: f_in
   complex(kind=CmplxKind), pointer :: ylm(:)
   complex(kind=CmplxKind), allocatable :: f_tmp(:)
!
   allocate(f_tmp(1:nr))
!
   vint = CZERO
!
   nj3 => getNumK3()
   kj3 => getK3()
!
   m1m= (-1)**mofk(k1)
   Loop_kl3: do kl3 = k3,1,-1
      if ( ffr3(kl3) == 0 ) then
         cycle Loop_kl3
      endif
      v02 = CZERO
   Loop_kl2: do kl2 = k2,1,-1
      if ( ffr2(kl2) == 0 ) then
         cycle Loop_kl2
      endif
      do ir = 1,nr
         f_tmp(ir) =  fr1(ir)*fr2(ir,kl2)*fr3(ir,kl3)
      enddo
      v0 = CZERO
      do n = 1,fcount(id)
         do i0 = Nqp,1,-1
         do j0 = Nqp,1,-1
         Loop_k0: do k0 = Nqp,1,-1
            ylm => Ylm_rmag(:,i0,j0,k0,n,id)
            posi = rmag(1:3,i0,j0,k0,n,id)
            r = sqrt(posi(1)*posi(1)+posi(2)*posi(2)+posi(3)*posi(3))
            if (r<0.0000001d0 ) then
              cycle Loop_k0
            endif
!           ----------------------------------------------------------
            call hunt(nr,rr,r,ir)
!           ----------------------------------------------------------
            if (ir > nr-(nrint-1)/2) then
               irp=nr-nrint+1
            else if (2*ir+1 > nrint) then
               irp=ir-(nrint-1)/2
            else
               irp=1
            endif
!           ----------------------------------------------------------
            call PolyInterp( nrint, rr(irp:irp+nrint-1),              &
                             f_tmp(irp:irp+nrint-1), r, f_in, err)
!           ----------------------------------------------------------
            v0 = v0 + weight(i0)*weight(j0)*weight(k0)*vj(i0,j0,k0,n,id)* &
                      f_in*m1m*conjg(ylm(k1))*ylm(kl2)*ylm(kl3)
         enddo Loop_k0
         enddo
         enddo
      enddo
      v02 = v02 + v0
   enddo Loop_kl2
      vint = vint + v02
   enddo Loop_kl3
!
   deallocate(f_tmp)
!
   end function getVolumeIntegral_3kl
!  ===================================================================
end module IsoparametricIntegrationModule 
