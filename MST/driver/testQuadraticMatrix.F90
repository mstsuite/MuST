program testQuadraticMatrix
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   use QuadraticMatrixModule, only : initQuadraticMatrix, endQuadraticMatrix, &
                                     solveQuadraticEquation, getEigenValue
   use MatrixDeterminantModule, only : MtxDet
!
   implicit none
!
   integer (kind=IntKind), parameter :: ndim = 3
   integer (kind=IntKind), parameter :: np = 51
   integer (kind=IntKind) :: i, j, k, info, nv
!
   real (kind=RealKind), parameter :: a = -1.0d0
   real (kind=RealKind), parameter :: b =  1.0d0
   real (kind=RealKind) :: e, de
!
   complex (kind=CmplxKind) :: s0(ndim,ndim), s1(ndim,ndim), s2(ndim,ndim)
   complex (kind=CmplxKind) :: s(ndim,ndim), det
   complex (kind=CmplxKind), pointer :: pv(:)
!
   call initQuadraticMatrix(ndim)
!
   do j = 1, ndim
      do i = 1, ndim
         s0(i,j) =  5.0*i-2.0*j
         s1(i,j) = -2.89+0.33*i-2.03*j
         s2(i,j) =  0.15*s1(i,j)*s1(i,j)   !  0.52-0.01*i*i+0.05*j*j
      enddo
   enddo
!
   call solveQuadraticEquation(s0,s1,s2,info)
!
   if (info /= 0) then
      stop 'Error in s0, s1, s2'
   endif
!
   pv => getEigenValue(nv)
   write(6,'(a)')'Index           Eigenvalues                Determinant'
   do k = 1, ndim*2
      do j = 1, ndim
         do i = 1, ndim
            s(i,j) = s0(i,j)+s1(i,j)*pv(k)+s2(i,j)*pv(k)*pv(k)
         enddo
      enddo
      call MtxDet(ndim,s,det)
      write(6,'(i4,4x,2d15.6,8x,2d15.6)')k,pv(k),det
   enddo
   nullify(pv)
!
   write(6,'(//,a)')'Index   Delta Value                Determinant'
   de = (b-a)/real(np-1,kind=RealKind)
   do k = 1, np
      e = de*(k-1) + a
      do j = 1, ndim
         do i = 1, ndim
            s(i,j) = s0(i,j)+s1(i,j)*e+s2(i,j)*e*e
         enddo
      enddo
      call MtxDet(ndim,s,det)
      write(6,'(i4,4x,d15.6,8x,2d15.6)')k,e,det
   enddo
!
   call endQuadraticMatrix()
!
end program testQuadraticMatrix
