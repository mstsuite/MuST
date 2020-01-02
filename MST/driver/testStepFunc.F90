program testStepFunc
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : THIRD, ONE, TWO, PI, PI4, CZERO, ZERO, TEN2m8
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use IntegrationModule, only : calIntegration
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics,       &
                                        endSphericalHarmonics
   use BesselModule, only : SphericalBessel
!
   use GauntFactorsModule, only : initGauntFactors, endGauntFactors
!
   use DataServiceCenterModule, only : initDataServiceCenter, &
                                       endDataServiceCenter, &
                                       isDataStorageExisting, &
                                       getDataStorage, RealMark
!
   use MPPModule, only : initMPP, endMPP
!
   use GroupCommModule, only : initGroupComm
!
   use InputModule, only : initInput, endInput, getKeyValue
!
   use SystemModule, only : initSystem, endSystem
   use SystemModule, only : getNumAtoms, getAtomPosition
!
   use ScfDataModule, only : initScfData, endScfData
   use ScfDataModule, only : getPotentialTypeParam
!
   use PotentialTypeModule, only : initPotentialType, endPotentialType, &
                                   isFullPotential
!
   use PolyhedraModule, only : initPolyhedra, endPolyhedra
   use PolyhedraModule, only : genPolyhedron
   use PolyhedraModule, only : getInscrSphRadius
   use PolyhedraModule, only : getOutscrSphRadius
   use PolyhedraModule, only : printPolyhedron
!
   use StepFunctionModule, only : initStepFunction, endStepFunction
   use StepFunctionModule, only : getVolumeIntegration
   use StepFunctionModule, only : printStepFunction
   use StepFunctionModule, only : testStepFunction
   use StepFunctionModule, only : getStepFunction
   use StepFunctionModule, only : getNumGaussRs
   use StepFunctionModule, only : getGaussR
   use StepFunctionModule, only : truncate
   use StepFunctionModule, only : getGaussRWeight
   use StepFunctionModule, only : getNumCriticalRs
   use StepFunctionModule, only : getCriticalR
!
   use IntegerFactorsModule, only : lofj, mofj
!
   implicit none
!
   character (len=10) :: fname
!
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: i, ir, ig, ic, jl, lp, mp
   integer (kind=IntKind) :: nr, ng, nc, ngc, n_int
   integer (kind=IntKind) :: lmax, ngaussr, ngausst, fstatus
   integer (kind=IntKind) :: lmax_max, jmax_max, kmax_max
   integer (kind=IntKind), allocatable :: lmax_step(:), ngr(:), ngt(:)
!
   real (kind=RealKind), pointer :: bravais(:,:)
   real (kind=RealKind), allocatable :: AtomPosition(:,:)
   real (kind=RealKind), allocatable :: r_mesh(:)
   real (kind=RealKind), allocatable :: rgc(:), wgc(:)
   real (kind=RealKind) :: rm, rs, hr, fac, v_mt
!
   real (kind=RealKind), pointer :: rg(:)
   real (kind=RealKind), pointer :: wg(:)
   real (kind=RealKind), pointer :: rc(:)
!
   complex (kind=CmplxKind), allocatable :: func(:,:), g(:)
   complex (kind=CmplxKind), allocatable :: func1(:,:), bjl(:)
   complex (kind=CmplxKind), pointer :: sg(:,:)
   complex (kind=CmplxKind) :: v, vc_mt
   complex (kind=CmplxKind) :: x, kappa
!
!  -------------------------------------------------------------------
   call initMPP()
   call initGroupComm()
   call initDataServiceCenter()
   call initInput()
   call readInputs(def_id,info_id)
   call initScfData(def_id)
   call initPotentialType(getPotentialTypeParam())
   call initSystem(def_id)
!  -------------------------------------------------------------------
!
   if (isDataStorageExisting('Bravais Vector')) then
!     ----------------------------------------------------------------
      bravais => getDataStorage('Bravais Vector',3,3,RealMark)
!     ----------------------------------------------------------------
   else
!     ----------------------------------------------------------------
      call ErrorHandler('testStepFunction','Bravais vector data does not exist')
!     ----------------------------------------------------------------
   endif
!
   NumAtoms = getNumAtoms()
!
   if (NumAtoms < 1) then
      call ErrorHandler('main','invalid NumAtoms',NumAtoms)
   endif
!  -------------------------------------------------------------------
   call initPolyhedra(NumAtoms,bravais,'main',0)
!  -------------------------------------------------------------------
!
   allocate(AtomPosition(1:3,1:NumAtoms))
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
   enddo
   do i=1,NumAtoms
!     ----------------------------------------------------------------
      call genPolyhedron(i,i,NumAtoms,AtomPosition)
!     ----------------------------------------------------------------
      call printPolyhedron(i)
!     ----------------------------------------------------------------
   enddo
!
   allocate( lmax_step(NumAtoms), ngr(NumAtoms), ngt(NumAtoms) )
   fstatus = getKeyValue(info_id,'Default Lmax-Step Func',lmax)
   fstatus = getKeyValue(def_id,'No. Gauss Pts. along r',ngaussr)
   fstatus = getKeyValue(def_id,'No. Gauss Pts. along theta',ngausst)
   do i=1,NumAtoms
      lmax_step(i) = lmax
      ngr(i) = ngaussr
      ngt(i) = ngausst
   enddo
!
   lmax_max = 0
   do i=1,NumAtoms
      lmax_max = max(lmax_max,lmax_step(i))
   enddo
   kmax_max = (lmax_max+1)*(lmax_max+1)
   jmax_max = (lmax_max+1)*(lmax_max+2)/2
   call initSphericalHarmonics(2*lmax_max)
!
!  -------------------------------------------------------------------
   call initGauntFactors(lmax_max, 'none', 0)
!  -------------------------------------------------------------------
   call initStepFunction(NumAtoms,lmax_max,lmax_step,ngr,ngt,'none',0)
!  call initStepFunction(NumAtoms,lmax_step,'none',0)
!  -------------------------------------------------------------------
!
   n_int = 10
   nr = 31+n_int
   allocate( r_mesh(nr), func(nr,jmax_max), func1(nr,jmax_max) )
   allocate( bjl(0:lmax_max) )
   do i=1,NumAtoms
!     ----------------------------------------------------------------
      call printStepFunction(i)
!     ----------------------------------------------------------------
      rs=getOutscrSphRadius(i)
      rm=getInscrSphRadius(i)
!     ----------------------------------------------------------------
      hr= rm/dble(nr-n_int)
      do ir= 1,nr-n_int
         r_mesh(ir)=hr*ir
         func(ir,1)=TWO*sqrt(PI)
      enddo
      hr=(rs-rm)/dble(n_int)
      do ir= 1,n_int
         r_mesh(nr-n_int+ir)= hr*ir + r_mesh(nr-n_int)
         func(nr-n_int+ir,1)=TWO*sqrt(PI)
      enddo
!     ----------------------------------------------------------------
      v = getVolumeIntegration(i,nr,r_mesh,1,0,func,vc_mt,.true.)
!     ----------------------------------------------------------------
      write(6,'(a,2f15.8)')'Volume integration VP = ', v
      write(6,'(a,2f15.8)')'Volume integration MT = ', vc_mt
   enddo
!
   allocate( g(nr) )
   allocate( rgc(500), wgc(500) )
   do i=1,NumAtoms
!     ----------------------------------------------------------------
      rm=getInscrSphRadius(i)
      rs=getOutscrSphRadius(i)
!     ----------------------------------------------------------------
      hr=(rs-rm)/dble(nr-1)
      write(fname,'(a,i1)')'stepfunc',i
      open(unit=10,file=fname,status='unknown',form='formatted')
      do jl = 1,(lmax_step(i)+1)*(lmax_step(i)+2)/2
         do ir=1,nr
            r_mesh(ir)=hr*(ir-1)+rm
!           ----------------------------------------------------------
            func(ir,jl) = getStepFunction(i,jl,r_mesh(ir))
!           ----------------------------------------------------------
         enddo
         if (maxval(abs(func(:,jl))) > TEN2m8) then
            do ir=1,nr
               write(10,'(3i4,f15.8,2x,2f15.8)')ir,lofj(jl),mofj(jl),r_mesh(ir),func(ir,jl)
            enddo
         endif
      enddo
      close(unit=10)
      do ir=1,nr
         func(ir,1)=func(ir,1)*r_mesh(ir)*r_mesh(ir)
      enddo
!     ----------------------------------------------------------------
      call calIntegration(0,nr,r_mesh(1:nr),func(1:nr,1),g)
!     ----------------------------------------------------------------
      write(6,'(a,f15.8)')'Scheme 1 Volume = ',                       &
                          real(g(nr),RealKind)*TWO*sqrt(PI)+PI4*rm*rm*rm*THIRD
!
      ng =  getNumGaussRs(i)
      rg => getGaussR(i)
      wg => getGaussRWeight(i)
      sg => getStepFunction(i)
      v = CZERO
      do ig = 1,ng
         v = v + rg(ig)*rg(ig)*wg(ig)*sg(ig,1)
      enddo
      v = v*TWO*sqrt(PI) + PI4*rm*rm*rm*THIRD
      write(6,'(a,f15.8)')'Scheme 2 Volume = ',real(v,RealKind)
!
      fac = rm*(sqrt(TWO)-ONE)
      nc = getNumCriticalRs(i)
      rc => getCriticalR(i)
      v = CZERO
      do ic = 2, nc
         ngc = ceiling(ngr(i)*(rc(ic)-rc(ic-1))/fac)
!        -------------------------------------------------------------
         call genGaussianPoints(ngc,rc(ic-1),rc(ic),rgc,wgc)
!        -------------------------------------------------------------
         do ig = 1, ngc
            v = v + getStepFunction(i,1,rgc(ig))*wgc(ig)*rgc(ig)**2
         enddo
      enddo
      v = v*TWO*sqrt(PI) + PI4*rm*rm*rm*THIRD
      write(6,'(a,f15.8)')'Scheme 3 Volume = ',real(v,RealKind)
   enddo
!
   do i=1,NumAtoms
!     ----------------------------------------------------------------
      call testStepFunction(i)
!     ----------------------------------------------------------------
   enddo
!
! Test Truncation
!
   kappa = sqrt(cmplx(.30d0,.20d0,kind=CmplxKind))
   do i=1,NumAtoms
      rs=getOutscrSphRadius(i)
      rm=getInscrSphRadius(i)
!     ----------------------------------------------------------------
      hr= rm/dble(nr-n_int)
      do ir= 1,nr-n_int
         r_mesh(ir)=hr*ir
      enddo
      hr=(rs-rm)/dble(n_int)
      do ir= 1,n_int
         r_mesh(nr-n_int+ir)= hr*ir + r_mesh(nr-n_int)
      enddo
      do ir=1,nr
         x=kappa*r_mesh(ir)
!        -------------------------------------------------------------
         call SphericalBessel(lmax_max,x,bjl)
!        -------------------------------------------------------------
         do lp=0,lmax_max
            do mp=0,lp
               jl = (lp+1)*(lp+2)/2-lp+mp
               func(ir,jl)=x*bjl(lp)
            enddo
         enddo
      enddo
!
      call truncate(i, nr-n_int, nr, r_mesh, func, jmax_max, n_int+1, func1, jmax_max)
!
!     ----------------------------------------------------------------
      v = getVolumeIntegration(i,nr,r_mesh,kmax_max,jmax_max,0,func,v_mt)
!     ----------------------------------------------------------------
      write(6,'(a,2f15.8)')'Func Volume integration VP = ', v
      write(6,'(a,2f15.8)')'Func Volume integration MT = ', v_mt
      do ir=1,n_int+1
         do lp=0,lmax_max
            do mp=0,lp
               jl = (lp+1)*(lp+2)/2-lp+mp
               func(nr-n_int-1+ir,jl)= func1(ir,jl)
            enddo
         enddo
      enddo
!     ----------------------------------------------------------------
      v = getVolumeIntegration(i,nr,r_mesh,kmax_max,jmax_max,0,func,v_mt,.false.)
!     ----------------------------------------------------------------
      write(6,'(a,2f15.8)')'Func1 Volume integration VP = ', v
      write(6,'(a,2f15.8)')'Func1 Volume integration MT = ', v_mt
   enddo
!
   deallocate(AtomPosition, r_mesh, func, g)
   deallocate(lmax_step, ngr, ngt)
   deallocate(rgc, wgc)
   nullify ( rg, wg, sg )
!
!  -------------------------------------------------------------------
   call endSphericalHarmonics()
   call endStepFunction()
   call endGauntFactors()
   call endPolyhedra()
   call endMPP()
!  -------------------------------------------------------------------
   stop 'Ok'
end program testStepFunc
