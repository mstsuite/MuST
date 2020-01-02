program testStepFunc2
   use KindParamModule, only : IntKind, RealKind, CmplxKind
!
   use MathParamModule, only : THIRD, ONE, TWO, PI, PI4, CZERO, ZERO
!
   use ErrorHandlerModule, only : ErrorHandler
!
   use IntegrationModule, only : calIntegration
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
   use GroupCommModule, only : initGroupComm, getGroupID
!
   use ProcMappingModule, only : initProcMapping, createParallelization, endProcMapping
!
   use Atom2ProcModule, only : initAtom2Proc, endAtom2Proc
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
   use PolyhedraModule, only : printPolyhedron, printPolyhedraTable
!
   use StepFunctionModule, only : initStepFunction, endStepFunction
   use StepFunctionModule, only : getVolumeIntegration
   use StepFunctionModule, only : printStepFunction
   use StepFunctionModule, only : testStepFunction
   use StepFunctionModule, only : getStepFunction
   use StepFunctionModule, only : getNumGaussRs
   use StepFunctionModule, only : getGaussR
   use StepFunctionModule, only : getGaussRWeight
   use StepFunctionModule, only : getNumCriticalRs
   use StepFunctionModule, only : getCriticalR
!
   use SphericalHarmonicsModule, only : initSphericalHarmonics
!
   implicit none
!
   character (len=12) :: fname
!
   integer (kind=IntKind) :: def_id, info_id
   integer (kind=IntKind) :: NumAtoms
   integer (kind=IntKind) :: i, ir, ig, ic, jl, l, m, kl
   integer (kind=IntKind) :: nr, ng, nc, ngc
   integer (kind=IntKind) :: lmax, ngaussr, ngausst, fstatus
   integer (kind=IntKind) :: lmax_max, jmax_max, kmax_max, kmax_step
   integer (kind=IntKind), allocatable :: lmax_step(:), ngr(:), ngt(:)
   integer (kind=IntKind), parameter :: NA_limit = 10
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
   complex (kind=CmplxKind), allocatable :: func1(:,:), func2(:,:)
   complex (kind=CmplxKind), pointer :: sg(:,:)
   complex (kind=CmplxKind) :: v
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
   NumAtoms = getNumAtoms()
!  -------------------------------------------------------------------
   call initProcMapping(isFullPotential(), 'none', 0, NumAtoms)
   call createParallelization()
   call initAtom2Proc(NumAtoms, NumAtoms)
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
!
   if (NumAtoms < 1) then
      call ErrorHandler('main','invalid NumAtoms',NumAtoms)
   endif
!  -------------------------------------------------------------------
   call initPolyhedra(NumAtoms,bravais,'main',0)
!  -------------------------------------------------------------------
!
   allocate(AtomPosition(3,NumAtoms))
   do i=1,NumAtoms
      AtomPosition(1:3,i)=getAtomPosition(i)
   enddo
   do i=1,NumAtoms
!     ----------------------------------------------------------------
      call genPolyhedron(i,i,NumAtoms,AtomPosition)
!     ----------------------------------------------------------------
      if (NumAtoms < NA_limit) then
!        -------------------------------------------------------------
         call printPolyhedron(i)
!        -------------------------------------------------------------
      endif
   enddo
!
   call printPolyhedraTable(getGroupID('Unit Cell'),1,0)
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
   jmax_max = (lmax_max+1)*(lmax_max+2)/2
   kmax_max = (lmax_max+1)**2
!
!  -------------------------------------------------------------------
   call initSphericalHarmonics(2*lmax_max)
   call initGauntFactors(lmax_max, 'none', 0)
!  -------------------------------------------------------------------
   call initStepFunction(NumAtoms,lmax_max,lmax_step,ngr,ngt,'none',0)
!  call initStepFunction(NumAtoms,lmax_step,'none',0)
!  -------------------------------------------------------------------
!
   nr = 601
   allocate( r_mesh(nr), func(nr,jmax_max), func1(nr,kmax_max),       &
             func2(nr,kmax_max))
   func(1:nr,1:jmax_max) = CZERO
   do i=1, min(NumAtoms,NA_limit)
!     ----------------------------------------------------------------
      call printStepFunction(i)
!     ----------------------------------------------------------------
      rs=getOutscrSphRadius(i)
!     ----------------------------------------------------------------
      hr=rs/dble(nr)
      do ir=1,nr
         r_mesh(ir)=hr*ir
         func(ir,1)=TWO*sqrt(PI)
      enddo
!     ----------------------------------------------------------------
      v = getVolumeIntegration(i,nr,r_mesh,1,1,0,func,v_mt)
!     ----------------------------------------------------------------
      write(6,'(a,2d15.8)')'Volume integration = ',v
   enddo
!
   allocate( g(nr), rgc(500), wgc(500) )
!
   do i=1, min(NumAtoms,NA_limit)
!     ----------------------------------------------------------------
      rm=getInscrSphRadius(i)
      rs=getOutscrSphRadius(i)
!     ----------------------------------------------------------------
      hr=(rs-rm)/dble(nr-1)
      write(fname,'(i3)')100+i; fname(1:1)='_'
      fname='stepfunc'//fname(1:3)
      open(unit=10,file=fname,status='unknown',form='formatted')
      do jl = 1,(lmax_step(i)+1)*(lmax_step(i)+2)/2
         do ir=1,nr
            r_mesh(ir)=hr*(ir-1)+rm
!           ----------------------------------------------------------
            func(ir,jl) = getStepFunction(i,jl,r_mesh(ir))
!           ----------------------------------------------------------
            write(10,'(2i4,f15.8,2x,2f15.8)')jl,ir,r_mesh(ir),func(ir,jl)
         enddo
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
   do i=1, min(NumAtoms,NA_limit)
      func1(1:nr,1:kmax_max) = CZERO
      func2(1:nr,1:kmax_max) = CZERO
!     ----------------------------------------------------------------
      rs=getOutscrSphRadius(i)
!     ----------------------------------------------------------------
      hr=rs/dble(nr)
      do ir=1,nr
         r_mesh(ir)=hr*ir
      enddo
      do l = 0, lmax_step(i)
         fac=1-2*mod(l,2)
         kl = (l+1)*(l+1)
         jl = (l+1)*(l+2)/2
         do m = l, -l, -1
            if (m >= 0) then
               do ir=1,nr
!                 ----------------------------------------------------
                  func1(ir,kl) = getStepFunction(i,jl,r_mesh(ir))
!                 ----------------------------------------------------
               enddo
               jl = jl - 1
            else
               do ir=1,nr
                  func1(ir,kl) = conjg(func1(ir,kl-2*m))*fac
               enddo
            endif
            fac = -fac
            kl = kl - 1
         enddo
      enddo
      func2(1:nr,1) = TWO*sqrt(PI)
      kmax_step = (lmax_step(i)+1)**2
!     ----------------------------------------------------------------
      v = getVolumeIntegration(i,nr,r_mesh,kmax_step,0,func1)
!     ----------------------------------------------------------------
      write(6,'(a,2d15.8)')'Volume integration of func1 = ',v
!     ----------------------------------------------------------------
      v = getVolumeIntegration(i,nr,r_mesh,kmax_step,0,func1,kmax_step,0,func2)
!     ----------------------------------------------------------------
      write(6,'(a,2d15.8)')'Volume integration of func1, func2 = ',v
   enddo
!
   deallocate(AtomPosition, r_mesh, func, func1, func2, g)
   deallocate(lmax_step, ngr, ngt)
   deallocate(rgc, wgc)
   nullify ( rg, wg, sg )
!
   if (NumAtoms > NA_limit) then
!     ----------------------------------------------------------------
      call testStepFunction()
!     ----------------------------------------------------------------
   else
      do i=1,NumAtoms
!        -------------------------------------------------------------
         call testStepFunction(i)
!        -------------------------------------------------------------
      enddo
   endif
!
!  -------------------------------------------------------------------
   call endAtom2Proc()
   call endProcMapping()
   call endStepFunction()
   call endGauntFactors()
   call endPolyhedra()
!  -------------------------------------------------------------------
   stop 'Ok'
end program testStepFunc2
