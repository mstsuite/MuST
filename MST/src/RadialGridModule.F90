module RadialGridModule
   use KindParamModule, only : IntKind
   use KindParamModule, only : RealKind
   use ErrorHandlerModule, only : ErrorHandler, StopHandler, WarningHandler
   use PublicTypeDefinitionsModule, only : GridStruct
!
public ::           &
   initRadialGrid,  &
   endRadialGrid,   &
   genRadialGrid,   &
   resetRadialGrid, &
   printRadialGrid, &
   getGrid,         &
   getRindex,       &
   getRmesh,        &
   getNumRmesh,     &
   getMaxNumRmesh,  &
   getRadialGridRadius, &
   getRadialGridPoint,  &
   pushRadialGridToAccel, &
   deleteRadialGridOnAccel
!
   interface genRadialGrid
      module procedure genRadialGrid0, genRadialGrid1, genRadialGrid2, &
                       genRadialGrid3
   end interface
!
private
   type (GridStruct), allocatable, target :: Grid(:)
!
   integer (kind=IntKind) :: NumGrids
!
   logical :: Initialized = .false.
!
   character (len=50) :: stop_routine
!
   integer (kind=IntKind) :: print_level
!
   integer (kind=IntKind) :: MaxNumRmesh
!
   integer (kind=IntKind), parameter :: n_extra=10
!
!  ===================================================================
!  rstart_default = the starting point of the radial grid from the origin
!  xstart_default is set to a value used in the previous version of 
!  KKR/LSMS code and rstart_defalt = exp(xstart_default)
!  In the future release, we will choose rstart_default = 1.0d-5
!  ===================================================================
   real (kind=RealKind), parameter :: xstart_default=-.1113096740000D+02
   real (kind=RealKind), parameter :: rstart_default=exp(xstart_default)
!
contains
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine initRadialGrid(n,istop,iprint)
!  ===================================================================
   implicit none
!
   character (len=*), intent(in) :: istop
!
   integer (kind=IntKind), intent(in) :: n, iprint
   integer (kind=IntKind) :: i
!
   NumGrids = n
   allocate( Grid(n) )
!
   do i=1,n
      Grid(i)%jmt = 0
      Grid(i)%jend = 0
      Grid(i)%nmult = 0
   enddo
!
   Initialized = .true.
!
   stop_routine = istop
   print_level = iprint
!
   MaxNumRmesh = 0
!
   end subroutine initRadialGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine endRadialGrid()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: id
!
   do id=1,NumGrids
      deallocate( Grid(id)%x_mesh, Grid(id)%r_mesh )
   enddo
   deallocate( Grid )
   Initialized = .false.
!
   MaxNumRmesh = 0
!
   end subroutine endRadialGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genRadialGrid0(id,rmt,rinsc,rws,rend,ndivin,ndivout)
!  ===================================================================
   use MathParamModule, only : HALF, ONE, TEN2m6
!
   implicit none
!
   character (len=13), parameter :: sname='genRadialGrid'
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ndivin   ! including rmt
   integer (kind=IntKind), intent(in), optional :: ndivout  ! including rend but 
                                                            ! excluding rmt
!
   real (kind=RealKind), intent(in) :: rmt
   real (kind=RealKind), intent(in) :: rinsc
   real (kind=RealKind), intent(in) :: rws
   real (kind=RealKind), intent(in) :: rend
!
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: jout, jinsc
!
   integer (kind=IntKind), parameter :: ndivout_default = 20
!
   real (kind=RealKind) :: xstart, h, xinsc
!   real (kind=RealKind) :: rstart
   real (kind=RealKind) :: xg
!    
!  ===================================================================
!  genRadialGrid: making grid points along r using double logrithm
!
!        Assuming the following parameters are known
!           rmt     = the muffin-tin sphere radius
!           rinsc   = the inscribed sphere radius of the atomic cell
!           rend    = the bounding sphere radius of the atomic cell
!           ndivin  = an estimate of the number of grid points between 
!                     r0 and rmt including rmt but excluding r0=exp[xstart]
!                     it's the lower limit of jmt.
!           ndivout = the number of grid points between rmt and rcirc
!                     including rcirc but excluding rmt
!            
!        genGrid setup the following parameters
!           xend    = log(rend)
!           xinsc     = log(rinsc)
!           xmt     = log(rmt)
!           xstart  = the starting point
!           hout    = x-step bwteen xmt and xcirc
!           hin     = x-step bwteen xstart and xmt
!           jend    = the number of x or r mesh points
!                     the number of grid points between r0 and rend,
!                     including both r0 and rend
!           jinsc   = the number of grid points between r0 and rinsc
!                     including both rinsc and r0
!           jmt     = the number of grid points between r0 and rmt
!                     including both rmt and r0
!           x_mesh  = x mesh points
!           r_mesh  = r mesh points
!  ===================================================================
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'id out of range',id)
!     ----------------------------------------------------------------
   else if (rend < rmt ) then
!     ----------------------------------------------------------------
      call ErrorHandler('genRadialGrid','rend < rmt',rend,rmt)
!     ----------------------------------------------------------------
   else if (rend < rinsc) then
!     ----------------------------------------------------------------
      call ErrorHandler('genRadialGrid','rend < rinsc',rend,rinsc)
!     ----------------------------------------------------------------
   else if (rinsc < rmt) then
!     ----------------------------------------------------------------
      call ErrorHandler('genRadialGrid','rinsc < rmt',rinsc,rmt)
!     ----------------------------------------------------------------
   else if ( ndivin < 2 ) then
      call ErrorHandler(sname,'Bad parameters: ndivin < 2', ndivin)
   endif
!
   xstart=log(rstart_default)
!
!  Grid(id)%rinsc=rinsc 
!  Grid(id)%rws=rws 
!  Grid(id)%xinsc=log(rinsc) 
!  Grid(id)%jinsc=ndivin
   xinsc=log(rinsc) 
   jinsc=ndivin
   Grid(id)%nmult=1
!
   if (abs(rinsc - rmt) > TEN2m6) then 
      Grid(id)%rmt=rmt
      Grid(id)%xmt=log(rmt)
!     h = (Grid(id)%xinsc-xstart)/real(Grid(id)%jinsc-1,RealKind)
      h = (xinsc-xstart)/real(jinsc-1,RealKind)
!     if (Grid(id)%xinsc >= Grid(id)%xmt ) then 
      if (xinsc >= Grid(id)%xmt ) then 
!        jout = ceiling((Grid(id)%xinsc - Grid(id)%xmt)/h)
         jout = ceiling((xinsc - Grid(id)%xmt)/h)
      else
!        jout = floor((Grid(id)%xinsc - Grid(id)%xmt)/h)
         jout = floor((xinsc - Grid(id)%xmt)/h)
      endif
!     Grid(id)%jmt = Grid(id)%jinsc - jout
      Grid(id)%jmt = jinsc - jout
!     Grid(id)%hin = (Grid(id)%xinsc - Grid(id)%xmt)/real(jout,RealKind)
      Grid(id)%hin = (xinsc - Grid(id)%xmt)/real(jout,RealKind)
!     xstart = Grid(id)%xinsc-(Grid(id)%jinsc-1)*Grid(id)%hin
      xstart = xinsc-(jinsc-1)*Grid(id)%hin
   else
!     Grid(id)%rmt=Grid(id)%rinsc
!     Grid(id)%xmt=Grid(id)%xinsc
      Grid(id)%rmt=rinsc
      Grid(id)%xmt=xinsc
!     Grid(id)%jmt=Grid(id)%jinsc
      Grid(id)%jmt=jinsc
!     Grid(id)%hin=(Grid(id)%xinsc-xstart)/real(Grid(id)%jinsc-1,RealKind)
      Grid(id)%hin=(xinsc-xstart)/real(jinsc-1,RealKind)
   endif
!
   if (abs(rend - rinsc) > TEN2m6) then
      Grid(id)%rend=rend
      Grid(id)%xend=log(rend)
   else
!     Grid(id)%rend=Grid(id)%rinsc
      Grid(id)%rend=rinsc
!     Grid(id)%xend=Grid(id)%xinsc
      Grid(id)%xend=xinsc
   endif
!
   if ( .not.present(ndivout) ) then
      if (abs(rend - rinsc) < TEN2m6) then
!        Grid(id)%jend=Grid(id)%jinsc
         Grid(id)%jend=jinsc
      else
!        Grid(id)%jend=Grid(id)%jinsc+ceiling((Grid(id)%xend-Grid(id)%xinsc)/Grid(id)%hin)
         Grid(id)%jend=jinsc+ceiling((Grid(id)%xend-xinsc)/Grid(id)%hin)
      endif
      Grid(id)%hout=Grid(id)%hin
   else if (abs(rend - rinsc) > TEN2m6) then
      if ( ndivout > 0 ) then
!        Grid(id)%jend=Grid(id)%jinsc+ndivout
         Grid(id)%jend=jinsc+ndivout
!        Grid(id)%hout=(Grid(id)%xend-Grid(id)%xinsc)/real(ndivout,kind=RealKind)
         Grid(id)%hout=(Grid(id)%xend-xinsc)/real(ndivout,kind=RealKind)
      else
!        Grid(id)%jend=Grid(id)%jinsc+                                &
!                      max(ceiling((Grid(id)%xend-Grid(id)%xinsc)/Grid(id)%hin),ndivout_default)
         Grid(id)%jend=jinsc+max(ceiling((Grid(id)%xend-xinsc)/Grid(id)%hin),ndivout_default)
!        Grid(id)%hout=(Grid(id)%xend-Grid(id)%xinsc)/real(Grid(id)%jend-Grid(id)%jinsc,kind=RealKind)
         Grid(id)%hout=(Grid(id)%xend-xinsc)/real(Grid(id)%jend-jinsc,kind=RealKind)
      endif
   else
!      Grid(id)%jend = Grid(id)%jinsc
       Grid(id)%jend = jinsc
       Grid(id)%hout = Grid(id)%hin
   endif
!
!  ===================================================================
!  Re-evaluate xstart.
!  ===================================================================
   Grid(id)%xstart=xstart
   Grid(id)%rstart=exp(Grid(id)%xstart)
   Grid(id)%jend_plus_n=Grid(id)%jend+n_extra
!  -------------------------------------------------------------------
   allocate(Grid(id)%x_mesh(Grid(id)%jend_plus_n))
   allocate(Grid(id)%r_mesh(Grid(id)%jend_plus_n))
!  -------------------------------------------------------------------
!
!  ===================================================================
!  set up x-mesh and r-mesh........................................
!  ===================================================================
!  do j=1,Grid(id)%jinsc
   do j=1,jinsc
      xg=Grid(id)%xstart+(j-1)*Grid(id)%hin
      Grid(id)%x_mesh(j)=xg
      Grid(id)%r_mesh(j)=exp(xg)
   enddo
!  jout=Grid(id)%jend_plus_n-Grid(id)%jinsc
   jout=Grid(id)%jend_plus_n-jinsc
   do j=1,jout
!     xg=Grid(id)%xinsc+j*Grid(id)%hout
      xg=xinsc+j*Grid(id)%hout
!     Grid(id)%x_mesh(j+Grid(id)%jinsc)=xg
      Grid(id)%x_mesh(j+jinsc)=xg
!     Grid(id)%r_mesh(j+Grid(id)%jinsc)=exp(xg)
      Grid(id)%r_mesh(j+jinsc)=exp(xg)
   enddo
!
   MaxNumRmesh = max(MaxNumRmesh, Grid(id)%jend_plus_n)
!
   end subroutine genRadialGrid0
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genRadialGrid1(id,rmt,rinsc,rws,rend,ndivin,ndivout,nmult)
!  ===================================================================
   use MathParamModule, only : HALF, ONE, TEN2m6
!
   implicit none
!
   character (len=13), parameter :: sname='genRadialGrid'
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ndivin   ! including rmt
   integer (kind=IntKind), intent(in) :: ndivout  ! including rend but 
                                                  ! excluding rmt
   integer (kind=IntKind), intent(in) :: nmult    ! = hin/hout
!
   real (kind=RealKind), intent(in) :: rmt
   real (kind=RealKind), intent(in) :: rinsc
   real (kind=RealKind), intent(in) :: rws
   real (kind=RealKind), intent(in) :: rend
!
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: jout, jinsc
!
   real (kind=RealKind) :: xstart, xinsc
!   real (kind=RealKind) :: rstart
   real (kind=RealKind) :: xg
!    
!  ===================================================================
!  genRadialGrid: making grid points along r using double logrithm
!
!        Assuming the following parameters are known
!           rend    = the bounding sphere radius of the atomic cell
!           rinsc   = the inscribed sphere radius of the atomic cell
!           rmt     = the muffin-tin sphere radius
!           ndivout = the number of grid points between rmt and rcirc
!                     including rcirc but excluding rmt
!           ndivin  = an estimate of the number of grid points between 
!                     r0 and rmt including rmt but excluding r0=exp[xstart]
!                     it's the lower limit of jmt.
!           nmult   = an estimate of the multiplier as the x-step outside 
!                     compared with the x-step inside
!            
!        genGrid setup the following parameters
!           xend    = log(rend)
!           xinsc     = log(rinsc)
!           xmt     = log(rmt)
!           xstart  = the starting point
!           hout    = x-step bwteen xmt and xcirc
!           hin     = x-step bwteen xstart and xmt
!           jend    = the number of x or r mesh points
!                     the number of grid points between r0 and rend,
!                     including both r0 and rend
!           jinsc   = the number of grid points between r0 and rinsc
!                     including both rinsc and r0
!           jmt     = the number of grid points between r0 and rmt
!                     including both rmt and r0
!           x_mesh  = x mesh points
!           r_mesh  = r mesh points
!  ===================================================================
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'id out of range',id)
!     ----------------------------------------------------------------
   else if (rend < rmt) then
!     ----------------------------------------------------------------
      call ErrorHandler('genRadialGrid','rend < rmt',rend,rmt)
!     ----------------------------------------------------------------
   else if( ndivin<=0 .or. (ndivin<=0 .and. nmult<=0) .or. &
            ndivin == 1) then
      print *,'ndivout = ',ndivout
      print *,'ndivin  = ',ndivin
      print *,'nmult   = ',nmult
      call ErrorHandler(sname,'Bad parameters: ndivin, ndivout, nmult')
   else if (rend - rmt < TEN2m6 .and. ndivin <= 0) then
      call ErrorHandler(sname,'if rend = rmt, we should have ndivin > 1')
   endif
!
   xstart=log(rstart_default)
!
   Grid(id)%rend=rend
   Grid(id)%xend=log(rend)
!  Grid(id)%rinsc=rinsc 
!  Grid(id)%xinsc=log(rinsc) 
!  Grid(id)%rws=rws
!  Grid(id)%jinsc=ndivin
   xinsc=log(rinsc) 
   jinsc=ndivin
   Grid(id)%rmt=rmt
   Grid(id)%xmt=log(rmt)
!
   if (ndivin<=0) then
      call WarningHandler("genRadialGrid1 ::", &
           "ndivin is an input parameter:: default to 1001 now!!!")
!     Grid(id)%jinsc=1001
      jinsc=1001
   endif
   if (abs(rend - rinsc) < TEN2m6) then
      if (ndivin<=0) then
         call ErrorHandler("genRadialGrid1 ::", &
              "ndivin has to be an input parameter!!!")
      endif
!     Grid(id)%jend=Grid(id)%jinsc
      Grid(id)%jend=jinsc
!     Grid(id)%jmt=Grid(id)%jinsc
      Grid(id)%jmt=jinsc
!     Grid(id)%rmt=Grid(id)%rinsc
      Grid(id)%rmt=rinsc
!     Grid(id)%xmt=Grid(id)%xinsc
      Grid(id)%xmt=xinsc
      Grid(id)%nmult=1
!     Grid(id)%hin=(Grid(id)%xinsc-xstart)/(ndivin-ONE)
      Grid(id)%hin=(xinsc-xstart)/(ndivin-ONE)
      Grid(id)%hout=Grid(id)%hin
   else if ( abs(rinsc - rmt) < TEN2m6 ) then
!     Grid(id)%jmt=Grid(id)%jinsc
      Grid(id)%jmt=jinsc
!     Grid(id)%rmt=Grid(id)%rinsc
      Grid(id)%rmt=rinsc
!     Grid(id)%xmt=Grid(id)%xinsc
      Grid(id)%xmt=xinsc
      if ( ndivout > 0 ) then
         Grid(id)%jend=Grid(id)%jmt+ndivout
         if ( nmult>0 ) then
            Grid(id)%hout=(Grid(id)%xend-Grid(id)%xmt)/dble(ndivout)
            Grid(id)%nmult=nmult
            Grid(id)%hin=nmult*Grid(id)%hout
         else if ( nmult==0 ) then
            Grid(id)%hout=(Grid(id)%xend-Grid(id)%xmt)/dble(ndivout)
!           Grid(id)%hin=(Grid(id)%xmt-xstart)/Grid(id)%jinsc
            Grid(id)%hin=(Grid(id)%xmt-xstart)/real(jinsc,kind=RealKind)
            Grid(id)%nmult=int(Grid(id)%hin/Grid(id)%hout)
            Grid(id)%hin=Grid(id)%nmult*Grid(id)%hout
         else
            Grid(id)%hout=(Grid(id)%rend-Grid(id)%rmt)/dble(ndivout)
!           Grid(id)%hin=(Grid(id)%xmt-xstart)/Grid(id)%jinsc
            Grid(id)%hin=(Grid(id)%xmt-xstart)/real(jinsc,kind=RealKind)
            Grid(id)%nmult = -ndivout
         endif
      else
         Grid(id)%hin=(Grid(id)%xmt-xstart)/ndivin
         if ( nmult>0 ) then
            Grid(id)%jend = &
            nmult*int((Grid(id)%xend-Grid(id)%xmt)/Grid(id)%hin+.5)
            Grid(id)%hout= (Grid(id)%xend-Grid(id)%xmt)/Grid(id)%jend
            Grid(id)%hin = nmult*Grid(id)%hout
            Grid(id)%nmult = nmult
            Grid(id)%jend = Grid(id)%jend+Grid(id)%jmt
         else if ( nmult==0 ) then
            Grid(id)%hout= (Grid(id)%xend-Grid(id)%xmt)/              &
                    int((Grid(id)%xend-Grid(id)%xmt)/Grid(id)%hin +.5)
            Grid(id)%hin = Grid(id)%hout
            Grid(id)%nmult = 1 
            Grid(id)%jend = Grid(id)%jmt+ &
                    int((Grid(id)%xend-Grid(id)%xmt)/Grid(id)%hin +.5)
         else
            Grid(id)%hout=(Grid(id)%rend-Grid(id)%rmt)/dble(50)
!           Grid(id)%hin=(Grid(id)%xmt-xstart)/Grid(id)%jinsc
            Grid(id)%hin=(Grid(id)%xmt-xstart)/real(jinsc,kind=RealKind)
            Grid(id)%nmult = -50
            Grid(id)%jend = Grid(id)%jmt+50
         endif
      endif
   else
!      Grid(id)%hin=(Grid(id)%xinsc-xstart)/Grid(id)%jinsc
       Grid(id)%hin=(xinsc-xstart)/real(jinsc,kind=RealKind)
!      Grid(id)%hout= (Grid(id)%xinsc-Grid(id)%xmt)/              &
!               int(abs(Grid(id)%xinsc-Grid(id)%xmt)/Grid(id)%hin +.5)
       Grid(id)%hout= (xinsc-Grid(id)%xmt)/int(abs(xinsc-Grid(id)%xmt)/Grid(id)%hin +.5)
!      jout = (Grid(id)%xinsc-Grid(id)%xmt)/Grid(id)%hout
       jout = (xinsc-Grid(id)%xmt)/Grid(id)%hout
!      Grid(id)%jmt = Grid(id)%jinsc-jout
       Grid(id)%jmt = jinsc-jout
       Grid(id)%hin = Grid(id)%hout
       if ( nmult>=0 ) then
          Grid(id)%nmult = -max(nmult,50)
       else
          Grid(id)%nmult = min(-50,nmult)
       endif
       if ( ndivout>0 ) then
          Grid(id)%nmult = -max(-Grid(id)%nmult,ndivout)
       else
          Grid(id)%nmult = min(-50,ndivout)
       endif
!      Grid(id)%jend = Grid(id)%jinsc - Grid(id)%nmult
       Grid(id)%jend = jinsc - Grid(id)%nmult
!      Grid(id)%hout = (Grid(id)%rend-Grid(id)%rinsc)/(-Grid(id)%nmult)
       Grid(id)%hout = (Grid(id)%rend-rinsc)/real(-Grid(id)%nmult,kind=RealKind)
         
!        if rmt < rinsc assure that they both fall on the radial grid
!        the inner grid is generated relative to rinsc not the rmt
!
!         call ErrorHandler("genRadialGrid1 ::", &
!              "Rmt/=Rinsc have not been implemented yet")
!      endif
   endif
!
!  ===================================================================
!  Re-evaluate xstart.
!  ===================================================================
!  Grid(id)%xstart=Grid(id)%xinsc-(Grid(id)%jinsc-ONE)*Grid(id)%hin
   Grid(id)%xstart=xinsc-(jinsc-ONE)*Grid(id)%hin
   Grid(id)%rstart=exp(Grid(id)%xstart)
   Grid(id)%jend_plus_n=Grid(id)%jend+n_extra*abs(Grid(id)%nmult)
!  -------------------------------------------------------------------
   allocate(Grid(id)%x_mesh(Grid(id)%jend_plus_n))
   allocate(Grid(id)%r_mesh(Grid(id)%jend_plus_n))
!  -------------------------------------------------------------------
!
!  ===================================================================
!  set up x-mesh and r-mesh........................................
!  ===================================================================
!  do j=1,Grid(id)%jinsc
   do j=1,jinsc
      xg=Grid(id)%xstart+(j-1)*Grid(id)%hin
      Grid(id)%x_mesh(j)=xg
      Grid(id)%r_mesh(j)=exp(xg)
   enddo
!  jout=Grid(id)%jend_plus_n-Grid(id)%jinsc
   jout=Grid(id)%jend_plus_n-jinsc
   if ( Grid(id)%nmult >=0 ) then
      do j=1,jout
!        xg=Grid(id)%xinsc+j*Grid(id)%hout
         xg=xinsc+j*Grid(id)%hout
         Grid(id)%x_mesh(j+Grid(id)%jmt)=xg
         Grid(id)%r_mesh(j+Grid(id)%jmt)=exp(xg)
      enddo
   else
      do j=1,jout
!        xg=Grid(id)%rinsc+j*Grid(id)%hout
         xg=rinsc+j*Grid(id)%hout
         Grid(id)%x_mesh(j+Grid(id)%jmt)=xg
         Grid(id)%r_mesh(j+Grid(id)%jmt)=xg
      enddo
   endif
!
   if ( rinsc - rmt < TEN2m6 ) then
!     Grid(id)%jinsc=Grid(id)%jmt
      jinsc=Grid(id)%jmt
   else
!     ----------------------------------------------------------------
!     call hunt(Grid(id)%jend_plus_n,Grid(id)%r_mesh,rinsc,Grid(id)%jinsc)
      call hunt(Grid(id)%jend_plus_n,Grid(id)%r_mesh,rinsc,jinsc)
!     ----------------------------------------------------------------
   endif
!
   MaxNumRmesh = max(MaxNumRmesh, Grid(id)%jend_plus_n)
!
   end subroutine genRadialGrid1
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genRadialGrid2(id,xstart,rmt,rinsc,rws,rend,jinsc,ndivout)
!  ===================================================================
   use MathParamModule, only : HALF, ONE, TEN2m6
!
!  *******************************************************************
!  generate single logrithmic grid
!  *******************************************************************
!
   implicit none
!
   character (len=13), parameter :: sname='genRadialGrid'
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: jinsc      ! including rmt
   integer (kind=IntKind), optional :: ndivout
   integer (kind=IntKind) :: j, jout
!
   real (kind=RealKind), intent(in) :: xstart
   real (kind=RealKind), intent(in) :: rmt
   real (kind=RealKind), intent(in) :: rinsc
   real (kind=RealKind), intent(in) :: rws
   real (kind=RealKind), intent(in) :: rend
   real (kind=RealKind) :: xg, xinsc
!    
!  ===================================================================
!  genRadialGrid: making grid points along r using single logrithm
!
!        Assuming the following parameters are known
!           rend    = the bounding sphere radius of the grid
!           rinsc   = the inscribed sphere radius of the grid
!           rmt     = the muffin-tin or ASA sphere radius
!           jmt     = the number of grid points between r0 and rmt
!                     including both rmt and r0
!            
!        genGrid setup the following parameters
!           xend    = log(rend)
!           xinsc   = log(rinsc)
!           xmt     = log(rmt)
!           xstart  = the starting point
!           hin     = x-step bwteen xstart and xmt
!           hout    = hin
!           jend    = the number of x or r mesh points
!                     the number of grid points between r0 and rend+epsilon,
!                     including both r0 and rend+epsilon
!           jinsc   = the number of x or r mesh points
!                     the number of grid points between r0 and rinsc,
!                     including both r0 and rinsc+epsilon
!           x_mesh  = x mesh points
!           r_mesh  = r mesh points
!  ===================================================================
!
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'id out of range',id)
!     ----------------------------------------------------------------
   endif
!
   Grid(id)%xstart=xstart
   Grid(id)%rstart=exp(xstart)
!
   if (Grid(id)%rstart >= rmt) then
!     ----------------------------------------------------------------
      call ErrorHandler('genRadialGrid','rstart >= rmt',Grid(id)%rstart,rmt)
!     ----------------------------------------------------------------
   else if (rend < rmt) then
!     ----------------------------------------------------------------
      call ErrorHandler('genRadialGrid','rend < rmt',rend,rmt)
!     ----------------------------------------------------------------
   else if ( jinsc <= 0) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'Bad parameters: jinsc', jinsc)
!     ----------------------------------------------------------------
   endif
!
!  Grid(id)%jinsc=jinsc
!  Grid(id)%rinsc=rinsc
!  Grid(id)%xinsc=log(rinsc)
!  Grid(id)%rws  =rws
   Grid(id)%rmt  =rmt
   xinsc=log(rinsc) 
   Grid(id)%xmt=log(rmt)
   Grid(id)%nmult=1
   if ( .not.present(ndivout) ) then
      if (abs(rend - rinsc) < TEN2m6) then
!        Grid(id)%hin=(Grid(id)%xinsc-xstart)/real(jinsc-1,RealKind)
         Grid(id)%hin=(xinsc-xstart)/real(jinsc-1,RealKind)
!        Grid(id)%jend=Grid(id)%jinsc
         Grid(id)%jend=jinsc
!        Grid(id)%jmt=Grid(id)%jinsc
         Grid(id)%jmt=jinsc
!        Grid(id)%rmt=Grid(id)%rinsc
         Grid(id)%rmt=rinsc
!        Grid(id)%xmt=Grid(id)%xinsc
         Grid(id)%xmt=xinsc
      else
         if ( abs(rinsc - rmt) < TEN2m6 ) then
!           Grid(id)%hin=(Grid(id)%xinsc-xstart)/real(jinsc-1,RealKind)
            Grid(id)%hin=(xinsc-xstart)/real(jinsc-1,RealKind)
!           Grid(id)%jmt=Grid(id)%jinsc
            Grid(id)%jmt=jinsc
!           Grid(id)%rmt=Grid(id)%rinsc
            Grid(id)%rmt=rinsc
!           Grid(id)%xmt=Grid(id)%xinsc
            Grid(id)%xmt=xinsc
            jout=ceiling( (log(rend)-Grid(id)%xmt)/Grid(id)%hin )
!           Grid(id)%jend=Grid(id)%jinsc+jout
            Grid(id)%jend=jinsc+jout
         else
!
!        If rmt < rinsc then redefine xstart and the number of grid points
!        such that both, rmt and rinsc, fall on the radial grid.
!
!           Grid(id)%hin=(Grid(id)%xinsc-xstart)/real(jinsc-1,RealKind)
            Grid(id)%hin=(xinsc-xstart)/real(jinsc-1,RealKind)
            jout=max(floor( (log(rinsc)-Grid(id)%xmt)/Grid(id)%hin ),1)
            Grid(id)%hin = (log(rinsc)-Grid(id)%xmt)/jout
!           Grid(id)%xstart=Grid(id)%xinsc-(Grid(id)%jinsc-ONE)*Grid(id)%hin
            Grid(id)%xstart=xinsc-(jinsc-ONE)*Grid(id)%hin
!           Grid(id)%jmt = Grid(id)%jinsc-jout
            Grid(id)%jmt = jinsc-jout
            Grid(id)%rstart=exp(xstart)
            jout=ceiling( (log(rend)-log(rinsc))/Grid(id)%hin )
!           Grid(id)%jend=Grid(id)%jinsc+jout
            Grid(id)%jend=jinsc+jout
         endif
      endif
      Grid(id)%hout=Grid(id)%hin
   else if (ndivout<0) then
      if (abs(rend - rinsc) < TEN2m6) then
!        Grid(id)%hin=(Grid(id)%xinsc-xstart)/real(jinsc-1,RealKind)
         Grid(id)%hin=(xinsc-xstart)/real(jinsc-1,RealKind)
!        Grid(id)%jend=Grid(id)%jinsc
         Grid(id)%jend=jinsc
!        Grid(id)%jmt=Grid(id)%jinsc
         Grid(id)%jmt=jinsc
!        Grid(id)%rmt=Grid(id)%rinsc
         Grid(id)%rmt=rinsc
!        Grid(id)%xmt=Grid(id)%xinsc
         Grid(id)%xmt=xinsc
         Grid(id)%hout=Grid(id)%hin
         Grid(id)%nmult=1
      else
         if ( abs(rinsc - rmt) < TEN2m6 ) then
!           Grid(id)%hin=(Grid(id)%xinsc-xstart)/real(jinsc-1,RealKind)
            Grid(id)%hin=(xinsc-xstart)/real(jinsc-1,RealKind)
!           Grid(id)%jmt=Grid(id)%jinsc
            Grid(id)%jmt=jinsc
!           Grid(id)%rmt=Grid(id)%rinsc
            Grid(id)%rmt=rinsc
!           Grid(id)%xmt=Grid(id)%xinsc
            Grid(id)%xmt=xinsc
            Grid(id)%hout=(Grid(id)%rend-rinsc)/(-ndivout)
!           Grid(id)%jend=Grid(id)%jinsc-ndivout
            Grid(id)%jend=jinsc-ndivout
            Grid(id)%nmult=-ndivout
         else
!
!        If rmt < rinsc then redefine xstart and the number of grid points
!        such that both, rmt and rinsc, fall on the radial grid.
!
!           Grid(id)%hin=(Grid(id)%xinsc-xstart)/real(jinsc-1,RealKind)
            Grid(id)%hin=(xinsc-xstart)/real(jinsc-1,RealKind)
            jout=max(floor( (log(rinsc)-Grid(id)%xmt)/Grid(id)%hin ),1)
            Grid(id)%hin = (log(rinsc)-Grid(id)%xmt)/jout
!           Grid(id)%xstart=Grid(id)%xinsc-(Grid(id)%jinsc-ONE)*Grid(id)%hin
            Grid(id)%xstart=xinsc-(jinsc-ONE)*Grid(id)%hin
!           Grid(id)%jmt = Grid(id)%jinsc-jout
            Grid(id)%jmt = jinsc-jout
            Grid(id)%rstart=exp(xstart)
            Grid(id)%hout=(Grid(id)%rend-rinsc)/(-ndivout)
!           Grid(id)%jend=Grid(id)%jinsc-ndivout
            Grid(id)%jend=jinsc-ndivout
            Grid(id)%nmult=-ndivout
         endif
      endif      
   endif
!
   Grid(id)%jend_plus_n=Grid(id)%jend+n_extra
!
!  -------------------------------------------------------------------
   allocate(Grid(id)%x_mesh(Grid(id)%jend_plus_n))
   allocate(Grid(id)%r_mesh(Grid(id)%jend_plus_n))
!  -------------------------------------------------------------------
!
!  ===================================================================
!  set up x-mesh and r-mesh........................................
!  ===================================================================
   
!  do j=1,Grid(id)%jinsc
   do j=1,jinsc
      xg=Grid(id)%xstart+(j-1)*Grid(id)%hin
      Grid(id)%x_mesh(j)=xg
      Grid(id)%r_mesh(j)=exp(xg)
   enddo
!  jout=Grid(id)%jend_plus_n-Grid(id)%jinsc
   jout=Grid(id)%jend_plus_n-jinsc
   if ( Grid(id)%nmult >=0 ) then
      do j=1,jout
!        xg=Grid(id)%xinsc+j*Grid(id)%hout
         xg=xinsc+j*Grid(id)%hout
         Grid(id)%x_mesh(j+Grid(id)%jmt)=xg
         Grid(id)%r_mesh(j+Grid(id)%jmt)=exp(xg)
      enddo
   else
      do j=1,jout
!        xg=Grid(id)%rinsc+j*Grid(id)%hout
         xg=rinsc+j*Grid(id)%hout
         Grid(id)%x_mesh(j+Grid(id)%jmt)=xg
         Grid(id)%r_mesh(j+Grid(id)%jmt)=xg
      enddo
   endif
!
   Grid(id)%rend=Grid(id)%r_mesh(Grid(id)%jend)
   Grid(id)%xend=Grid(id)%x_mesh(Grid(id)%jend)
!
   MaxNumRmesh = max(MaxNumRmesh, Grid(id)%jend_plus_n)
!
   end subroutine genRadialGrid2
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine genRadialGrid3(id,rmt,rend,ndivin,rstart)
!  ===================================================================
   use MathParamModule, only : HALF, ONE, TEN2m6
!
   implicit none
!
   character (len=13), parameter :: sname='genRadialGrid'
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind), intent(in) :: ndivin   ! including rmt
!
   real (kind=RealKind), intent(in) :: rmt
   real (kind=RealKind), intent(in) :: rend
   real (kind=RealKind), intent(in), optional :: rstart   ! The starting point
                                                          ! of the radial grid
!
   integer (kind=IntKind) :: j, ndivout
!
   real (kind=RealKind) :: xend, xg
!    
!  ===================================================================
!  genRadialGrid3: making grid points along r using double logrithm
!
!        Assuming the following parameters are known
!           rmt     = the muffin-tin sphere radius
!           rend    = an end point estimate of the radial grid.
!                     Note: rend may not be on the radial grid point
!           ndivin  = the number of grid points from 
!                     rstart to rmt, including both rmt and rstart
!           rstart  = the starting point of the radial grid (optional)
!                     It will take the default value if it is not given from
!                     the input
!            
!        genGrid setup the following parameters
!           xend    = log(rend)
!           xmt     = log(rmt)
!           xstart  = the starting point
!           hout    = x-step bwteen xmt and xcirc
!           hin     = x-step bwteen xstart and xmt
!           jend    = the number of x or r mesh points
!                     the number of grid points between r0 and rend,
!                     including both r0 and rend
!           jmt     = the number of grid points between r0 and rmt
!                     including both rmt and r0
!           x_mesh  = x mesh points
!           r_mesh  = r mesh points
!  ===================================================================
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler(sname,'id out of range',id)
!     ----------------------------------------------------------------
   else if (rend < rmt ) then
!     ----------------------------------------------------------------
      call ErrorHandler('genRadialGrid','rend < rmt',rend,rmt)
!     ----------------------------------------------------------------
   else if ( ndivin < 2 ) then
      call ErrorHandler(sname,'Bad parameters: ndivin < 2', ndivin)
   endif
!
   if (present(rstart)) then
      if (rstart > 1.0d-5) then
         call WarningHandler(sname,'rstart > 0.00001',rstart)
      endif
      Grid(id)%rstart = rstart
      Grid(id)%xstart = log(rstart)
   else
      Grid(id)%rstart = rstart_default
      Grid(id)%xstart = log(rstart_default)
   endif
!
   Grid(id)%rmt = rmt
   Grid(id)%xmt = log(rmt)
   Grid(id)%jmt = ndivin
   Grid(id)%hin = (Grid(id)%xmt-Grid(id)%xstart)/real(Grid(id)%jmt-1,RealKind)
   Grid(id)%nmult = 1
   Grid(id)%hout = Grid(id)%hin
!
   xend = log(rend)
   ndivout = ceiling((xend-Grid(id)%xmt)/Grid(id)%hout)
   Grid(id)%xend = Grid(id)%xmt+ndivout*Grid(id)%hout
   Grid(id)%rend = exp(Grid(id)%xend)
   if (Grid(id)%rend < rend) then
      call ErrorHandler(sname,'Grid(id)%rend < rend',Grid(id)%rend,rend)
   endif
   Grid(id)%jend = Grid(id)%jmt + ndivout
   Grid(id)%jend_plus_n = Grid(id)%jend+n_extra
!
!  -------------------------------------------------------------------
   allocate(Grid(id)%x_mesh(Grid(id)%jend_plus_n))
   allocate(Grid(id)%r_mesh(Grid(id)%jend_plus_n))
!  -------------------------------------------------------------------
!
!  ===================================================================
!  set up x-mesh and r-mesh........................................
!  ===================================================================
   do j = 1, Grid(id)%jmt
      xg = Grid(id)%xstart + (j-1)*Grid(id)%hin
      Grid(id)%x_mesh(j) = xg
      Grid(id)%r_mesh(j) = exp(xg)
   enddo
   do j = Grid(id)%jmt+1, Grid(id)%jend_plus_n
      xg = Grid(id)%xmt + (j-Grid(id)%jmt)*Grid(id)%hout
      Grid(id)%x_mesh(j) = xg
      Grid(id)%r_mesh(j) = exp(xg)
   enddo
!
   MaxNumRmesh = max(MaxNumRmesh, Grid(id)%jend_plus_n)
!
   end subroutine genRadialGrid3
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine resetRadialGrid(id,rmt,rinsc,rws,rend)
!  ===================================================================
   use MathParamModule, only : HALF, ONE, TEN2m6
!
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
!
!   integer (kind=IntKind) :: ndivin   ! including rmt
!   integer (kind=IntKind) :: ndivout  ! including rend but 
                                                  ! excluding rmt
!   integer (kind=IntKind) :: nmult    ! = hin/hout
!
   real (kind=RealKind), intent(in) :: rmt
   real (kind=RealKind), intent(in) :: rinsc
   real (kind=RealKind), intent(in) :: rws
   real (kind=RealKind), intent(in) :: rend
!
   integer (kind=IntKind) :: j
   integer (kind=IntKind) :: jout
!
   real (kind=RealKind) :: xg
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('resetRadialGrid','need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler('resetRadialGrid','id out of range',id)
!     ----------------------------------------------------------------
   else if (rend < rmt) then
!     ----------------------------------------------------------------
      call ErrorHandler('resetRadialGrid','rend < rmt',rend,rmt)
!     ----------------------------------------------------------------
   else if(Grid(id)%jmt<=0 .or. Grid(id)%jend<=0 .or. Grid(id)%nmult<=0) then
      print *,'jmt = ',Grid(id)%jmt
      print *,'jend  = ',Grid(id)%jend
      print *,'nmult   = ',Grid(id)%nmult
!     ----------------------------------------------------------------
      call ErrorHandler('resetRadialGrid','invalid initial radial grid')
!     ----------------------------------------------------------------
   endif
!
   Grid(id)%rend=rend
   Grid(id)%xend=log(rend)
!  Grid(id)%rinsc=rinsc
!  Grid(id)%xinsc=log(rinsc)
   Grid(id)%rmt=rmt
!  Grid(id)%rws=rws
   Grid(id)%xmt=log(rmt)
   if (rend - rmt < TEN2m6) then
      Grid(id)%jend_plus_n=Grid(id)%jmt+Grid(id)%jend_plus_n-Grid(id)%jend
      Grid(id)%jend=Grid(id)%jmt
!     Grid(id)%jinsc=Grid(id)%jmt
!     Grid(id)%rinsc=rmt
!     Grid(id)%xinsc=log(rmt)
      Grid(id)%nmult=1
      Grid(id)%hin=(Grid(id)%xmt-Grid(id)%xstart)/(Grid(id)%jmt-ONE)
      Grid(id)%hout=Grid(id)%hin
   else
!     if ( rinsc - rmt < TEN2m6 ) then
!        Grid(id)%jinsc=Grid(id)%jmt
!        Grid(id)%rinsc=rmt
!        Grid(id)%xinsc=log(rmt)
!     else
!        call ErrorHandler("resetRadialGrid", &
!             "Rmt/=Rinsc have not been implemented yet")
!     endif
      Grid(id)%hout=(Grid(id)%xend-Grid(id)%xmt)/dble(Grid(id)%jend- &
                                                      Grid(id)%jmt)
      Grid(id)%hin=Grid(id)%nmult*Grid(id)%hout
   endif
!
   Grid(id)%xstart=Grid(id)%xmt-(Grid(id)%jmt-ONE)*Grid(id)%hin
   Grid(id)%rstart=exp(Grid(id)%xstart)
!
   if (.not.associated(Grid(id)%x_mesh)) then
!     ----------------------------------------------------------------
      allocate(Grid(id)%x_mesh(Grid(id)%jend_plus_n))
!     ----------------------------------------------------------------
   endif
   if (.not.associated(Grid(id)%r_mesh)) then
!     ----------------------------------------------------------------
      allocate(Grid(id)%r_mesh(Grid(id)%jend_plus_n))
!     ----------------------------------------------------------------
   endif
!
!  ===================================================================
!  set up x-mesh and r-mesh........................................
!  ===================================================================
   do j=1,Grid(id)%jmt
      xg=Grid(id)%xstart+(j-1)*Grid(id)%hin
      Grid(id)%x_mesh(j)=xg
      Grid(id)%r_mesh(j)=exp(xg)
   enddo
   jout=Grid(id)%jend_plus_n-Grid(id)%jmt
   do j=1,jout
      xg=Grid(id)%xmt+j*Grid(id)%hout
      Grid(id)%x_mesh(j+Grid(id)%jmt)=xg
      Grid(id)%r_mesh(j+Grid(id)%jmt)=exp(xg)
   enddo
!
   end subroutine resetRadialGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine printRadialGrid(id)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('printRadialGrid','need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler('printRadialGrid','invalid grid index',id)
!     ----------------------------------------------------------------
   endif
!
   write(6,'(/,80(''-''))')
   write(6,'(/,24x,a)')' *******************************'
   write(6,'(24x,a)')  ' * Output from printRadialGrid *'
   write(6,'(24x,a,/)')' *******************************'
   write(6,'(a)')'============================================================'
   write(6,'(  ''Local atom index'',t40,''='',i5)')    id
   write(6,'(  ''jmt        '',t40,''='',i5)')    Grid(id)%jmt
!  write(6,'(  ''jinsc      '',t40,''='',i5)')    Grid(id)%jinsc
!  write(6,'(  ''jws        '',t40,''='',i5)')    Grid(id)%jws  
   write(6,'(  ''jend       '',t40,''='',i5)')    Grid(id)%jend 
   write(6,'(  ''jend_plus_n'',t40,''='',i5)')    Grid(id)%jend_plus_n
   write(6,'(  ''hin        '',t40,''='',d20.13)')Grid(id)%hin
   write(6,'(  ''hout       '',t40,''='',d20.13)')Grid(id)%hout
   write(6,'(  ''nmult      '',t40,''='',i5)')    Grid(id)%nmult
   write(6,'(  ''xstart     '',t40,''='',d20.13)')Grid(id)%xstart
   write(6,'(  ''xmt        '',t40,''='',d20.13)')Grid(id)%xmt
   write(6,'(  ''xend       '',t40,''='',d20.13)')Grid(id)%xend
   write(6,'(  ''rstart     '',t40,''='',d20.13)')Grid(id)%rstart
   write(6,'(  ''rmt        '',t40,''='',d20.13)')Grid(id)%rmt
!  write(6,'(  ''rinsc      '',t40,''='',d20.13)')Grid(id)%rinsc
   write(6,'(  ''rend       '',t40,''='',d20.13)')Grid(id)%rend 
!  write(6,'(  ''rws        '',t40,''='',d20.13)')Grid(id)%rws 
   write(6,'(a)')'============================================================'
!
   end subroutine printRadialGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getGrid(id) result(gp)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   type (GridStruct), pointer :: gp
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getGrid','need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler('genGrid','id out of range',id)
!     ----------------------------------------------------------------
   endif
!
   gp=>Grid(id)
   end function getGrid
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRindex(id,r) result(i)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: i
!
   real (kind=RealKind), intent(in) :: r
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getRindex','need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler('genRindex','id out of range',id)
!     ----------------------------------------------------------------
   endif
   if (r <= Grid(id)%rstart) then
      i=1
   else if (r >= Grid(id)%r_mesh(Grid(id)%jend_plus_n)) then
      i=Grid(id)%jend_plus_n
   else
!     ----------------------------------------------------------------
      call hunt(Grid(id)%jend_plus_n,Grid(id)%r_mesh,r,i)
!     ----------------------------------------------------------------
   endif
!
   end function getRindex
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getNumRmesh(id) result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getNumRmesh','need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler('getNumRmesh','id out of range',id)
!     ----------------------------------------------------------------
   endif
!
   n = Grid(id)%jend
!
   end function getNumRmesh
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getMaxNumRmesh() result(n)
!  ===================================================================
   implicit none
   integer (kind=IntKind) :: n
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getMaxNumRmesh','need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (MaxNumRmesh < 1) then
!     ----------------------------------------------------------------
      call WarningHandler('getMaxNumRmesh','No mesh has been generated')
!     ----------------------------------------------------------------
   endif
!
   n = MaxNumRmesh
!
   end function getMaxNumRmesh
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRmesh(id) result(r)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: n
!
   real (kind=RealKind), pointer :: r(:)
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getRmesh','need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler('getRmesh','id out of range',id)
!     ----------------------------------------------------------------
   endif
!
   n = Grid(id)%jend_plus_n
   r => Grid(id)%r_mesh(1:n)
!
   end function getRmesh
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRadialGridRadius(id) result(r)
!  ===================================================================
   implicit none
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: n
!
   real (kind=RealKind) :: r
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getRadialGridRadius','need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler('getRadialGridRadius','id out of range',id)
!     ----------------------------------------------------------------
   endif
!
   n = Grid(id)%jend_plus_n
   r = Grid(id)%r_mesh(n)
!
   end function getRadialGridRadius
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   function getRadialGridPoint(id,r) result(jr)
!  ===================================================================
   use MathParamModule, only : ZERO
   implicit none
!
   integer (kind=IntKind), intent(in) :: id
   integer (kind=IntKind) :: j, jr
!
   real (kind=RealKind), intent(in) :: r
!
   if (.not.Initialized) then
!     ----------------------------------------------------------------
      call ErrorHandler('getRadialGridPoint','need to call initRadialGrid first')
!     ----------------------------------------------------------------
   else if (id < 1 .or. id > NumGrids) then
!     ----------------------------------------------------------------
      call ErrorHandler('getRadialGridPoint','id out of range',id)
!     ----------------------------------------------------------------
   else if (r <= ZERO) then
!     ----------------------------------------------------------------
      call ErrorHandler('getRadialGridPoint','r <= 0',r)
!     ----------------------------------------------------------------
   endif
!
   jr = Grid(id)%jend
   LOOP_j: do j = Grid(id)%jend, 1, -1
      if (Grid(id)%r_mesh(j) < r) then
         jr = j
         exit LOOP_j
      endif
   enddo LOOP_j
!
   if (jr < Grid(id)%jend) then
      if (Grid(id)%r_mesh(jr+1)-r <= r-Grid(id)%r_mesh(jr)) then
         jr = jr+1
      endif
   endif
!
   end function getRadialGridPoint
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine pushRadialGridToAccel()
!  ===================================================================
   implicit none
!
   integer (kind=IntKind) :: n, id, jend_plus_n
   integer (kind=IntKind) :: iend_packed(NumGrids)
!
   real (kind=RealKind), pointer :: r_mesh(:)
   real (kind=RealKind), allocatable :: hin_packed(:)
   real (kind=RealKind), allocatable :: hout_packed(:)
   real (kind=RealKind), allocatable :: r_mesh_packed(:)
!
   n = 0
   do id = 1, NumGrids
      iend_packed(id) = Grid(id)%jend_plus_n
      n = n + Grid(id)%jend_plus_n
   enddo
   allocate(hin_packed(NumGrids), hout_packed(NumGrids), r_mesh_packed(n))
!
   n = 0
   do id = 1, NumGrids
      jend_plus_n = Grid(id)%jend_plus_n
      r_mesh => Grid(id)%r_mesh
      hin_packed(id) = Grid(id)%hin
      hout_packed(id) = Grid(id)%hout
!     ----------------------------------------------------------------
      call dcopy(jend_plus_n,r_mesh,1,r_mesh_packed(n+1),1)
!     ----------------------------------------------------------------
      n = n + jend_plus_n
   enddo
#ifdef TMP_ACCEL
!  -------------------------------------------------------------------
   call initialize_radialgrid(NumGrids,iend_packed)
   call push_radialgrid(hin_packed, hout_packed, r_mesh_packed)
!  -------------------------------------------------------------------
#endif
   deallocate(hin_packed, hout_packed, r_mesh_packed)
!
   end subroutine pushRadialGridToAccel
!  ===================================================================
!
!  *******************************************************************
!
!  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   subroutine deleteRadialGridOnAccel()
!  ===================================================================
   implicit none
!
#ifdef TMP_ACCEL
!  -------------------------------------------------------------------
   call delete_radialgrid()
!  -------------------------------------------------------------------
#endif
!
   end subroutine deleteRadialGridOnAccel
!  ===================================================================
end module RadialGridModule
