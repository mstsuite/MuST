module PublicTypeDefinitionsModule
   use KindParamModule, only : IntKind, RealKind, CmplxKind
   implicit none
!
public
!
!  AtomModule data types
!
   type LizLmaxStruct
      real (kind=RealKind) :: rad
      real (kind=RealKind) :: rad_s
      integer (kind=IntKind) :: nmax
      integer (kind=IntKind) :: NumShells
      integer (kind=IntKind), pointer :: lmax_shell(:)
   end type LizLmaxStruct
!
!  InputModule
!
   type InputTableStruct
      character (len=160) :: TableName
      integer (kind=IntKind) :: TableIndex
      integer (kind=IntKind) :: NumData
      integer (kind=IntKind) :: UnitNumber
      logical :: Open
      character (len=50), pointer :: KeyName(:)
      character (len=160), pointer :: KeyValue(:)
      type (InputTableStruct), pointer :: next
   end type InputTableStruct
!
!  MixingModule
!
   type MixListStruct
      integer (kind=IntKind) :: size
      real (kind=RealKind) :: rms
      real (kind=RealKind) :: weight
      real (kind=RealKind), pointer :: mesh(:)
      real (kind=RealKind), pointer :: vector_old(:)
      real (kind=RealKind), pointer :: vector_new(:)
      type (MixListStruct), pointer :: next
   end type MixListStruct
!
!  NeighborModule
!
   type NeighborStruct
      integer (kind=IntKind) :: NumAtoms             ! no. of neighbor atoms
                                                     ! in cluster. Note: center
                                                     ! atom is not counted.
      integer (kind=IntKind) :: NumReceives          ! No. receiving calls
                                                     ! needed for retrieving
                                                     ! the neighboring atoms'
                                                     ! data
      integer (kind=IntKind), pointer :: Z(:)        ! Atomic number of 
                                                     ! neigbor atom
      integer (kind=IntKind), pointer :: Lmax(:)     ! L-cut off for each atom
                                                     ! in the cluster
      integer (kind=IntKind), pointer :: ProcIndex(:)  ! the index of the CPU
                                                     ! that the neighbor atom
                                                     ! is mapped onto.
      integer (kind=IntKind), pointer :: LocalIndex(:) ! the local index of the
                                                     ! neighbor atom on the CPU
                                                     ! it is mapped onto.
      integer (kind=IntKind), pointer :: GlobalIndex(:) ! the global index of
                                                     ! the neighbor atom.
      real (kind=RealKind), pointer :: Position(:,:) ! the coordinates of the
                                                     ! neighbor atom.
      integer (kind=IntKind), pointer :: IndexMN(:)  ! the count of contributing
                                                     ! atoms to Tau
      real (kind=RealKind), pointer :: Rmunu(:,:,:)  !

      integer (kind=IntKind) :: NumShells            ! no. of shells in LIZ
      integer (kind=IntKind),pointer:: ShellIndex(:) ! index to which shell 
                                                     ! this neighbor atom
                                                     ! belongs
      integer (kind=IntKind),pointer:: NAsOnShell(:) ! Number Atoms on each shell
      real (kind=RealKind), pointer :: ShellRad(:)   ! radii of shells

   end type NeighborStruct
!
   type CommTableStruct
      integer (kind=IntKind) :: NumRequests          ! No. of times that my atom is on other
                                                     ! atoms's neighbor list or no. of times
                                                     ! that other atoms are on my neighbor list
      integer (kind=IntKind) :: NumComms             ! No. of sending/receiving calls needed
      integer (kind=IntKind), allocatable :: RemoteAtomList(:) ! The list of atoms that need
                                                               ! my t-matrix or whose t-matrix are needed
      integer (kind=IntKind), allocatable :: RemoteProcList(:) ! The list of processors that
                                                               ! need me to send my t-matrix or
                                                               ! that will send their t-matrix to me
   end type CommTableStruct
!
!  RadialGridModule
!
   type GridStruct
!     integer (kind=IntKind) :: jws         ! the r-mesh point closest to rws
!     real (kind=RealKind) :: rws           ! the Wigner-Seitz sphere radius
!     real (kind=RealKind) :: rinsc         ! the inscribed sphere radius 
!     real (kind=RealKind) :: xinsc         ! log(rinsc)
      integer (kind=IntKind) :: jinsc       ! Num. of r-mesh inside the inscribed
                                            ! sphere radius.
      real (kind=RealKind) :: rend          ! the terminal sphere radius
                                            ! = rmt, rws, or rcirc
      real (kind=RealKind) :: rmt           ! the muffin-tin sphere radius
      real (kind=RealKind) :: rstart        ! the near origin starting point
      real (kind=RealKind) :: xend          ! log(rend)
      real (kind=RealKind) :: xmt           ! log(rmt)
      real (kind=RealKind) :: xstart        ! the near origin starting point
      real (kind=RealKind) :: hout          ! x-step bwteen xmt and xcirc
      real (kind=RealKind) :: hin           ! x-step bwteen xstart and xmt
!
      integer (kind=IntKind) :: nmult       ! = hin/hout
      integer (kind=IntKind) :: jend        ! the number of x or r mesh points
                                            ! inside rend (including rend)
      integer (kind=IntKind) :: jmt         ! Num. of r-mesh on and inside
                                            ! the Muffin-tin radius.
      integer (kind=IntKind) :: jend_plus_n ! = jend+n_extra
!
      real (kind=RealKind), pointer :: x_mesh(:)
      real (kind=RealKind), pointer :: r_mesh(:)
   end type GridStruct
!
!  StrConstModule 
!
   type ScmBlockStruct
      integer (kind=IntKind) :: lmaxi
      integer (kind=IntKind) :: lmaxj
      complex (kind=CmplxKind), pointer :: strcon_matrix(:,:)
   end type ScmBlockStruct
!
   type AtomOnUniformGridStruct
      integer (kind=IntKind) :: NumLocalAtoms               ! Determined by the parallelization over
                                                            ! atoms
      real (kind=RealKind), allocatable :: AtomPosition(:,:)
      integer (kind=IntKind), allocatable :: AtomBox(:,:,:) ! Stores the starting and ending index
                                                            ! of the uniform grid points associated
                                                            ! with each atom
      integer (kind=IntKind), allocatable :: NumGridPointsInAtomBox(:)
      integer (kind=IntKind), allocatable :: NumGridPointsInCell(:) ! Number of grid points in atomic cell
      integer (kind=IntKind), allocatable :: InCellGridPointABIndex(:,:) ! AtomBox index of each grid point in atomic cell
      integer (kind=IntKind) :: NumGridPointsOnCellBound    ! number of grid points on cell boundary for all local atoms
      integer (kind=IntKind), allocatable :: CBGridPointIndex(:) ! Global index of each grid point on cell boundary
      integer (kind=IntKind), allocatable :: NumNeighboringAtoms(:) ! number of associated atoms of each grid point on cell boundary
      integer (kind=IntKind), allocatable :: NumTargetProcs(:)
      integer (kind=IntKind), allocatable :: NumSourceProcs(:)
      integer (kind=IntKind), allocatable :: TargetProc(:,:)
      integer (kind=IntKind), allocatable :: SourceProc(:,:)
   end type AtomOnUniformGridStruct
!
   type UniformGridStruct
      integer (kind=IntKind) :: nga    ! no. of grid points along a-dim
      integer (kind=IntKind) :: ngb    ! no. of grid points along b-dim
      integer (kind=IntKind) :: ngc    ! no. of grid points along c-dim
      integer (kind=IntKind) :: ng     ! ng = nga*ngb*ngc
      integer (kind=IntKind) :: nproc_a ! no. of processor mesh along a-dim
      integer (kind=IntKind) :: nproc_b ! no. of processor mesh along b-dim
      integer (kind=IntKind) :: nproc_c ! no. of processor mesh along c-dim
      integer (kind=IntKind) :: gstart(3)
      integer (kind=IntKind) :: gend(3)
      integer (kind=IntKind) :: NumLocalGridPoints
      real (kind=RealKind) :: vec_origin(3)
      real (kind=RealKind) :: cell(3,3)
      real (kind=RealKind) :: grid_step_a(3)
      real (kind=RealKind) :: grid_step_b(3)
      real (kind=RealKind) :: grid_step_c(3)
      type (AtomOnUniformGridStruct), pointer :: AtomOnGrid
   end type UniformGridStruct
!
   type LloydStruct
      logical:: lloydOn
      logical:: lloydOnly
      complex(kind=CmplxKind) :: q_lloyd(2)
   end type LloydStruct
!
   type SROStruct
      complex (kind=CmplxKind), pointer :: tau_l_sro(:,:,:) 
      complex (kind=CmplxKind), pointer :: kau_l_sro(:,:,:)
   end type SROStruct
!
   type MatrixBlockStruct
      integer (kind=IntKind) :: lmax_kkr
      integer (kind=IntKind) :: kmax_kkr
      integer (kind=IntKind) :: row_index
      integer (kind=IntKind) :: global_index
      complex (kind=CmplxKind), pointer :: tau_l(:,:,:)
      complex (kind=CmplxKind), pointer :: kau_l(:,:,:)
      type (SROStruct), pointer :: NeighMat_l(:) ! Contains all the matrices for the neighbors
   end type MatrixBlockStruct
!
   type MatrixBandStruct
      integer (kind=IntKind) :: global_index
      integer (kind=IntKind) :: column_index  ! The column index in the band matrix, not the "big" KKR matrix
      integer (kind=IntKind) :: lmax_kkr
      integer (kind=IntKind) :: kmax_kkr
      type (MatrixBlockStruct), pointer :: MatrixBlock(:)
   end type MatrixBandStruct
!
end module PublicTypeDefinitionsModule
