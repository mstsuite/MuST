module nvtx
  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter :: NVTX_COLOR_GREEN = 2
  integer, parameter :: NVTX_COLOR_BLUE = 3
  integer, parameter :: NVTX_COLOR_RED = 1

  interface
    subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
      use iso_c_binding
      character(kind=C_CHAR) :: name(*)
    end subroutine nvtxRangePushA
    
    subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
    end subroutine nvtxRangePop
    
    subroutine nvtxMarkA(name) bind(C, name='nvtxMarkA')
      use iso_c_binding
      character(kind=C_CHAR) :: name(*)
    end subroutine nvtxMarkA
  end interface
end module nvtx
