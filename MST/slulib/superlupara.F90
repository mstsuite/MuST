
module superlupara_mod

!----------------------------------------------------
! This module contains some parameter using in SUPERLU for Fortran user.
!----------------------------------------------------


implicit none
public superlu_ptr

!----------------------------------------------------
! kind of integer to hold a SuperLU pointer.  Use default integer.
! This might need to be changed on systems with large memory.
! If changed, be sure to change it in superlu_wrappers.c too.

! integer, parameter :: superlu_ptr = kind(0)
integer, parameter :: superlu_ptr = 8


!----------------------------------------------------
! The following parameters are defined:

! These values come from superlu_defs.h.  If the values in there change with
! the version of SuperLU, then they need to be changed here, too.

integer, parameter, public :: NO                      = 0, & ! yes_no_t
                      YES                     = 1, &
                      DOFACT                  = 0, & ! fact_t
                      SamePattern             = 1, &
                      SamePattern_SameRowPerm = 2, &
                      FACTORED                = 3, &
                      NOROWPERM               = 0, & ! rowperm_t
                      LargeDiag               = 1, &
                      MY_PERMR                = 2, &
                      NATURAL                 = 0, & ! colperm_t
                      MMD_ATA                 = 1, &
                      MMD_AT_PLUS_A           = 2, &
                      COLAMD                  = 3, &
                      MY_PERMC                = 4, &
                      NOTRANS                 = 0, & ! trans_t
                      TRANS                   = 1, &
                      CONJ                    = 2, &
                      NOEQUIL                 = 0, & ! DiagScale_t  Need?
                      ROW                     = 1, &
                      COL                     = 2, &
                      BOTH                    = 3, &
                      NOREFINE                = 0, & ! IterRefine_t
                      DOUBLE                  = 1, &
                      EXTRA                   = 2, &
                      LUSUP                   = 0, & ! MemType  Need?
                      UCOL                    = 1, &
                      LSUB                    = 2, &
                      USUB                    = 3, &
                      SYSTEM                  = 0, & ! LU_space_t  Need?
                      USER                    = 1
integer, parameter, public :: SLU_NC                  = 0, & ! Stype_t
                      SLU_NR                  = 1, &
                      SLU_SC                  = 2, &
                      SLU_SR                  = 3, &
                      SLU_NCP                 = 4, &
                      SLU_DN                  = 5, &
                      SLU_NR_loc              = 6, &
                      SLU_S                   = 0, & ! Dtype_t
                      SLU_D                   = 1, &
                      SLU_C                   = 2, &
                      SLU_Z                   = 3, &
                      SLU_GE                  = 0, & ! Mtype_t
                      SLU_TRLU                = 1, &
                      SLU_TRUU                = 2, &
                      SLU_TRL                 = 3, &
                      SLU_TRU                 = 4, &
                      SLU_SYL                 = 5, &
                      SLU_SYU                 = 6, &
                      SLU_HEL                 = 7, &
                      SLU_HEU                 = 8


!----------------------------------------------------

end module superlupara_mod
