/*
 * -- Distributed SuperLU header file (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 */

#ifndef __SUPERLU_C2F_WRAP /* allow multiple inclusions */
#define __SUPERLU_C2F_WRAP


/* Used in function f2c_comm */
/* #define NO_MPI2 */

/* Fortran name mangling choices */
#define LOWERCASE 1
#define UNDERSCORE 2
#define DOUBLE_UNDERSCORE 3
#define UPPERCASE 4

/* Pick the form of mangling that the Fortran compiler uses.
   This should actually be done on the compile line with -DFNAME=whatever */

#ifdef NoChange
#define FNAME LOWERCASE 
#elif Add_
#define FNAME UNDERSCORE 
#elif Add__
#define FNAME DOUBLE_UNDERSCORE 
#elif UpCase
#define FNAME UPPERCASE 
#endif

/* Fortran name mangling */

#if FNAME==LOWERCASE

#define f_create_gridinfo_handle         f_create_gridinfo_handle
#define f_create_options_handle          f_create_options_handle
#define f_create_ScalePermstruct_handle  f_create_scalepermstruct_handle
#define f_create_LUstruct_handle         f_create_lustruct_handle
#define f_create_SOLVEstruct_handle      f_create_solvestruct_handle
#define f_create_SuperMatrix_handle      f_create_supermatrix_handle
#define f_destroy_gridinfo_handle        f_destroy_gridinfo_handle
#define f_destroy_options_handle         f_destroy_options_handle
#define f_destroy_ScalePermstruct_handle f_destroy_scalepermstruct_handle
#define f_destroy_LUstruct_handle        f_destroy_lustruct_handle
#define f_destroy_SOLVEstruct_handle     f_destroy_solvestruct_handle
#define f_destroy_SuperMatrix_handle     f_destroy_supermatrix_handle
#define f_create_SuperLUStat_handle      f_create_superlustat_handle
#define f_destroy_SuperLUStat_handle     f_destroy_superlustat_handle
#define f_get_gridinfo                   f_get_gridinfo
#define f_get_SuperMatrix                f_get_supermatrix
#define f_set_SuperMatrix                f_set_supermatrix
#define f_get_CompRowLoc_Matrix          f_get_comprowloc_matrix 
#define f_set_CompRowLoc_Matrix          f_set_comprowloc_matrix
#define f_get_superlu_options            f_get_superlu_options
#define f_set_superlu_options            f_set_superlu_options
#define f_set_default_options            f_set_default_options
#define f_superlu_gridinit               f_superlu_gridinit
#define f_superlu_gridexit               f_superlu_gridexit
#define f_ScalePermstructInit            f_scalepermstructinit
#define f_ScalePermstructFree            f_scalepermstructfree
#define f_PStatInit                      f_pstatinit
#define f_PStatFree                      f_pstatfree
#define f_PStatPrint                     f_pstatprint
#define f_LUstructInit                   f_lustructinit
#define f_LUstructFree                   f_lustructfree
#define f_Destroy_LU                     f_destroy_lu
#define f_create_SuperLUStat             f_create_superlustat
#define f_destroy_SuperLUStat            f_destroy_superlustat
#define f_zCreate_CompRowLoc_Matrix_dist f_zcreate_comprowloc_matrix_dist
#define f_Destroy_CompRowLoc_Matrix_dist f_destroy_comprowloc_matrix_dist
#define f_Destroy_SuperMatrix_Store_dist f_destroy_supermatrix_store_dist
#define f_zSolveFinalize                 f_zsolvefinalize
#define f_pzgssvx                        f_pzgssvx
#define f_zcreate_dist_matrix            f_zcreate_dist_matrix
#define f_check_malloc                   f_check_malloc
/* #define f_zPrint_CompRowLoc_Matrix_dist  f_zprint_comprowloc_matrix_dist */

#elif FNAME==UNDERSCORE

#define f_create_gridinfo_handle         f_create_gridinfo_handle_
#define f_create_options_handle          f_create_options_handle_
#define f_create_ScalePermstruct_handle  f_create_scalepermstruct_handle_
#define f_create_LUstruct_handle         f_create_lustruct_handle_
#define f_create_SOLVEstruct_handle      f_create_solvestruct_handle_
#define f_create_SuperMatrix_handle      f_create_supermatrix_handle_
#define f_destroy_gridinfo_handle        f_destroy_gridinfo_handle_
#define f_destroy_options_handle         f_destroy_options_handle_
#define f_destroy_ScalePermstruct_handle f_destroy_scalepermstruct_handle_
#define f_destroy_LUstruct_handle        f_destroy_lustruct_handle_
#define f_destroy_SOLVEstruct_handle     f_destroy_solvestruct_handle_
#define f_destroy_SuperMatrix_handle     f_destroy_supermatrix_handle_
#define f_create_SuperLUStat_handle      f_create_superlustat_handle_
#define f_destroy_SuperLUStat_handle     f_destroy_superlustat_handle_
#define f_get_gridinfo                   f_get_gridinfo_
#define f_get_SuperMatrix                f_get_supermatrix_
#define f_set_SuperMatrix                f_set_supermatrix_
#define f_get_CompRowLoc_Matrix          f_get_comprowloc_matrix_ 
#define f_set_CompRowLoc_Matrix          f_set_comprowloc_matrix_
#define f_get_superlu_options            f_get_superlu_options_
#define f_set_superlu_options            f_set_superlu_options_
#define f_set_default_options            f_set_default_options_
#define f_superlu_gridinit               f_superlu_gridinit_
#define f_superlu_gridexit               f_superlu_gridexit_
#define f_ScalePermstructInit            f_scalepermstructinit_
#define f_ScalePermstructFree            f_scalepermstructfree_
#define f_PStatInit                      f_pstatinit_
#define f_PStatFree                      f_pstatfree_
#define f_PStatPrint                     f_pstatprint_
#define f_LUstructInit                   f_lustructinit_
#define f_LUstructFree                   f_lustructfree_
#define f_Destroy_LU                     f_destroy_lu_
#define f_create_SuperLUStat             f_create_superlustat_
#define f_destroy_SuperLUStat            f_destroy_superlustat_
#define f_zCreate_CompRowLoc_Matrix_dist f_zcreate_comprowloc_matrix_dist_
#define f_Destroy_CompRowLoc_Matrix_dist f_destroy_comprowloc_matrix_dist_
#define f_Destroy_SuperMatrix_Store_dist f_destroy_supermatrix_store_dist_
#define f_zSolveFinalize                 f_zsolvefinalize_
#define f_pzgssvx                        f_pzgssvx_
#define f_zcreate_dist_matrix            f_zcreate_dist_matrix_
#define f_check_malloc                   f_check_malloc_
/* #define f_zPrint_CompRowLoc_Matrix_dist  f_zprint_comprowloc_matrix_dist_ */

#elif FNAME==DOUBLE_UNDERSCORE

#define f_create_gridinfo_handle         f_create_gridinfo_handle__
#define f_create_options_handle          f_create_options_handle__
#define f_create_ScalePermstruct_handle  f_create_scalepermstruct_handle__
#define f_create_LUstruct_handle         f_create_lustruct_handle__
#define f_create_SOLVEstruct_handle      f_create_solvestruct_handle__
#define f_create_SuperMatrix_handle      f_create_supermatrix_handle__
#define f_destroy_gridinfo_handle        f_destroy_gridinfo_handle__
#define f_destroy_options_handle         f_destroy_options_handle__
#define f_destroy_ScalePermstruct_handle f_destroy_scalepermstruct_handle__
#define f_destroy_LUstruct_handle        f_destroy_lustruct_handle__
#define f_destroy_SOLVEstruct_handle     f_destroy_solvestruct_handle__
#define f_destroy_SuperMatrix_handle     f_destroy_supermatrix_handle__
#define f_create_SuperLUStat_handle      f_create_superlustat_handle__
#define f_destroy_SuperLUStat_handle     f_destroy_superlustat_handle__
#define f_get_gridinfo                   f_get_gridinfo__
#define f_get_SuperMatrix                f_get_supermatrix__
#define f_set_SuperMatrix                f_set_supermatrix__
#define f_get_SuperMatrix                f_get_supermatrix__
#define f_set_SuperMatrix                f_set_supermatrix__
#define f_get_superlu_options            f_get_superlu_options__
#define f_set_superlu_options            f_set_superlu_options__
#define f_set_default_options            f_set_default_options__
#define f_superlu_gridinit               f_superlu_gridinit__
#define f_superlu_gridexit               f_superlu_gridexit__
#define f_ScalePermstructInit            f_scalepermstructinit__
#define f_ScalePermstructFree            f_scalepermstructfree__
#define f_PStatInit                      f_pstatinit__
#define f_PStatFree                      f_pstatfree__
#define f_PStatPrint                     f_pstatprint__
#define f_LUstructInit                   f_lustructinit__
#define f_LUstructFree                   f_lustructfree__
#define f_Destroy_LU                     f_destroy_lu__
#define f_create_SuperLUStat             f_create_superlustat__
#define f_destroy_SuperLUStat            f_destroy_superlustat__
#define f_zCreate_CompRowLoc_Matrix_dist f_zcreate_comprowloc_matrix_dist__
#define f_Destroy_CompRowLoc_Matrix_dist f_destroy_comprowloc_matrix_dist__
#define f_Destroy_SuperMatrix_Store_dist f_destroy_supermatrix_store_dist__
#define f_zSolveFinalize                 f_zsolvefinalize__
#define f_pzgssvx                        f_pzgssvx__
#define f_zcreate_dist_matrix            f_zcreate_dist_matrix__
#define f_check_malloc                   f_check_malloc__
/* #define f_zPrint_CompRowLoc_Matrix_dist  f_zprint_comprowloc_matrix_dist__ */

#elif FNAME==UPPERCASE

#define f_create_gridinfo_handle         F_CREATE_GRIDINFO_HANDLE
#define f_create_options_handle          F_CREATE_OPTIONS_HANDLE
#define f_create_ScalePermstruct_handle  F_CREATE_SCALEPERMSTRUCT_HANDLE
#define f_create_LUstruct_handle         F_CREATE_LUSTRUCT_HANDLE
#define f_create_SOLVEstruct_handle      F_CREATE_SOLVESTRUCT_HANDLE
#define f_create_SuperMatrix_handle      F_CREATE_SUPERMATRIX_HANDLE
#define f_destroy_gridinfo_handle        F_DESTROY_GRIDINFO_HANDLE
#define f_destroy_options_handle         F_DESTROY_OPTIONS_HANDLE
#define f_destroy_ScalePermstruct_handle F_DESTROY_SCALEPERMSTRUCT_HANDLE
#define f_destroy_LUstruct_handle        F_DESTROY_LUSTRUCT_HANDLE
#define f_destroy_SOLVEstruct_handle     F_DESTROY_SOLVESTRUCT_HANDLE
#define f_destroy_SuperMatrix_handle     F_DESTROY_SUPERMATRIX_HANDLE
#define f_create_SuperLUStat_handle      F_CREATE_SUPERLUSTAT_HANDLE
#define f_destroy_SuperLUStat_handle     F_DESTROY_SUPERLUSTAT_HANDLE
#define f_get_gridinfo                   F_GET_GRIDINFO
#define f_get_SuperMatrix                F_GET_SUPERMATRIX
#define f_set_SuperMatrix                F_SET_SUPERMATRIX
#define f_get_CompRowLoc_Matrix          F_GET_COMPROWLOC_MATRIX
#define f_set_CompRowLoc_Matrix          F_SET_COMPROWLOC_MATRIX
#define f_get_superlu_options            F_GET_SUPERLU_OPTIONS
#define f_set_superlu_options            F_SET_SUPERLU_OPTIONS
#define f_set_default_options            F_SET_DEFAULT_OPTIONS
#define f_superlu_gridinit               F_SUPERLU_GRIDINIT
#define f_superlu_gridexit               F_SUPERLU_GRIDEXIT
#define f_ScalePermstructInit            F_SCALEPERMSTRUCTINIT
#define f_ScalePermstructFree            F_SCALEPERMSTRUCTFREE
#define f_PStatInit                      F_PSTATINIT
#define f_PStatFree                      F_PSTATFREE
#define f_PStatPrint                     F_PSTATPRINT
#define f_LUstructInit                   F_LUSTRUCTINIT
#define f_LUstructFree                   F_LUSTRUCTFREE
#define f_Destroy_LU                     F_DESTROY_LU
#define f_create_SuperLUStat             F_CREATE_SUPERLUSTAT
#define f_destroy_SuperLUStat            F_DESTROY_SUPERLUSTAT
#define f_zCreate_CompRowLoc_Matrix_dist F_ZCREATE_COMPROWLOC_MATRIX_DIST
#define f_Destroy_CompRowLoc_Matrix_dist F_DESTROY_COMPROWLOC_MATRIX_DIST
#define f_Destroy_SuperMatrix_Store_dist F_DESTROY_SUPERMATRIX_STORE_DIST
#define f_zSolveFinalize                 F_ZSOLVEFINALIZE
#define f_pzgssvx                        F_PZGSSVX
#define f_zcreate_dist_matrix            F_ZCREATE_DIST_MATRIX
#define f_check_malloc                   F_CHECK_MALLOC
/* #define f_zPrint_CompRowLoc_Matrix_dist  F_ZPRINT_COMPROWLOC_MATRIX_DIST */
#endif

#endif /* __SUPERLU_C2F_WRAP */
