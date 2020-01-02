#include "superlu_ddefs.h"

/* kind of integer to hold a pointer.  Use int.
   This might need to be changed on systems with large memory.
   If changed, be sure to change it in superlutype.f90 too */

typedef int fptr;

/* Used in function f2c_comm */
#define NO_MPI2


/* Fortran name mangling choices */
#define LOWERCASE 1
#define UNDERSCORE 2
#define DOUBLE_UNDERSCORE 3
#define UPPERCASE 4

/* Pick the form of mangling that the Fortran compiler uses.
   This should actually be done on the compile line with -DFNAME=whatever */

#define FNAME LOWERCASE

/* Fortran name mangling */

#if FNAME==LOWERCASE

#define f_create_gridinfo                f_create_gridinfo
#define f_create_options                 f_create_options
#define f_create_ScalePermstruct         f_create_scalepermstruct
#define f_create_LUstruct                f_create_lustruct
#define f_create_SOLVEstruct             f_create_solvestruct
#define f_create_SuperMatrix             f_create_supermatrix
#define f_destroy_gridinfo               f_destroy_gridinfo
#define f_destroy_options                f_destroy_options
#define f_destroy_ScalePermstruct        f_destroy_scalepermstruct
#define f_destroy_LUstruct               f_destroy_lustruct
#define f_destroy_SOLVEstruct            f_destroy_solvestruct
#define f_destroy_SuperMatrix            f_destroy_supermatrix
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
#define f_LUstructInit                   f_lustructinit
#define f_LUstructFree                   f_lustructfree
#define f_Destroy_LU                     f_destroy_lu
#define f_create_SuperLUStat             f_create_superlustat
#define f_destroy_SuperLUStat            f_destroy_superlustat
#define f_dCreate_CompRowLoc_Matrix_dist f_dcreate_comprowloc_mat_dist
#define f_Destroy_CompRowLoc_Matrix_dist f_destroy_comprowloc_mat_dist
#define f_Destroy_SuperMatrix_Store_dist f_destroy_supermat_store_dist
#define f_dSolveFinalize                 f_dsolvefinalize
#define f_pdgssvx                        f_pdgssvx
#define f_dcreate_dist_matrix            f_dcreate_dist_matrix
#define f_check_malloc                   f_check_malloc
#define f_dPrint_CompRowLoc_Matrix_dist  f_dprint_comprowloc_mat_dist

#elif FNAME==UNDERSCORE

#define f_create_gridinfo                f_create_gridinfo_
#define f_create_options                 f_create_options_
#define f_create_ScalePermstruct         f_create_scalepermstruct_
#define f_create_LUstruct                f_create_lustruct_
#define f_create_SOLVEstruct             f_create_solvestruct_
#define f_create_SuperMatrix             f_create_supermatrix_
#define f_destroy_gridinfo               f_destroy_gridinfo_
#define f_destroy_options                f_destroy_options_
#define f_destroy_ScalePermstruct        f_destroy_scalepermstruct_
#define f_destroy_LUstruct               f_destroy_lustruct_
#define f_destroy_SOLVEstruct            f_destroy_solvestruct_
#define f_destroy_SuperMatrix            f_destroy_supermatrix_
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
#define f_LUstructInit                   f_lustructinit_
#define f_LUstructFree                   f_lustructfree_
#define f_Destroy_LU                     f_destroy_lu_
#define f_create_SuperLUStat             f_create_superlustat_
#define f_destroy_SuperLUStat            f_destroy_superlustat_
#define f_dCreate_CompRowLoc_Matrix_dist f_dcreate_comprowloc_mat_dist_
#define f_Destroy_CompRowLoc_Matrix_dist f_destroy_comprowloc_mat_dist_
#define f_Destroy_SuperMatrix_Store_dist f_destroy_supermat_store_dist_
#define f_dSolveFinalize                 f_dsolvefinalize_
#define f_pdgssvx                        f_pdgssvx_
#define f_dcreate_dist_matrix            f_dcreate_dist_matrix_
#define f_check_malloc                   f_check_malloc_
#define f_dPrint_CompRowLoc_Matrix_dist  f_dprint_comprowloc_mat_dist_

#elif FNAME==DOUBLE_UNDERSCORE

#define f_create_gridinfo                f_create_gridinfo__
#define f_create_options                 f_create_options__
#define f_create_ScalePermstruct         f_create_scalepermstruct__
#define f_create_LUstruct                f_create_lustruct__
#define f_create_SOLVEstruct             f_create_solvestruct__
#define f_create_SuperMatrix             f_create_supermatrix__
#define f_destroy_gridinfo               f_destroy_gridinfo__
#define f_destroy_options                f_destroy_options__
#define f_destroy_ScalePermstruct        f_destroy_scalepermstruct__
#define f_destroy_LUstruct               f_destroy_lustruct__
#define f_destroy_SOLVEstruct            f_destroy_solvestruct__
#define f_destroy_SuperMatrix            f_destroy_supermatrix__
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
#define f_LUstructInit                   f_lustructinit__
#define f_LUstructFree                   f_lustructfree__
#define f_Destroy_LU                     f_destroy_lu__
#define f_create_SuperLUStat             f_create_superlustat__
#define f_destroy_SuperLUStat            f_destroy_superlustat__
#define f_dCreate_CompRowLoc_Matrix_dist f_dcreate_comprowloc_mat_dis__
#define f_Destroy_CompRowLoc_Matrix_dist f_destroy_comprowloc_mat_dis__
#define f_Destroy_SuperMatrix_Store_dist f_destroy_supermat_store_dist__
#define f_dSolveFinalize                 f_dsolvefinalize__
#define f_pdgssvx                        f_pdgssvx__
#define f_dcreate_dist_matrix            f_dcreate_dist_matrix__
#define f_check_malloc                   f_check_malloc__
#define f_dPrint_CompRowLoc_Matrix_dist  f_dprint_comprowloc_mat_dist__

#elif FNAME==UPPERCASE

#define f_create_gridinfo                F_CREATE_GRIDINFO
#define f_create_options                 F_CREATE_OPTIONS
#define f_create_ScalePermstruct         F_CREATE_SCALEPERMSTRUCT
#define f_create_LUstruct                F_CREATE_LUSTRUCT
#define f_create_SOLVEstruct             F_CREATE_SOLVESTRUCT
#define f_create_SuperMatrix             F_CREATE_SUPERMATRIX
#define f_destroy_gridinfo               F_DESTROY_GRIDINFO
#define f_destroy_options                F_DESTROY_OPTIONS
#define f_destroy_ScalePermstruct        F_DESTROY_SCALEPERMSTRUCT
#define f_destroy_LUstruct               F_DESTROY_LUSTRUCT
#define f_destroy_SOLVEstruct            F_DESTROY_SOLVESTRUCT
#define f_destroy_SuperMatrix            F_DESTROY_SUPERMATRIX
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
#define f_LUstructInit                   F_LUSTRUCTINIT
#define f_LUstructFree                   F_LUSTRUCTFREE
#define f_Destroy_LU                     F_DESTROY_LU
#define f_create_SuperLUStat             F_CREATE_SUPERLUSTAT
#define f_destroy_SuperLUStat            F_DESTROY_SUPERLUSTAT
#define f_dCreate_CompRowLoc_Matrix_dist F_DCREATE_COMPROWLOC_MAT_DIST
#define f_Destroy_CompRowLoc_Matrix_dist F_DESTROY_COMPROWLOC_MAT_DIST
#define f_Destroy_SuperMatrix_Store_dist F_DESTROY_SUPERMAT_STORE_DIST
#define f_dSolveFinalize                 F_DSOLVEFINALIZE
#define f_pdgssvx                        F_PDGSSVX
#define f_dcreate_dist_matrix            F_DCREATE_DIST_MATRIX
#define f_check_malloc                   F_CHECK_MALLOC
#define f_dPrint_CompRowLoc_Matrix_dist  F_DPRINT_COMPROWLOC_MAT_DIST
#endif

/* some MPI implementations may require conversion between a Fortran
   communicator and a C communicator.  This routine is used to perform the
   conversion.  It may need different forms for different MPI libraries. */

/* NO_MPI2 should be defined on the compiler command line if the MPI
   library does not provide MPI_Comm_f2c */


MPI_Comm f2c_comm(int *f_comm)
{
#ifndef NO_MPI2

/* MPI 2 provides a standard way of doing this */
   return MPI_Comm_f2c((MPI_Fint)(*f_comm));
#else

/* will probably need some special cases here */
/* when in doubt, just return the input */
   return (MPI_Comm)(*f_comm);
#endif
}


/* functions that create memory for a struct and return a handle */

void f_create_gridinfo(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(gridinfo_t));
}

void f_create_options(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(superlu_options_t));
}

void f_create_ScalePermstruct(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(ScalePermstruct_t));
}

void f_create_LUstruct(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(LUstruct_t));
}

void f_create_SOLVEstruct(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(SOLVEstruct_t));
}

void f_create_SuperMatrix(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(SuperMatrix));
}

void f_create_SuperLUStat(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(SuperLUStat_t));
}

/* functions that free the memory allocated by the above functions */

void f_destroy_gridinfo(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_options(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_ScalePermstruct(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_LUstruct(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_SOLVEstruct(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_SuperMatrix(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_SuperLUStat(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

/* functions that get or set values in a C struct.
   This is not the complete set of structs for which a user might want
   to get/set a component, and there may be missing components. */

void f_get_gridinfo(fptr *grid, int *iam, int *nprow, int *npcol)
{
  *iam=((gridinfo_t *) *grid)->iam;
  *npcol=((gridinfo_t *) *grid)->npcol;
  *nprow=((gridinfo_t *) *grid)->nprow;
}

void f_get_SuperMatrix(fptr *A, int *nrow, int *ncol)
{
   *nrow = ((SuperMatrix *) *A)->nrow;
   *ncol = ((SuperMatrix *) *A)->ncol;
}

void f_set_SuperMatrix(fptr *A, int *nrow, int *ncol)
{
   ((SuperMatrix *) *A)->nrow = *nrow;
   ((SuperMatrix *) *A)->ncol = *ncol;
}

void f_get_CompRowLoc_Matrix(fptr *A, int *m, int *n, int *nnz_loc,
                                      int *m_loc, int *fst_row)
{
  *m=((SuperMatrix *) *A)->nrow;
  *n=((SuperMatrix *) *A)->ncol;
  *m_loc=((NRformat_loc *) ((SuperMatrix *) *A)->Store)->m_loc;
  *nnz_loc=((NRformat_loc *) ((SuperMatrix *) *A)->Store)->nnz_loc;
  *fst_row=((NRformat_loc *) ((SuperMatrix *) *A)->Store)->fst_row;
}

void f_set_CompRowLoc_Matrix(fptr *A, int *m, int *n, int *nnz_loc,
                                      int *m_loc, int *fst_row)
{
  ((SuperMatrix *) *A)->nrow = *m;
  ((SuperMatrix *) *A)->ncol = *n;
  ((NRformat_loc *) ((SuperMatrix *) *A)->Store)->m_loc = *m_loc;
  ((NRformat_loc *) ((SuperMatrix *) *A)->Store)->nnz_loc = *nnz_loc;
  ((NRformat_loc *) ((SuperMatrix *) *A)->Store)->fst_row = *fst_row;
}

void f_get_superlu_options(fptr *opt, int *Fact, int *Trans, int *Equil,
                           int *RowPerm, int *ColPerm, int *ReplaceTinyPivot,
                           int *IterRefine, int *SolveInitialized,
                           int *RefineInitialized)
{
   *Fact = (int) ((superlu_options_t *) *opt)->Fact;
   *Trans = (int) ((superlu_options_t *) *opt)->Trans;
   *Equil = (int) ((superlu_options_t *) *opt)->Equil;
   *RowPerm = (int) ((superlu_options_t *) *opt)->RowPerm;
   *ColPerm = (int) ((superlu_options_t *) *opt)->ColPerm;
   *ReplaceTinyPivot = (int) ((superlu_options_t *) *opt)->ReplaceTinyPivot;
   *IterRefine = (int) ((superlu_options_t *) *opt)->IterRefine;
   *SolveInitialized = (int) ((superlu_options_t *) *opt)->SolveInitialized;
   *RefineInitialized = (int) ((superlu_options_t *) *opt)->RefineInitialized;
}

void f_set_superlu_options(fptr *opt, int *Fact, int *Trans, int *Equil,
                           int *RowPerm, int *ColPerm, int *ReplaceTinyPivot,
                           int *IterRefine, int *SolveInitialized,
                           int *RefineInitialized)
{
   ((superlu_options_t *) *opt)->Fact = (fact_t) *Fact;
   ((superlu_options_t *) *opt)->Trans = (trans_t) *Trans;
   ((superlu_options_t *) *opt)->Equil = (yes_no_t) *Equil;
   ((superlu_options_t *) *opt)->RowPerm = (rowperm_t) *RowPerm;
   ((superlu_options_t *) *opt)->ColPerm = (colperm_t) *ColPerm;
   ((superlu_options_t *) *opt)->ReplaceTinyPivot = (yes_no_t) *ReplaceTinyPivot;
   ((superlu_options_t *) *opt)->IterRefine = (IterRefine_t) *IterRefine;
   ((superlu_options_t *) *opt)->SolveInitialized = (yes_no_t) *SolveInitialized;
   ((superlu_options_t *) *opt)->RefineInitialized = (yes_no_t) *RefineInitialized;
}

/* wrappers for SuperLU functions */

void f_set_default_options(fptr *options)
{
   set_default_options_dist((superlu_options_t *) *options);
}

void f_superlu_gridinit(int *Bcomm, int *nprow, int *npcol, fptr *grid)
{
  
   superlu_gridinit(f2c_comm(Bcomm), (int_t) *nprow, (int_t) *npcol,
                    (gridinfo_t *) *grid);
}

void f_superlu_gridexit(fptr *grid)
{
   superlu_gridexit((gridinfo_t *) *grid);
}

void f_ScalePermstructInit(int *m, int *n, fptr *ScalePermstruct)
{
   ScalePermstructInit((const int_t) *m, (const int_t) *n,
                       (ScalePermstruct_t *) *ScalePermstruct);
}

void f_ScalePermstructFree(fptr *ScalePermstruct)
{
   ScalePermstructFree((ScalePermstruct_t *) *ScalePermstruct);
}

void f_PStatInit(fptr *stat)
{
   PStatInit((SuperLUStat_t *) *stat);
}

void f_PStatFree(fptr *stat)
{
   PStatFree((SuperLUStat_t *) *stat);
}

void f_LUstructInit(int *m, int *n, fptr *LUstruct)
{
   LUstructInit((const int_t) *m, (const int_t) *n, (LUstruct_t *) *LUstruct);
}

void f_LUstructFree(fptr *LUstruct)
{
   LUstructFree((LUstruct_t *) *LUstruct);
}

void f_Destroy_LU(int *n, fptr *grid, fptr *LUstruct)
{
   Destroy_LU((int_t) *n, (gridinfo_t *) *grid, (LUstruct_t *) *LUstruct);
}

void f_dCreate_CompRowLoc_Matrix_dist(fptr *A, int *m, int *n, int *nnz_loc,
                                      int *m_loc, int *fst_row, double *nzval,
                                      int *colind, int *rowptr, int *stype,
                                      int *dtype, int *mtype)
{
   dCreate_CompRowLoc_Matrix_dist((SuperMatrix *) *A, (int_t) *m, (int_t) *n,
                                  (int_t) *nnz_loc, (int_t) *m_loc,
                                  (int_t) *fst_row, (double *) nzval,
                                  (int_t *) colind, (int_t *) rowptr,
                                  (Stype_t) *stype, (Dtype_t) *dtype,
                                  (Mtype_t) *mtype);
}

void f_Destroy_CompRowLoc_Matrix_dist(fptr *A)
{
   Destroy_CompRowLoc_Matrix_dist((SuperMatrix *) *A);
}

void f_Destroy_SuperMatrix_Store_dist(fptr *A)
{
   Destroy_SuperMatrix_Store_dist((SuperMatrix *) *A);
}

void f_dSolveFinalize(fptr *options, fptr *SOLVEstruct)
{
   dSolveFinalize((superlu_options_t *) *options,
                  (SOLVEstruct_t *) *SOLVEstruct);
}

void f_pdgssvx(fptr *options, fptr *A, fptr *ScalePermstruct, double *B,
               int *ldb, int *nrhs, fptr *grid, fptr *LUstruct,
               fptr *SOLVEstruct, double *berr, fptr *stat, int *info)
{
   pdgssvx((superlu_options_t *) *options, (SuperMatrix *) *A,
           (ScalePermstruct_t *) *ScalePermstruct, B, *ldb, *nrhs,
           (gridinfo_t *) *grid, (LUstruct_t *) *LUstruct,
           (SOLVEstruct_t *) *SOLVEstruct, berr, (SuperLUStat_t *) *stat, info);
}

/* Create the distributed matrix */

void f_dcreate_dist_matrix(fptr *A, int *m, int *n, int *nnz, double *nzval,
		    int *colind, int *rowptr, fptr *grid, int *m_loc, int *fst_row)
{
   dCreate_CompRowLoc_Matrix_dist((SuperMatrix *) *A, (int_t) *m, (int_t) *n, 
			(int_t) *nnz, (int_t) *m_loc, (int_t) *fst_row,
			(double *) nzval, (int_t *) colind, (int_t *) rowptr,
			SLU_NR_loc, SLU_D, SLU_GE);

}

/* Check the malloc */

void f_check_malloc(int *iam)
{
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC((int_t) *iam, "Check Malloc");
#endif
}

void f_dprint_comprowloc_mat_dist(fptr *A)
{
   dPrint_CompRowLoc_Matrix_dist((SuperMatrix *) *A);
}
