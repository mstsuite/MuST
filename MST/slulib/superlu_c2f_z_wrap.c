#include "superlu_zdefs.h"
#include "superlu_c2f_z_wrap.h"

/* kind of integer to hold a pointer.  Use int.
   This might need to be changed on systems with large memory.
   If changed, be sure to change it in superlupara.f90 too */

typedef long int fptr;  /* 64-bit */

/* some MPI implementations may require conversion between a Fortran
   communicator and a C communicator.  This routine is used to perform the
   conversion.  It may need different forms for different MPI libraries. */

/* NO_MPI2 should be defined on the compiler command line if the MPI
   library does not provide MPI_Comm_f2c */

MPI_Comm f2c_comm(int *f_comm)
{
#ifdef NO_MPI2

/* MPI 2 provides a standard way of doing this */
   return MPI_Comm_f2c((MPI_Fint)(*f_comm));
#else

/* will probably need some special cases here */
/* when in doubt, just return the input */
   return (MPI_Comm)(*f_comm);
#endif
}


/* functions that create memory for a struct and return a handle */

void f_create_gridinfo_handle(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(gridinfo_t));
}

void f_create_options_handle(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(superlu_options_t));
}

void f_create_ScalePermstruct_handle(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(ScalePermstruct_t));
}

void f_create_LUstruct_handle(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(LUstruct_t));
}

void f_create_SOLVEstruct_handle(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(SOLVEstruct_t));
}

void f_create_SuperMatrix_handle(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(SuperMatrix));
}

void f_create_SuperLUStat_handle(fptr *handle)
{
   *handle = (fptr) SUPERLU_MALLOC(sizeof(SuperLUStat_t));
}

/* functions that free the memory allocated by the above functions */

void f_destroy_gridinfo_handle(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_options_handle(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_ScalePermstruct_handle(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_LUstruct_handle(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_SOLVEstruct_handle(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_SuperMatrix_handle(fptr *handle)
{
   SUPERLU_FREE((void *)*handle);
}

void f_destroy_SuperLUStat_handle(fptr *handle)
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
                           int *RefineInitialized, int *PrintStat)
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
   *PrintStat = (int) ((superlu_options_t *) *opt)->PrintStat;
}

void f_set_superlu_options(fptr *opt, int *Fact, int *Trans, int *Equil,
                           int *RowPerm, int *ColPerm, int *ReplaceTinyPivot,
                           int *IterRefine, int *SolveInitialized,
                           int *RefineInitialized, int *PrintStat)
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
   ((superlu_options_t *) *opt)->PrintStat = (yes_no_t) *PrintStat;
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

void f_zCreate_CompRowLoc_Matrix_dist(fptr *A, int *m, int *n, int *nnz_loc,
                                      int *m_loc, int *fst_row, doublecomplex *nzval,
                                      int *colind, int *rowptr, int *stype,
                                      int *dtype, int *mtype)
{
   zCreate_CompRowLoc_Matrix_dist((SuperMatrix *) *A, (int_t) *m, (int_t) *n,
                                  (int_t) *nnz_loc, (int_t) *m_loc,
                                  (int_t) *fst_row, (doublecomplex *) nzval,
                                  (int_t *) colind, (int_t *) rowptr,
                                  (Stype_t) *stype, (Dtype_t) *dtype,
                                  (Mtype_t) *mtype);
}

void f_Destroy_CompRowLoc_Matrix_dist(fptr *A)
{
/*   Destroy_CompRowLoc_Matrix_dist_test((SuperMatrix *) *A); */
   Destroy_CompRowLoc_Matrix_dist((SuperMatrix *) *A);
}

void f_Destroy_SuperMatrix_Store_dist(fptr *A)
{
   Destroy_SuperMatrix_Store_dist((SuperMatrix *) *A);
}

void f_zSolveFinalize(fptr *options, fptr *SOLVEstruct)
{
   zSolveFinalize((superlu_options_t *) *options,
                  (SOLVEstruct_t *) *SOLVEstruct);
}

void f_pzgssvx(fptr *options, fptr *A, fptr *ScalePermstruct, doublecomplex *B,
               int *ldb, int *nrhs, fptr *grid, fptr *LUstruct,
               fptr *SOLVEstruct, double *berr, fptr *stat, int *info)
{
   pzgssvx((superlu_options_t *) *options, (SuperMatrix *) *A,
           (ScalePermstruct_t *) *ScalePermstruct, B, *ldb, *nrhs,
           (gridinfo_t *) *grid, (LUstruct_t *) *LUstruct,
           (SOLVEstruct_t *) *SOLVEstruct, berr, (SuperLUStat_t *) *stat, info);
/* PStatPrint((superlu_options_t *) *options, (SuperLUStat_t *) *stat,
	      (gridinfo_t *) *grid); */
}

/* Create the distributed matrix */

void f_zcreate_dist_matrix(fptr *A, int *m, int *n, int *nnz, double *nzval,
		    int *colind, int *rowptr, fptr *grid, int *m_loc, int *fst_row)
{
   zCreate_CompRowLoc_Matrix_dist((SuperMatrix *) *A, (int_t) *m, (int_t) *n, 
			(int_t) *nnz, (int_t) *m_loc, (int_t) *fst_row,
			(doublecomplex *) nzval, (int_t *) colind, (int_t *) rowptr,
			SLU_NR_loc, SLU_Z, SLU_GE);
}

/* Check the malloc */

void f_check_malloc(int *iam)
{
#if ( DEBUGlevel>=1 )
    CHECK_MALLOC((int_t) *iam, "Check Malloc");
#endif
}
