#ifndef HEADER_lp_matrix
#define HEADER_lp_matrix

#include "lp_types.h"


/* Sparse matrix element (ordered columnwise) */
typedef struct _MATitem
{
  int  rownr;
  int  colnr;
  REAL value;
} MATitem;


/* Compiler option development features */
/*#define DebugInv*/                /* Report array values at inversion */


/* Matrix column access macros to be able to easily change storage model */
#define CAM_Record                0
#define CAM_Vector                1
#if 0
 #define MatrixColAccess           CAM_Record
#else
 #define MatrixColAccess           CAM_Vector
#endif

#if MatrixColAccess==CAM_Record
#define SET_MAT_ijA(item,i,j,A)   mat->col_mat[item].rownr = i; \
                                  mat->col_mat[item].colnr = j; \
                                  mat->col_mat[item].value = A
#define COL_MAT_COLNR(item)       (mat->col_mat[item].colnr)
#define COL_MAT_ROWNR(item)       (mat->col_mat[item].rownr)
#define COL_MAT_VALUE(item)       (mat->col_mat[item].value)
#define COL_MAT_COPY(left,right)  mat->col_mat[left] = mat->col_mat[right]
#define COL_MAT_MOVE(to,from,rec) MEMMOVE(&(mat->col_mat[to]),&(mat->col_mat[from]),rec)
#define matRowColStep             (sizeof(MATitem)/sizeof(int))
#define matValueStep              (sizeof(MATitem)/sizeof(REAL))

#else /* if MatrixColAccess==CAM_Vector */
#define SET_MAT_ijA(item,i,j,A)   mat->col_mat_rownr[item] = i; \
                                  mat->col_mat_colnr[item] = j; \
                                  mat->col_mat_value[item] = A
#define COL_MAT_COLNR(item)       (mat->col_mat_colnr[item])
#define COL_MAT_ROWNR(item)       (mat->col_mat_rownr[item])
#define COL_MAT_VALUE(item)       (mat->col_mat_value[item])
#define COL_MAT_COPY(left,right)  COL_MAT_COLNR(left) = COL_MAT_COLNR(right); \
                                  COL_MAT_ROWNR(left) = COL_MAT_ROWNR(right); \
                                  COL_MAT_VALUE(left) = COL_MAT_VALUE(right)
#define COL_MAT_MOVE(to,from,rec) MEMMOVE(&COL_MAT_COLNR(to),&COL_MAT_COLNR(from),rec); \
                                  MEMMOVE(&COL_MAT_ROWNR(to),&COL_MAT_ROWNR(from),rec); \
                                  MEMMOVE(&COL_MAT_VALUE(to),&COL_MAT_VALUE(from),rec)
#define matRowColStep             1
#define matValueStep              1

#endif


/* Matrix row access macros to be able to easily change storage model */
#define RAM_Index                 0
#define RAM_FullCopy              1
#define MatrixRowAccess           RAM_Index

#if MatrixRowAccess==RAM_Index
#define ROW_MAT_COLNR(item)       COL_MAT_COLNR(mat->row_mat[item])
#define ROW_MAT_ROWNR(item)       COL_MAT_ROWNR(mat->row_mat[item])
#define ROW_MAT_VALUE(item)       COL_MAT_VALUE(mat->row_mat[item])

#elif MatrixColAccess==CAM_Record
#define ROW_MAT_COLNR(item)       (mat->row_mat[item].colnr)
#define ROW_MAT_ROWNR(item)       (mat->row_mat[item].rownr)
#define ROW_MAT_VALUE(item)       (mat->row_mat[item].value)

#else /* if MatrixColAccess==CAM_Vector */
#define ROW_MAT_COLNR(item)       (mat->row_mat_colnr[item])
#define ROW_MAT_ROWNR(item)       (mat->row_mat_rownr[item])
#define ROW_MAT_VALUE(item)       (mat->row_mat_value[item])

#endif


typedef struct _MATrec
{
  /* Owner reference */
  lprec     *lp;

  /* Active dimensions */
  int       rows;
  int       columns;
  REAL      epsvalue;           /* Zero element rejection threshold */

  /* Allocated memory */
  int       rows_alloc;
  int       columns_alloc;
  int       mat_alloc;          /* The allocated size for matrix sized structures */

  /* Sparse problem matrix storage */
#if MatrixColAccess==CAM_Record  
  MATitem   *col_mat;           /* mat_alloc : The sparse data storage */
#else /*MatrixColAccess==CAM_Vector*/
  int       *col_mat_colnr;
  int       *col_mat_rownr;
  REAL      *col_mat_value;
#endif  
  int       *col_end;           /* columns_alloc+1 : col_end[i] is the index of the
                                   first element after column i; column[i] is stored
                                   in elements col_end[i-1] to col_end[i]-1 */
#if MatrixRowAccess==RAM_Index
  int       *row_mat;           /* mat_alloc : From index 0, row_mat contains the
                                   row-ordered index of the elements of col_mat */
#elif MatrixColAccess==CAM_Record
  MATitem   *row_mat;           /* mat_alloc : From index 0, row_mat contains the
                                   row-ordered copy of the elements in col_mat */
#else /*if MatrixColAccess==CAM_Vector*/
  int       *row_mat_colnr;
  int       *row_mat_rownr;
  REAL      *row_mat_value;
#endif
  int       *row_end;           /* rows_alloc+1 : row_end[i] is the index of the
                                   first element in row_mat after row i */
  MYBOOL    row_end_valid;      /* TRUE if row_end & row_mat are valid */
  MYBOOL    is_roworder;        /* TRUE if the current (temporary) matrix order is row-wise */

} MATrec;


#ifdef __cplusplus
__EXTERN_C {
#endif

/* Sparse matrix routines */
STATIC MATrec *mat_create(lprec *lp, int rows, int columns, REAL epsvalue);
STATIC void mat_free(MATrec **matrix);
STATIC MYBOOL inc_matrow_space(MATrec *mat, int deltarows);
STATIC int mat_rowcompact(MATrec *mat);
STATIC int mat_colcompact(MATrec *mat, int prevcols);
STATIC MYBOOL inc_matcol_space(MATrec *mat, int deltacols);
STATIC MYBOOL inc_mat_space(MATrec *mat, int mindelta);
STATIC int mat_shiftrows(MATrec *mat, int *bbase, int delta);
STATIC int mat_shiftcols(MATrec *mat, int *bbase, int delta);
STATIC int mat_appendrow(MATrec *mat, int count, REAL *row, int *colno, REAL mult, MYBOOL checkrowmode);
STATIC int mat_appendcol(MATrec *mat, int count, REAL *column, int *rowno, REAL mult, MYBOOL checkrowmode);
MYBOOL mat_get_data(lprec *lp, int matindex, MYBOOL isrow, int **rownr, int **colnr, REAL **value);
MYBOOL mat_set_rowmap(MATrec *mat, int row_mat_index, int rownr, int colnr, int col_mat_index);
STATIC MYBOOL mat_indexrange(MATrec *mat, int index, MYBOOL isrow, int *startpos, int *endpos);
STATIC MYBOOL mat_validate(MATrec *mat);
STATIC int mat_findelm(MATrec *mat, int row, int column);
STATIC int mat_findins(MATrec *mat, int row, int column, int *insertpos, MYBOOL validate);
STATIC void mat_multcol(MATrec *mat, int col_nr, REAL mult);
STATIC REAL mat_getitem(MATrec *mat, int row, int column);
STATIC int mat_nonzeros(MATrec *mat);
STATIC int mat_collength(MATrec *mat, int colnr);
STATIC int mat_rowlength(MATrec *mat, int rownr);
STATIC void mat_multrow(MATrec *mat, int row_nr, REAL mult);
STATIC MYBOOL mat_setvalue(MATrec *mat, int Row, int Column, REAL Value, MYBOOL doscale);
STATIC MYBOOL mat_setrow(MATrec *mat, int rowno, int count, REAL *row, int *colno, MYBOOL doscale, MYBOOL checkrowmode);
STATIC MYBOOL mat_setcol(MATrec *mat, int colno, int count, REAL *column, int *rowno, MYBOOL doscale, MYBOOL checkrowmode);
STATIC int mat_checkcounts(MATrec *mat, int *rownum, int *colnum, MYBOOL freeonexit);
STATIC MYBOOL mat_transpose(MATrec *mat, MYBOOL col1_row0);

/* Refactorization and recomputation routine */
MYBOOL __WINAPI invert(lprec *lp, MYBOOL shiftbounds, MYBOOL final);

/* Sparse matrix operators */
STATIC MYBOOL fimprove(lprec *lp, REAL *pcol, int *nzidx, REAL roundzero);
STATIC void ftran(lprec *lp, REAL *rhsvector, int *nzidx, REAL roundzero);
STATIC MYBOOL bimprove(lprec *lp, REAL *rhsvector, int *nzidx, REAL roundzero);
STATIC void btran(lprec *lp, REAL *rhsvector, int *nzidx, REAL roundzero);
STATIC int prod_Ax(lprec *lp, int varset, REAL *input, int *nzinput, int range, REAL roundzero, REAL multfactor, REAL *output, int *nzoutput);
STATIC int prod_xA(lprec *lp, int varset, REAL *input, int *nzinput, int range, REAL roundzero, REAL ofscalar, REAL *output, int *nzoutput);
STATIC void prod_xA2(lprec *lp, REAL *prow, int prange, REAL proundzero, int *pnzprow,
                                REAL *drow, int drange, REAL droundzero, int *dnzdrow, REAL ofscalar);
STATIC MYBOOL fsolve(lprec *lp, int varin, REAL *pcol, int *nzidx, REAL roundzero, REAL ofscalar, MYBOOL prepareupdate);
STATIC MYBOOL bsolve(lprec *lp, int row_nr, REAL *rhsvector, int *nzidx, REAL roundzero, REAL ofscalar);
STATIC void bsolve_xA2(lprec *lp, int row_nr1, REAL *vector1, REAL roundzero1, int *nzvector1,
                                  int row_nr2, REAL *vector2, REAL roundzero2, int *nzvector2);

#ifdef __cplusplus
}
#endif

#endif /* HEADER_lp_matrix */
