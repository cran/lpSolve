
#include <string.h>
#include "commonlib.h"
#include "lp_lib.h"
#include "lp_report.h"
#include "lp_price.h"
#include "lp_pricePSE.h"
#include "lp_matrix.h"

#ifdef FORTIFY
# include "lp_fortify.h"
#endif


/* -------------------------------------------------------------------------
   Basic matrix routines in lp_solve v5.0+
   -------------------------------------------------------------------------
    Author:        Michel Berkelaar (to lp_solve v3.2),
                   Kjell Eikland
    Contact:       kjell.eikland@broadpark.no
    License terms: LGPL.

    Requires:      lp_lib.h, lp_pricerPSE.h, lp_matrix.h

    Release notes:
    v5.0.0  1 January 2004      First integrated and repackaged version.
    v5.0.1  7 May 2004          Added matrix transpose function.
    v5.1.0  20 July 2004        Reworked with flexible matrix storage model.

   ------------------------------------------------------------------------- */

STATIC MATrec *mat_create(lprec *lp, int rows, int columns, REAL epsvalue)
{
  MATrec *newmat;

  newmat = (MATrec *) calloc(1, sizeof(*newmat));
  newmat->lp = lp;

  newmat->rows_alloc = 0;
  newmat->columns_alloc = 0;
  newmat->mat_alloc = 0;

  inc_matrow_space(newmat, rows);
  newmat->rows = rows;
  inc_matcol_space(newmat, columns);
  newmat->columns = columns;
  inc_mat_space(newmat, 0);

  newmat->epsvalue = epsvalue;

  return( newmat );
}

STATIC void mat_free(MATrec **matrix)
{
#if MatrixColAccess==CAM_Record
  FREE((*matrix)->col_mat);
#else /*if MatrixColAccess==CAM_Vector*/
  FREE((*matrix)->col_mat_colnr);
  FREE((*matrix)->col_mat_rownr);
  FREE((*matrix)->col_mat_value);
#endif
  FREE((*matrix)->col_end);

#if MatrixRowAccess==RAM_Index
  FREE((*matrix)->row_mat);
#elif MatrixColAccess==CAM_Record
  FREE((*matrix)->row_mat);
#else /*if MatrixRowAccess==COL_Vector*/
  FREE((*matrix)->row_mat_colnr);
  FREE((*matrix)->row_mat_rownr);
  FREE((*matrix)->row_mat_value);
#endif
  FREE((*matrix)->row_end);

  FREE(*matrix);
}

STATIC MYBOOL inc_mat_space(MATrec *mat, int mindelta)
{
  int spaceneeded;

  if(mindelta <= 0)
    mindelta = MAX(mat->rows, mat->columns) + 1;
  if(mat->mat_alloc == 0)
    spaceneeded = mindelta;
  else
    spaceneeded = mat_nonzeros(mat) + mindelta;

  if(spaceneeded >= mat->mat_alloc) {
    /* Let's allocate at least MAT_START_SIZE entries */
    if(mat->mat_alloc < MAT_START_SIZE)
      mat->mat_alloc = MAT_START_SIZE;

    /* Increase the size by RESIZEFACTOR each time it becomes too small */
    while(spaceneeded >= mat->mat_alloc)
      mat->mat_alloc += mat->mat_alloc / RESIZEFACTOR;

#if MatrixColAccess==CAM_Record
    mat->col_mat = (MATitem *) realloc(mat->col_mat, (mat->mat_alloc) * sizeof(*(mat->col_mat)));
#else /*if MatrixColAccess==CAM_Vector*/
    allocINT(mat->lp,  &(mat->col_mat_colnr), mat->mat_alloc, AUTOMATIC);
    allocINT(mat->lp,  &(mat->col_mat_rownr), mat->mat_alloc, AUTOMATIC);
    allocREAL(mat->lp, &(mat->col_mat_value), mat->mat_alloc, AUTOMATIC);
#endif

#if MatrixRowAccess==RAM_Index
    allocINT(mat->lp, &(mat->row_mat), mat->mat_alloc, AUTOMATIC);
#elif MatrixColAccess==CAM_Record
    mat->row_mat = (MATitem *) realloc(mat->row_mat, (mat->mat_alloc) * sizeof(*(mat->row_mat)));
#else /*if MatrixColAccess==CAM_Vector*/
    allocINT(mat->lp,  &(mat->row_mat_colnr), mat->mat_alloc, AUTOMATIC);
    allocINT(mat->lp,  &(mat->row_mat_rownr), mat->mat_alloc, AUTOMATIC);
    allocREAL(mat->lp, &(mat->row_mat_value), mat->mat_alloc, AUTOMATIC);
#endif
  }
  return(TRUE);
}

STATIC MYBOOL inc_matrow_space(MATrec *mat, int deltarows)
{
  int    rowsum, oldrowsalloc /*, rowdelta = DELTAROWALLOC */;
  MYBOOL status = TRUE;

  /* Adjust lp row structures */
  if(mat->rows+deltarows > mat->rows_alloc) {
    /* deltarows = MAX(deltarows, rowdelta); */
    deltarows = MAX(DELTAROWALLOC, mat->rows+deltarows);
    oldrowsalloc = mat->rows_alloc;
    mat->rows_alloc += deltarows;
    rowsum = mat->rows_alloc + 1;

    status = allocINT(mat->lp, &mat->row_end, rowsum, AUTOMATIC);
    mat->row_end_valid = FALSE;
  }
  return( status );
}

STATIC MYBOOL inc_matcol_space(MATrec *mat, int deltacols)
{
  int    i,colsum, oldcolsalloc, newcolcount;
  MYBOOL status = TRUE;

  newcolcount = mat->columns+deltacols;

  if(newcolcount >= mat->columns_alloc) {

    /* Update memory allocation and sizes */
    oldcolsalloc = mat->columns_alloc;
    deltacols = MAX(DELTACOLALLOC, newcolcount);
    mat->columns_alloc += deltacols;
    colsum = mat->columns_alloc + 1;
    status = allocINT(mat->lp, &mat->col_end, colsum, AUTOMATIC);

    /* Update column pointers */
    if(oldcolsalloc == 0)
      mat->col_end[0] = 0;
    for(i = MIN(oldcolsalloc, mat->columns) + 1; i < colsum; i++)
      mat->col_end[i] = mat->col_end[i-1];
    mat->row_end_valid = FALSE;
  }
  return( status );
}

STATIC int mat_collength(MATrec *mat, int colnr)
{
  return( mat->col_end[colnr] - mat->col_end[colnr-1] );
}

STATIC int mat_rowlength(MATrec *mat, int rownr)
{
  if(mat_validate(mat)) {
    if(rownr <= 0)
      return( mat->row_end[0] );
    else
      return( mat->row_end[rownr] - mat->row_end[rownr-1] );
  }
  else
    return( 0 );
}

STATIC int mat_nonzeros(MATrec *mat)
{
  return( mat->col_end[mat->columns] );
}

STATIC MYBOOL mat_indexrange(MATrec *mat, int index, MYBOOL isrow, int *startpos, int *endpos)
{
#ifdef Paranoia
  if(isrow && ((index < 0) || (index > mat->rows)))
    return( FALSE );
  else if(!isrow && ((index < 1) || (index > mat->columns)))
    return( FALSE );
#endif

  if(isrow && mat_validate(mat)) {
    if(index == 0)
      *startpos = 0;
    else
      *startpos = mat->row_end[index-1];
    *endpos = mat->row_end[index];
  }
  else {
    *startpos = mat->col_end[index-1];
    *endpos = mat->col_end[index];
  }
  return( TRUE );
}

STATIC int mat_shiftrows(MATrec *mat, int *bbase, int delta)
{
  int     j, k, i, ii, thisrow, *colend, base;
  MYBOOL  preparecompact = FALSE;
  int     *rownr;

  if(delta == 0)
    return( 0 );
  base = abs(*bbase);

  if(delta > 0) {

    /* Insert row by simply incrementing existing row indeces */
    if(base <= mat->rows) {
      k = mat_nonzeros(mat);
      rownr = &COL_MAT_ROWNR(0);
      for(ii = 0; ii < k; ii++, rownr += matRowColStep) {
        if(*rownr >= base)
          *rownr += delta;
      }
    }

    /* Set defaults (actual basis set in separate procedure) */
    for(i = 0; i < delta; i++) {
      ii = base + i;
      mat->row_end[ii] = 0;
    }
  }
  else if(base <= mat->rows) {

    /* Check if we should prepare for compacting later
       (this is in order to speed up multiple row deletions) */
    preparecompact = (MYBOOL) (*bbase < 0);
    if(preparecompact)
      *bbase = my_flipsign((*bbase));

    /* First make sure we don't cross the row count border */
    if(base-delta-1 > mat->rows)
      delta = base - mat->rows - 1;

    /* Then scan over all entries shifting and updating rows indeces */
    if(preparecompact) {
      k = 0;
      for(j = 1, colend = mat->col_end + 1;
          j <= mat->columns; j++, colend++) {
        i = k;
        k = *colend;
        rownr = &COL_MAT_ROWNR(i);
        for(; i < k; i++, rownr += matRowColStep) {
          thisrow = *rownr;
          if(thisrow < base)
            continue;
          else if(thisrow >= base-delta)
            *rownr += delta;
          else
            *rownr = -1;
        }
      }
    }
    else {
      k = 0;
      ii = 0;
      for(j = 1, colend = mat->col_end + 1;
          j <= mat->columns; j++, colend++) {
        i = k;
        k = *colend;
        rownr = &COL_MAT_ROWNR(i);
        for(; i < k; i++, rownr += matRowColStep) {
          thisrow = *rownr;
          if(thisrow >= base) {
            if(thisrow >= base-delta)
              *rownr += delta;
            else
              continue;
          }
          if(ii != i) {
            COL_MAT_COPY(ii, i);
          }
          ii++;
        }
        *colend = ii;
      }
    }
  }
  return( 0 );
}

STATIC int mat_rowcompact(MATrec *mat)
{
  int     i, ii, j, k, nn, *colend;
  int     *rownr;

  nn = 0;
  k  = 0;
  ii = 0;
  for(j = 1, colend = mat->col_end + 1;
      j <= mat->columns; j++, colend++) {
    i = k;
    k = *colend;
    rownr = &COL_MAT_ROWNR(i);
    for(; i < k;
        i++, rownr += matRowColStep) {
      if(*rownr < 0) {
        nn++;
        continue;
      }
      if(ii != i) {
        COL_MAT_COPY(ii, i);
      }
      ii++;
    }
    *colend = ii;
  }
  return(nn);
}

STATIC int mat_colcompact(MATrec *mat, int prevcols)
{
  int  i, ii, j, k, n, nn, *colend, *newcolend, *colnr, newcolnr;

  nn = 0;
  k  = 0;
  ii = 0;
  newcolnr = 1;
  for(j = 1, colend = newcolend = mat->col_end + 1;
      j <= prevcols; j++, colend++) {
    n = 0;
    i = k;
    k = *colend;
    for(colnr = &COL_MAT_COLNR(i); i < k;
        i++, colnr += matRowColStep) {
      if(*colnr < 0) {
        n++;
        nn++;
        continue;
      }
      if(ii != i) {
        COL_MAT_COPY(ii, i);
        COL_MAT_COLNR(ii) = newcolnr;
      }
      ii++;
    }
    *newcolend = ii;
    if(n == 0) {
      newcolend++;
      newcolnr++;
    }
  }
  return(nn);
}

STATIC int mat_shiftcols(MATrec *mat, int *bbase, int delta)
{
  int     i, ii, k, n, base;


  k = 0;
  if(delta == 0)
    return( k );
  base = abs(*bbase);

  if(delta > 0) {
    /* Shift pointers right */
    for(ii = mat->columns; ii > base; ii--) {
      i = ii + delta;
      mat->col_end[i] = mat->col_end[ii];
    }
    /* Set defaults */
    for(i = 0; i < delta; i++) {
      ii = base + i;
      mat->col_end[ii] = mat->col_end[ii-1];
    }
  }
  else {

    /* Check if we should prepare for compacting later
       (this is in order to speed up multiple column deletions) */
    MYBOOL preparecompact = (MYBOOL) (*bbase < 0);
    if(preparecompact)
      *bbase = my_flipsign((*bbase));

    /* First make sure we don't cross the column count border */
    if(base-delta-1 > mat->columns)
      delta = base - mat->columns - 1;

    /* Then scan over all entries shifting and updating column indeces */
    if(preparecompact) {
      int *colnr;

      i = mat->col_end[base-1];
      k = mat->col_end[base-delta-1];
      for(colnr = &COL_MAT_COLNR(i); i < k;
          i++, colnr += matRowColStep)
        *colnr = -1;
    }
    else {
      /* Delete sparse matrix data, if required */
      if(base <= mat->columns) {

        i = mat->col_end[base-1];          /* Beginning of data to be deleted */
        ii = mat->col_end[base-delta-1];   /* Beginning of data to be shifted left */
        k = ii-i;                          /* Number of entries to be deleted */
        if(k > 0) {
          n = mat_nonzeros(mat);
          n -= ii;
          COL_MAT_MOVE(i, ii, n);
        }

        /* Update indexes */
        for(i = base; i <= mat->columns + delta; i++) {
          ii = i - delta;
          mat->col_end[i] = mat->col_end[ii] - k;
        }
      }
    }
  }
  return( k );
}

STATIC MYBOOL mat_setcol(MATrec *mat, int colno, int count, REAL *column, int *rowno, MYBOOL doscale, MYBOOL checkrowmode)
{
  int    i, jj = 0, elmnr, orignr, newnr, firstrow;
  MYBOOL *addto = NULL, isA, isNZ;
  REAL   value;

  /* Check if we are in row order mode and should add as row instead;
     the matrix will be transposed at a later stage */
  if(checkrowmode && mat->is_roworder)
    return( mat_setrow(mat, colno, count, column, rowno, doscale, FALSE) );

  /* Initialize and validate */
  isA = (MYBOOL) (mat == mat->lp->matA);
  isNZ = (MYBOOL) (rowno != NULL);
  if(!isNZ)
    count = mat->lp->rows;
  else if((count < 0) || (count > mat->rows))
    return( FALSE );
  if(isNZ && (count > 0)) {
    if(count > 1)
      sortREALByINT(column, rowno, count, 0, TRUE);
    if((rowno[0] < 0) || (rowno[count-1] > mat->rows))
      return( FALSE );
  }

  /* Optionally tally and map the new non-zero values */
  firstrow = mat->rows + 1;
  if(isNZ) {
    newnr = count;
    if(newnr) {
      firstrow = rowno[0];
      jj = rowno[newnr - 1];
    }
  }
  else {
    newnr = 0;
    if(!allocMYBOOL(mat->lp, &addto, mat->rows + 1, TRUE)) {
      mat->lp->spx_status = NOMEMORY;
      return( FALSE );
    }
    for(i = mat->rows; i >= 0; i--) {
      if(fabs(column[i]) > mat->epsvalue) {
        addto[i] = TRUE;
        firstrow = i;
        newnr++;
      }
    }
  }

  /* Shift existing column data and adjust position indeces */
  orignr = mat_collength(mat, colno);
  elmnr = newnr - orignr;
  i = mat_nonzeros(mat) - mat->col_end[colno];
  if((elmnr != 0) && (i > 0)) {
    COL_MAT_MOVE(mat->col_end[colno] + elmnr, mat->col_end[colno], i);
  }
  if(elmnr != 0)
    for(i = colno; i <= mat->columns; i++)
      mat->col_end[i] += elmnr;

  /* We are now ready to copy the new data */
  if(isNZ) {
    for(jj = mat->col_end[colno-1], i = 0; i < count; jj++, i++) {
      value = column[i];
#ifdef DoMatrixRounding
      value = roundToPrecision(value, mat->epsvalue);
#endif
      if(isA && doscale)
        value = scaled_mat(mat->lp, value, rowno[i], colno);
      if(isA)
        value = my_chsign(is_chsign(mat->lp, rowno[i]), value);
      SET_MAT_ijA(jj, rowno[i], colno, value);
    }
  }
  else {
    jj = mat->col_end[colno-1];
    for(i = firstrow; i <= mat->rows; i++) {
      if(!addto[i])
        continue;
      value = column[i];
#ifdef DoMatrixRounding
      value = roundToPrecision(value, mat->epsvalue);
#endif
      if(isA && doscale)
        value = scaled_mat(mat->lp, value, i, colno);
      if(isA)
        value = my_chsign(is_chsign(mat->lp, i), value);
      SET_MAT_ijA(jj, i, colno, value);
      jj++;
    }
  }

  /* Finish and return */
  FREE(addto);
  mat->row_end_valid = FALSE;
  return( TRUE );

}

STATIC MYBOOL mat_setrow(MATrec *mat, int rowno, int count, REAL *row, int *colno, MYBOOL doscale, MYBOOL checkrowmode)
{
  int    k, kk, i, ii, j, jj = 0, jj_j, elmnr, orignr, newnr, firstcol, rownr, colnr;
  MYBOOL *addto = NULL, isA, isNZ;
  REAL   value;

  /* Check if we are in row order mode and should add as column instead;
     the matrix will be transposed at a later stage */
  if(checkrowmode && mat->is_roworder)
    return( mat_setcol(mat, rowno, count, row, colno, doscale, FALSE) );

  /* Do initialization and validation */
  if(!mat_validate(mat))
    return( FALSE );
  isA = (MYBOOL) (mat == mat->lp->matA);
  isNZ = (MYBOOL) (colno != NULL);
  if(!isNZ)
    count = mat->columns;
  else if((count < 0) || (count > mat->columns))
    return( FALSE );
  if(isNZ && (count > 0)) {
    if(count > 1)
      sortREALByINT(row, colno, count, 0, TRUE);
    if((colno[0] < 1) || (colno[count-1] > mat->columns))
      return( FALSE );
  }

  /* Optionally tally and map the new non-zero values */
  firstcol = mat->columns + 1;
  if(isNZ) {
    newnr = count;
    if(newnr)
      firstcol = colno[0];
  }
  else {
    newnr = 0;
    if(!allocMYBOOL(mat->lp, &addto, mat->columns + 1, TRUE)) {
      mat->lp->spx_status = NOMEMORY;
      return( FALSE );
    }
    for(i = mat->columns; i >= 1; i--) {
      if(fabs(row[i]) > mat->epsvalue) {
        addto[i] = TRUE;
        firstcol = i;
        newnr++;
      }
    }
  }

  /* Pack initial entries if existing row data has a lower column
     start index than the first index of the new vector */
  orignr = mat_nonzeros(mat);
  k = newnr - mat_rowlength(mat, rowno);
  kk = 0;
  if(rowno == 0)
    ii = 0;
  else
    ii = mat->row_end[rowno-1];
  if((orignr == 0) || (ii >= orignr))
    j = 1;
  else
    j = ROW_MAT_COLNR(ii);
  jj = mat->col_end[firstcol - 1];
  if(jj >= orignr)
    colnr = firstcol;
  else
    colnr = COL_MAT_COLNR(jj);
  if(j < colnr) {
    elmnr = mat->col_end[j-1];
    jj = elmnr;
    for( ; j < colnr; j++) {
      /* Shift entries in current column */
      for( ; jj < mat->col_end[j]; jj++) {
        if(COL_MAT_ROWNR(jj) != rowno) {
          COL_MAT_COPY(elmnr, jj);
          elmnr++;
        }
      }
      /* Update next column start index */
      mat->col_end[j] = elmnr;
    }
    jj_j = jj - elmnr;  /* The shrinkage count */
  }
  else
    jj_j = 0;

  /* Make sure we have sufficient space for any additional entries and move existing data down;
     this ensures that we only have to relocate matrix elements up in the next stage */
  jj_j = MAX(0, newnr - jj_j);
  if(jj_j > 0) {
    if(!inc_mat_space(mat, jj_j)) {
      FREE(addto);
      return( FALSE );
    }
    COL_MAT_MOVE(jj+jj_j, jj, orignr-jj);
    jj += jj_j;
  }

  /* Handle case where the matrix was empty before */
  if(orignr == 0) {
    if(isNZ)
      elmnr = count;
    else
      elmnr = mat->columns;
    jj_j = 0;
    for(newnr = 0; newnr < elmnr; newnr++) {
      if(isNZ)
        colnr = colno[newnr];
      else
        colnr = newnr + 1;
      if(isNZ || addto[colnr]) {
        if(isNZ)
          value = row[newnr];
        else
          value = row[colnr];
#ifdef DoMatrixRounding
        value = roundToPrecision(value, mat->epsvalue);
#endif
        if(isA && doscale)
          value = scaled_mat(mat->lp, value, rowno, colnr);
        if(isA)
          value = my_chsign(is_chsign(mat->lp, rowno), value);
        SET_MAT_ijA(jj_j, rowno, colnr, value);
        jj_j++;
      }
      /* Update column start position if we have crossed a column */
      while(colnr >= firstcol) {
        mat->col_end[firstcol] = jj_j;
        firstcol++;
      }
    }

    /* Make sure we update tail empty column offsets */
    while(firstcol <= mat->columns) {
      mat->col_end[firstcol] = jj_j;
      firstcol++;
    }
    jj_j = 0;
  }

  /* Start from the top of the first non-zero column of the new row */
  elmnr = orignr + jj_j;
  if(jj < elmnr) {
    newnr = 0;
    j = jj - mat->col_end[firstcol - 1];
    colnr = firstcol;
    while((jj < elmnr) || (newnr < count)) {

      /* Update column start position if we have crossed a column */
      while(colnr > firstcol) {
        mat->col_end[firstcol] = kk;
        firstcol++;
      }

      /* See if we have a row equal to or greater than the target row */
      jj_j = jj - j;
      if(jj < elmnr) {
        rownr = COL_MAT_ROWNR(jj);
        colnr = COL_MAT_COLNR(jj);
      }
      else {
        rownr = rowno;
        colnr = mat->columns + 1;
      }

      if(isNZ) {
        if(newnr < count)
          kk = colno[newnr];
        else
          kk = mat->columns + 1;
      }
      else
        kk = newnr + 1;

      /* Test if there is an available new item ... */
      if((isNZ && (kk > colnr)) ||                    /* If this is not the case */
         (!isNZ && !addto[kk])) {
        /* DELETE if there is an existing value */
        if(!isNZ)
          newnr++;
        if(rownr == rowno) {
          kk = jj_j;
          j++;
          jj++;
          continue;
        }
        /* KEEP otherwise and move entry up */
      }
      else if((colnr > kk) ||                         /* Existing column index > new => INSERT */
              ((colnr == kk) && (rownr >= rowno)) ) { /* Same column index, existing row >= target row => INSERT/REPLACE */

        if(isNZ)
          value = row[newnr];
        else
          value = row[newnr+1];
        newnr++;
#ifdef DoMatrixRounding
        value = roundToPrecision(value, mat->epsvalue);
#endif
        if(isA && doscale)
          value = scaled_mat(mat->lp, value, rowno, colnr);
        if(isA)
          value = my_chsign(is_chsign(mat->lp, rowno), value);
        SET_MAT_ijA(jj_j, rowno, kk, value);

        /* Adjust if we have inserted an element */
        if((colnr > kk) || (rownr > rowno)) {
          j--;
          jj--;
        }
        colnr = kk;
        kk = jj_j;
        jj++;
        continue;
      }

      /* Shift the matrix element up by the active difference */
      if(jj_j != jj) {
        COL_MAT_COPY(jj_j, jj);
      }
      kk = jj_j;
      jj++;

    }

    /* Update pending / incomplete column start position */
    while(colnr > firstcol) {
      mat->col_end[firstcol] = kk;
      firstcol++;
    }

    /* Make sure we update tail column offsets */
    jj_j = jj - j;
    while(firstcol <= mat->columns) {
      mat->col_end[firstcol] = jj_j;
      firstcol++;
    }
  }

  mat->row_end_valid = FALSE;

  FREE(addto);
  return( TRUE );

}

STATIC int mat_appendrow(MATrec *mat, int count, REAL *row, int *colno, REAL mult, MYBOOL checkrowmode)
{
  int    i, j, jj = 0, stcol, elmnr, orignr, newnr, firstcol;
  MYBOOL *addto = NULL, isA, isNZ;
  REAL   value;

  /* Check if we are in row order mode and should add as column instead;
     the matrix will be transposed at a later stage */
  if(checkrowmode && mat->is_roworder)
    return( mat_appendcol(mat, count, row, colno, mult, FALSE) );

  /* Do initialization and validation */
  isA = (MYBOOL) (mat == mat->lp->matA);
  isNZ = (MYBOOL) (colno != NULL);
  if(isNZ && (count > 0)) {
    if(count > 1)
      sortREALByINT(row, colno, count, 0, TRUE);
    if((colno[0] < 1) || (colno[count-1] > mat->columns))
      return( 0 );
  }
  else if(row != NULL)
    row[0] = 0;

  /* Optionally tally and map the new non-zero values */
  firstcol = mat->columns + 1;
  if(isNZ) {
    newnr = count;
    if(newnr) {
      firstcol = colno[0];
      jj = colno[newnr - 1];
    }
  }
  else {
    newnr = 0;
    if(!allocMYBOOL(mat->lp, &addto, mat->columns + 1, TRUE)) {
      mat->lp->spx_status = NOMEMORY;
      return( newnr );
    }
    for(i = mat->columns; i >= 1; i--) {
      if(fabs(row[i]) > mat->epsvalue) {
        addto[i] = TRUE;
        firstcol = i;
        newnr++;
      }
    }
  }

  /* Make sure we have sufficient space */
  if(!inc_mat_space(mat, newnr)) {
    FREE(addto);
    return( 0 );
  }

  /* Insert the non-zero constraint values */
  orignr = mat_nonzeros(mat) - 1;
  elmnr = orignr + newnr;

  for(j = mat->columns; j >= firstcol; j--) {
    stcol = mat->col_end[j] - 1;
    mat->col_end[j] = elmnr + 1;

   /* Add a new non-zero entry */
    if(((isNZ) && (j == jj)) || ((addto != NULL) && (addto[j]))) {
      newnr--;
      if(isNZ) {
        value = row[newnr];
        if(newnr)
          jj = colno[newnr - 1];
        else
          jj = 0;
      }
      else
        value = row[j];
#ifdef DoMatrixRounding
      value = roundToPrecision(value, mat->epsvalue);
#endif
      value *= mult;
      if(isA)
        value = scaled_mat(mat->lp, value, mat->rows, j);
      SET_MAT_ijA(elmnr, mat->rows, j, value);
      elmnr--;
    }

   /* Shift previous column entries down */
    i = stcol - mat->col_end[j-1] + 1;
    if(i > 0) {
      orignr -= i;
      elmnr  -= i;
      COL_MAT_MOVE(elmnr+1, orignr+1, i);
    }
  }

  FREE(addto);

  return( newnr );

}

STATIC int mat_appendcol(MATrec *mat, int count, REAL *column, int *rowno, REAL mult, MYBOOL checkrowmode)
{
  int     i, row, elmnr, lastnr;
  REAL    value;
  MYBOOL  isA, isNZ;

  /* Check if we are in row order mode and should add as column instead;
     the matrix will be transposed at a later stage */
  if(checkrowmode && mat->is_roworder)
    return( mat_appendrow(mat, count, column, rowno, mult, FALSE) );

  /* Make sure we have enough space */
  if(!inc_mat_space(mat, mat->rows+1))
    return( 0 );

  /* Do initialization and validation */
  isA = (MYBOOL) (mat == mat->lp->matA);
  isNZ = (MYBOOL) (rowno != NULL);
  if(isNZ && (count > 0)) {
    if(count > 1)
      sortREALByINT(column, rowno, count, 0, TRUE);
    if((rowno[0] < 0))
      return( 0 );
  }
  if(rowno != NULL)
    count--;

  /* Append sparse regular constraint values */
  elmnr = mat->col_end[mat->columns - 1];
  if(column != NULL) {
    row = -1;
    for(i = ((isNZ || !mat->is_roworder) ? 0 : 1); i <= count ; i++) {
      value = column[i];
      if(fabs(value) > mat->epsvalue) {
        if(isNZ) {
          lastnr = row;
          row = rowno[i];
          /* Check if we have come to the Lagrangean constraints */
          if(row > mat->rows)
            break;
          if(row <= lastnr)
            return( -1 );
        }
        else
          row = i;
#ifdef DoMatrixRounding
        value = roundToPrecision(value, mat->epsvalue);
#endif
        if(mat->is_roworder)
          value *= mult;
        else if(isA) {
          value = my_chsign(is_chsign(mat->lp, row), value);
          value = scaled_mat(mat->lp, value, row, mat->columns);
        }

       /* Store the item and update counters */
        SET_MAT_ijA(elmnr, row, mat->columns, value);
        elmnr++;
      }
    }

   /* Fill dense Lagrangean constraints */
    if(get_Lrows(mat->lp) > 0)
      mat_appendcol(mat->lp->matL, get_Lrows(mat->lp), column+mat->rows, NULL, mult, checkrowmode);

  }

 /* Set end of data */
  mat->col_end[mat->columns] = elmnr;

  return( mat->col_end[mat->columns] - mat->col_end[mat->columns-1] );
}

STATIC int mat_checkcounts(MATrec *mat, int *rownum, int *colnum, MYBOOL freeonexit)
{
  int i, j, n;
  int *rownr;

  if(rownum == NULL)
    allocINT(mat->lp, &rownum, mat->rows + 1, TRUE);
  if(colnum == NULL)
    allocINT(mat->lp, &colnum, mat->columns + 1, TRUE);

  for(i = 1 ; i <= mat->columns; i++) {
    j = mat->col_end[i - 1];
    n = mat->col_end[i];
    rownr = &COL_MAT_ROWNR(j);
    for(; j < n;
        j++, rownr += matRowColStep) {
      colnum[i]++;
      rownum[*rownr]++;
    }
  }

  n = 0;
  if((mat->lp->do_presolve != PRESOLVE_NONE) &&
     (mat->lp->spx_trace || (mat->lp->verbose > NORMAL))) {
    for(j = 1; j <= mat->columns; j++)
      if(colnum[j] == 0) {
        n++;
        report(mat->lp, FULL, "mat_checkcounts: Variable %s is not used in any constraints\n",
                              get_col_name(mat->lp, j));
      }
    for(i = 0; i <= mat->rows; i++)
      if(rownum[i] == 0) {
        n++;
        report(mat->lp, FULL, "mat_checkcounts: Constraint %s empty\n",
                              get_row_name(mat->lp, i));
      }
  }

  if(freeonexit) {
    FREE(rownum);
    FREE(colnum);
  }

  return( n );

}

STATIC MYBOOL mat_validate(MATrec *mat)
/* Routine to make sure that row mapping arrays are valid */
{
  int     i, j, je, *rownum;
  int     *rownr, *colnr;

  if(!mat->row_end_valid) {

    MEMCLEAR(mat->row_end, mat->rows + 1);
    allocINT(mat->lp, &rownum, mat->rows + 1, TRUE);

    /* First tally row counts and then cumulate them */
    j = mat_nonzeros(mat);
    rownr = &COL_MAT_ROWNR(0);
    for(i = 0; i < j; i++, rownr += matRowColStep)
      mat->row_end[*rownr]++;
    for(i = 1; i <= mat->rows; i++)
      mat->row_end[i] += mat->row_end[i - 1];

    /* Calculate the column index for every non-zero */
    for(i = 1; i <= mat->columns; i++) {
      j = mat->col_end[i - 1];
      je = mat->col_end[i];
      rownr = &COL_MAT_ROWNR(j);
      colnr = &COL_MAT_COLNR(j);
      for(; j < je; j++, rownr += matRowColStep, colnr += matRowColStep) {
#ifdef Paranoiax
        if(*colnr < i) {
          report(mat->lp, SEVERE, "mat_validate: Internal storage inconsistency for row %d, column %d\n",
                                  *rownr, *colnr);
          mat->lp->spx_status = UNKNOWNERROR;
          return(FALSE);
        }
#endif
        *colnr = i;
        if(*rownr == 0)
          mat_set_rowmap(mat, rownum[*rownr],
                              *rownr, i, j);
        else
          mat_set_rowmap(mat, mat->row_end[*rownr - 1] + rownum[*rownr],
                              *rownr, i, j);
        rownum[*rownr]++;
      }
    }

    FREE(rownum);
    mat->row_end_valid = TRUE;
  }

  if(mat == mat->lp->matA)
    mat->lp->model_is_valid = TRUE;
  return( TRUE );
}

MYBOOL mat_get_data(lprec *lp, int matindex, MYBOOL isrow, int **rownr, int **colnr, REAL **value)
{
  MATrec *mat = lp->matA;

#if MatrixRowAccess == RAM_Index
  if(isrow)
    matindex = mat->row_mat[matindex];
  if(rownr != NULL)
    *rownr = &COL_MAT_ROWNR(matindex);
  if(colnr != NULL)
    *colnr = &COL_MAT_COLNR(matindex);
  if(value != NULL)
    *value = &COL_MAT_VALUE(matindex);

#else
  if(isrow) {
    if(rownr != NULL)
      *rownr = &ROW_MAT_ROWNR(matindex);
    if(colnr != NULL)
      *colnr = &ROW_MAT_COLNR(matindex);
    if(value != NULL)
      *value = &ROW_MAT_VALUE(matindex);
  }
  else {
    if(rownr != NULL)
      *rownr = &COL_MAT_ROWNR(matindex);
    if(colnr != NULL)
      *colnr = &COL_MAT_COLNR(matindex);
    if(value != NULL)
      *value = &COL_MAT_VALUE(matindex);
  }

#endif

  return( TRUE );
}


MYBOOL mat_set_rowmap(MATrec *mat, int row_mat_index, int rownr, int colnr, int col_mat_index)
{
#if MatrixRowAccess == RAM_Index
  mat->row_mat[row_mat_index] = col_mat_index;

#elif MatrixColAccess==CAM_Record
  mat->row_mat[row_mat_index].rownr = rownr;
  mat->row_mat[row_mat_index].colnr = colnr;
  mat->row_mat[row_mat_index].value = COL_MAT_VALUE(col_mat_index);

#else /* if MatrixColAccess==CAM_Vector */
  mat->row_mat_rownr[row_mat_index] = rownr;
  mat->row_mat_colnr[row_mat_index] = colnr;
  mat->row_mat_value[row_mat_index] = COL_MAT_VALUE(col_mat_index);

#endif

  return( TRUE );
}

/* Implement combined binary/linear sub-search for matrix look-up */
int mat_findelm(MATrec *mat, int row, int column)
{
  int low, high, mid, item;

#if 0
  if(mat->row_end_valid && (row > 0) &&
     (ROW_MAT_COLNR(mat->row_mat[(low = mat->row_end[row-1])]) == column))
    return(low);
#endif

  if((column < 1) || (column > mat->columns)) {
    report(mat->lp, IMPORTANT, "mat_findelm: Column %d out of range\n", column);
    return( -1 );
  }
  if((row < 0) || (row > mat->rows)) {
    report(mat->lp, IMPORTANT, "mat_findelm: Row %d out of range\n", row);
    return( -1 );
  }

  low = mat->col_end[column - 1];
  high = mat->col_end[column] - 1;
  if(low > high)
    return( -2 );

 /* Do binary search logic */
  mid = (low+high) / 2;
  item = COL_MAT_ROWNR(mid);
  while(high - low > LINEARSEARCH) {
    if(item < row) {
      low = mid + 1;
      mid = (low+high) / 2;
      item = COL_MAT_ROWNR(mid);
    }
    else if(item > row) {
      high = mid - 1;
      mid = (low+high) / 2;
      item = COL_MAT_ROWNR(mid);
    }
    else {
      low = mid;
      high = mid;
    }
  }

 /* Do linear scan search logic */
  if((high > low) && (high - low <= LINEARSEARCH)) {
    item = COL_MAT_ROWNR(low);
    while((low < high) && (item < row)) {
      low++;
      item = COL_MAT_ROWNR(low);
    }
    if(item == row)
      high = low;
  }

  if((low == high) && (row == item))
    return( low );
  else
    return( -2 );
}

int mat_findins(MATrec *mat, int row, int column, int *insertpos, MYBOOL validate)
{
  int low, high, mid, item, exitvalue, insvalue;

#if 0
  if(mat->row_end_valid && (row > 0) &&
     (ROW_MAT_COLNR(mat->row_mat[(low = mat->row_end[row-1])]) == column)) {
    insvalue = low;
    exitvalue = low;
    goto Done;
  }
#endif

  insvalue = -1;

  if((column < 1) || (column > mat->columns)) {
    if((column > 0) && !validate) {
      insvalue = mat->col_end[mat->columns];
      exitvalue = -2;
      goto Done;
    }
    report(mat->lp, IMPORTANT, "mat_findins: Column %d out of range\n", column);
    exitvalue = -1;
    goto Done;
  }
  if((row < 0) || (row > mat->rows)) {
    if((row >= 0) && !validate) {
      insvalue = mat->col_end[column];
      exitvalue = -2;
      goto Done;
    }
    report(mat->lp, IMPORTANT, "mat_findins: Row %d out of range\n", row);
    exitvalue = -1;
    goto Done;
  }

  low = mat->col_end[column - 1];
  insvalue = low;
  high = mat->col_end[column] - 1;
  if(low > high) {
    exitvalue = -2;
    goto Done;
  }

 /* Do binary search logic */
  mid = (low+high) / 2;
  item = COL_MAT_ROWNR(mid);
  while(high - low > LINEARSEARCH) {
    if(item < row) {
      low = mid + 1;
      mid = (low+high) / 2;
      item = COL_MAT_ROWNR(mid);
    }
    else if(item > row) {
      high = mid - 1;
      mid = (low+high) / 2;
      item = COL_MAT_ROWNR(mid);
    }
    else {
      low = mid;
      high = mid;
    }
  }

 /* Do linear scan search logic */
  if((high > low) && (high - low <= LINEARSEARCH)) {
    item = COL_MAT_ROWNR(low);
    while((low < high) && (item < row)) {
      low++;
      item = COL_MAT_ROWNR(low);
    }
    if(item == row)
      high = low;
  }

  insvalue = low;
  if((low == high) && (row == item))
    exitvalue = low;
  else {
    if((low < mat->col_end[column]) && (COL_MAT_ROWNR(low) < row))
      insvalue++;
    exitvalue = -2;
  }

Done:
  if(insertpos != NULL)
    (*insertpos) = insvalue;
  return( exitvalue );
}

STATIC REAL mat_getitem(MATrec *mat, int row, int column)
{
  int elmnr;

#ifdef DirectOverrideOF
  if((row == 0) && (mat == mat->lp->matA) && (mat->lp->OF_override != NULL))
    return( mat->lp->OF_override[column] );
  else
#endif
  {
    elmnr = mat_findelm(mat, row, column);
    if(elmnr >= 0)
      return( COL_MAT_VALUE(elmnr) );
    else
      return( 0 );
  }
}

STATIC void mat_multrow(MATrec *mat, int row_nr, REAL mult)
{
  int i, k1, k2;

#if 0
  if(row_nr == 0) {
    k2 = mat->col_end[0];
    for(i = 1; i <= mat->columns; i++) {
      k1 = k2;
      k2 = mat->col_end[i];
      if((k1 < k2) && (COL_MAT_ROWNR(k1) == row_nr))
        COL_MAT_VALUE(k1) *= mult;
    }
  }
  else if(mat_validate(mat)) {
    if(row_nr == 0)
      k1 = 0;
    else
#else
  if(mat_validate(mat)) {
    if(row_nr == 0)
      k1 = 0;
    else
#endif
    k1 = mat->row_end[row_nr-1];
    k2 = mat->row_end[row_nr];
    for(i = k1; i < k2; i++)
      ROW_MAT_VALUE(i) *= mult;
  }
}

STATIC void mat_multcol(MATrec *mat, int col_nr, REAL mult)
{
  int    i, ie;
  MYBOOL isA;

#ifdef Paranoia
  if(col_nr < 1 || col_nr > mat->columns) {
    report(mat->lp, IMPORTANT, "mult_column: Column %d out of range\n", col_nr);
    return;
  }
#endif
  if(mult == 1.0)
    return;

  isA = (MYBOOL) (mat == mat->lp->matA);

  ie = mat->col_end[col_nr];
  for(i = mat->col_end[col_nr - 1]; i < ie; i++)
    COL_MAT_VALUE(i) *= mult;
  if(isA && (get_Lrows(mat->lp) > 0))
    mat_multcol(mat->lp->matL, col_nr, mult);
}

STATIC MYBOOL mat_setvalue(MATrec *mat, int Row, int Column, REAL Value, MYBOOL doscale)
{
  int    elmnr, lastelm, i;
  MYBOOL isA;

  /* This function is inefficient if used to add new matrix entries in
     other places than at the end of the matrix. OK for replacing existing
     a non-zero value with another non-zero value */
  isA = (MYBOOL) (mat == mat->lp->matA);
  if(isA && mat->is_roworder)
    mat_transpose(mat, TRUE);

  /* Set small numbers to zero */
  if(fabs(Value) < mat->epsvalue)
    Value = 0;
#ifdef DoMatrixRounding
  else
    Value = roundToPrecision(Value, mat->epsvalue);
#endif

  /* Check if we need to update column space */
  if(Column > mat->columns) {
    if(isA)
      inc_col_space(mat->lp, Column - mat->columns);
    else
      inc_matcol_space(mat, Column - mat->columns);
  }

  /* Find out if we already have such an entry, or return insertion point */
  i = mat_findins(mat, Row, Column, &elmnr, FALSE);
  if(i == -1)
    return(FALSE);

  if(isA) {
    if((Row > 0) && mat->lp->is_basic[mat->rows+Column])
      mat->lp->basis_valid = FALSE;
    mat->lp->doInvert = TRUE;
  }

  if(i >= 0) {
    /* there is an existing entry */
    if(fabs(Value) > mat->epsvalue) { /* we replace it by something non-zero */
      if(isA) {
        Value = my_chsign(is_chsign(mat->lp, Row), Value);
        if(doscale && mat->lp->scaling_used)
          Value = scaled_mat(mat->lp, Value, Row,Column);
      }
      COL_MAT_VALUE(elmnr) = Value;
    }
    else { /* setting existing non-zero entry to zero. Remove the entry */
      /* This might remove an entire column, or leave just a bound. No
          nice solution for that yet */

      /* Shift up tail end of the matrix */
      lastelm = mat_nonzeros(mat);
#if 0
      for(i = elmnr; i < lastelm ; i++) {
        COL_MAT_COPY(i, i + 1);
      }
#else
      lastelm -= elmnr;
      COL_MAT_MOVE(elmnr, elmnr + 1, lastelm);
#endif
      for(i = Column; i <= mat->columns; i++)
        mat->col_end[i]--;

      mat->row_end_valid = FALSE;
    }
  }
  else if(fabs(Value) > mat->epsvalue) {
    /* no existing entry. make new one only if not nearly zero */
    /* check if more space is needed for matrix */
    if (!inc_mat_space(mat, 1))
      return(FALSE);

    if(Column > mat->columns) {
      i = mat->columns+1;
      if(isA)
        shift_coldata(mat->lp, i, Column - mat->columns);
      else
        mat_shiftcols(mat, &i, Column - mat->columns);
    }

    /* Shift down tail end of the matrix by one */
    lastelm = mat_nonzeros(mat);
#if 1 /* Does compiler optimization work better here? */
    for(i = lastelm; i > elmnr ; i--) {
      COL_MAT_COPY(i, i - 1);
    }
#else
    lastelm -= elmnr - 1;
    COL_MAT_MOVE(elmnr + 1, elmnr, lastelm);
#endif

    /* Set new element */
    if(isA) {
      Value = my_chsign(is_chsign(mat->lp, Row), Value);
      if(doscale)
        Value = scaled_mat(mat->lp, Value, Row, Column);
    }
    SET_MAT_ijA(elmnr, Row, Column, Value);

    /* Update column indexes */
    for(i = Column; i <= mat->columns; i++)
      mat->col_end[i]++;

    mat->row_end_valid = FALSE;
  }

  if(isA && (mat->lp->var_is_free != NULL) && (mat->lp->var_is_free[Column] > 0))
    return( mat_setvalue(mat, Row, mat->lp->var_is_free[Column], -Value, doscale) );
  return(TRUE);
}

STATIC int mat_findcolumn(MATrec *mat, int matindex)
{
  int j;

  for(j = 1; j <= mat->columns; j++) {
    if(matindex < mat->col_end[j])
      break;
  }
  return(j);
}

STATIC MYBOOL mat_transpose(MATrec *mat, MYBOOL col1_row0)
{
  int     i, j, nz, k;
  MYBOOL  status;

  if(col1_row0 && !mat->is_roworder)
    return( FALSE );

  status = mat_validate(mat);
  if(status) {

    /* Create a column-ordered sparse element list; "column" index must be shifted */
    nz = mat_nonzeros(mat);
    if(nz > 0) {
#if MatrixColAccess==CAM_Record
      MATitem *newmat;
      newmat = (MATitem *) malloc((mat->mat_alloc) * sizeof(*(mat->col_mat)));
      if(col1_row0) {  /* Transposition of fast constraint add mode */
        for(i = nz-1; i >= 0 ; i--) {
          newmat[i] = mat->col_mat[mat->row_mat[i]];
          newmat[i].row_nr = newmat[i].col_nr-1;
        }
      }
      else {           /* Transposition where row index 0 becomes columns+1 */
        j = mat->row_end[0];
        for(i = nz-1; i >= j ; i--) {
          k = i-j;
          newmat[k] = mat->col_mat[mat->row_mat[i]];
          newmat[k].row_nr = newmat[k].col_nr;
        }
        for(i = j-1; i >= 0 ; i--) {
          k = nz-j+i;
          newmat[k] = mat->col_mat[mat->row_mat[i]];
          newmat[k].row_nr = newmat[k].col_nr;
        }
      }
      swapPTR((void **) &mat->col_mat, (void **) &newmat);
      FREE(newmat);
#else /*if MatrixColAccess==CAM_Vector*/
      REAL *newValue = NULL;
      int  *newRownr = NULL;
      allocREAL(mat->lp, &newValue, mat->mat_alloc, FALSE);
      allocINT(mat->lp, &newRownr, mat->mat_alloc, FALSE);

      if(col1_row0) {  /* Transposition of fast constraint add mode */
        for(i = nz-1; i >= 0 ; i--) {
          newValue[i] = ROW_MAT_VALUE(i);
          newRownr[i] = ROW_MAT_COLNR(i) - 1;
        }
      }
      else {           /* Transposition where row index 0 becomes columns+1 */
        j = mat->row_end[0];
        for(i = nz-1; i >= j ; i--) {
          k = i-j;
          newValue[k] = ROW_MAT_VALUE(i);
          newRownr[k] = ROW_MAT_COLNR(i);
        }
        for(i = j-1; i >= 0 ; i--) {
          k = nz-j+i;
          newValue[k] = ROW_MAT_VALUE(i);
          newRownr[k] = ROW_MAT_COLNR(i);
        }
      }

      swapPTR((void **) &mat->col_mat_rownr, (void **) &newRownr);
      swapPTR((void **) &mat->col_mat_value, (void **) &newValue);
      FREE(newValue);
      FREE(newRownr);
#endif
    }

    /* Transfer row start to column start position; must adjust for different offsets */
    if(mat->rows == mat->rows_alloc)
      inc_matcol_space(mat, 1);
    if(col1_row0)
      mat->columns--;
    else {
      j = mat->row_end[0];
      for(i = mat->rows; i >= 1; i--)
        mat->row_end[i] -= j;
      mat->rows++;
      mat->row_end[mat->rows] = nz;
    }
    swapPTR((void **) &mat->row_end, (void **) &mat->col_end);

    /* Swap array sizes */
    swapINT(&mat->rows, &mat->columns);
    swapINT(&mat->rows_alloc, &mat->columns_alloc);

    /* Finally set current storage mode */
    mat->is_roworder = (MYBOOL) !mat->is_roworder;
    mat->row_end_valid = FALSE;
  }
  return(status);
}

/* ---------------------------------------------------------------------------------- */
/* High level matrix inverse and product routines in lp_solve                         */
/* ---------------------------------------------------------------------------------- */

/* ---------------------------------------------------------------------------------- */
/* Inverse handling                                                                   */
/* ---------------------------------------------------------------------------------- */
/*
        A brief description of the inversion logic in lp_solve
   -----------------------------------------------------------------

   In order to better understand the code in lp_solve in relation to
   standard matrix-based textbook descriptions, I (KE) will briefly
   explain the conventions and associated matrix algebra.  The matrix
   description of a linear program (as represented by lp_solve) goes
   like this:

           maximize         c'x
           subject to  r <=  Ax <= b
           where       l <=   x <= u

   The matrix A is partitioned into two column sets [B|N], where B is
   a square matrix of "basis" variables containing non-fixed
   variables of the linear program at any given stage and N is the
   submatrix of corresponding non-basic, fixed variables. The
   variables (columns) in N may be fixed at their lower or upper levels.

   Similarly, the c vector is partitioned into the basic and non-basic
   parts [z|n].

   lp_solve stores the objective vector c and the data matrix A in a
   common sparse format where c is put as the 0-th row of A. This may
   be called the "A~" form:

                       A~ = [ c ]
                            [ A ]

   Linear programming involves solving linear equations based on B,
   various updates and bookkeeping operations.  The relationship
   between the common storage of c and A (e.g. A~) vs. the inverse of
   B therefore needs to be understood.  In lp_solve, B is stored in
   an expanded, bordered format using the following (non-singular)
   representation:

                       B~ = [ 1 z ]
                            [ 0 B ]

   Any inversion engine used by lp_solve must therefore explicitly
   represent and handle the implications of this structure for
   associated matrix multiplications.

   The standard matrix formula for computing the inverse of a bordered
   matrix shows what the inversion of B~ actually produces:

                  Inv(B~) = [ 1 -z*Inv(B) ]
                            [ 0   Inv(B)  ]

   The A~ and B~ representations mean that it becomes necessary to be
   aware of the side effects of the presence of the top row when doing
   product operations such as b'N, btran and ftran.  A very nice
   thing about the matrix representation in lp_solve is that a very
   common update in the simplex algorithm (reduced costs) is obtained
   simply by setting 1 at the top of the vector being pre-multiplied
   with Inv(B~).

   However, if the objective vector (c) is changed, the representation
   requires B / B~ to be reinverted.  Also, when doing ftran, btran
   and x'A-type operations, you will patently get the incorrect result
   if you simply copy the operations described in textbooks.  First I'll
   show the results of an ftran operation:

                   Bx = a  ==>  x = ftran(a)

   In lp_solve, this operation solves:

                   [ 1 z ] [y] = [d]
                   [ 0 B ] [x]   [a]

   Using the Inv(B~) expression earlier, the ftran result is:

             [y] = [ 1 -z*Inv(B) ] [d] = [ d - z*Inv(B)*a ]
             [x]   [ 0   Inv(B)  ] [a]   [   Inv(B)*a     ]

   Similarily, doing the left solve - performing the btran calculation:

                   [x y] [ 1 z ] = [d a']
                         [ 0 B ]

   ... will produce the following result in lp_solve:

   [x y] = [d a'] [ 1 -z*Inv(B) ] = [ d | -d*z*Inv(B) + a'*Inv(B) ]
                  [ 0   Inv(B)  ]

   So, if you thought you were computing "a'*Inv(B)", look again.
   In order to produce the desired result, you have to set d to 0
   before the btran operation.

   Equipped with this understanding, I hope that you see that
   the approach in lp_solve is actually pretty convenient.  It
   also makes it far easier to extend functionality by drawing on
   formulas and expressions from LP literature that assume the
   conventional syntax and representation.

                                   Kjell Eikland -- November 2003.

*/

STATIC MYBOOL __WINAPI invert(lprec *lp, MYBOOL shiftbounds, MYBOOL final)
{
  MYBOOL *usedpos, resetbasis;
  REAL   test;
  int    k, i, j;
  int    singularities, usercolB;

 /* Make sure the tags are correct */
  if(!mat_validate(lp->matA)) {
    lp->spx_status = INFEASIBLE;
    return(FALSE);
  }

 /* Create the inverse management object at the first call to invert() */
  if(lp->invB == NULL)
    lp->bfp_init(lp, lp->rows, 0);
  else
    lp->bfp_preparefactorization(lp);
  singularities = 0;

 /* Must save spx_status since it is used to carry information about
    the presence and handling of singular columns in the matrix */
  if(userabort(lp, MSG_INVERT))
    return(FALSE);

#ifdef Paranoia
  if(lp->spx_trace)
    report(lp, DETAILED, "invert: iteration %6d, inv-length %4d, OF " MPSVALUEMASK " \n",
                         get_total_iter(lp), lp->bfp_colcount(lp), (double) -lp->rhs[0]);
#endif

 /* Store state of pre-existing basis, and at the same time check if
    the basis is I; in this case take the easy way out */
  if(!allocMYBOOL(lp, &usedpos, lp->sum + 1, TRUE)) {
    lp->spx_status = NOMEMORY;
    lp->bb_break = TRUE;
    return(FALSE);
  }
  usedpos[0] = TRUE;
  usercolB = 0;
  for(i = 1; i <= lp->rows; i++) {
    k = lp->var_basic[i];
    if(k > lp->rows)
      usercolB++;
    usedpos[k] = TRUE;
  }
#ifdef Paranoia
  if(!verifyBasis(lp))
    report(lp, SEVERE, "invert: Invalid basis detected (iteration %d).\n",
                       get_total_iter(lp));
#endif

 /* Tally matrix nz-counts and check if we should reset basis
    indicators to all slacks */
  resetbasis = (MYBOOL) ((usercolB > 0) && lp->bfp_canresetbasis(lp));
  k = 0;
  for(i = 1; i <= lp->rows; i++) {
    if(lp->var_basic[i] > lp->rows)
      k += mat_collength(lp->matA, lp->var_basic[i] - lp->rows);
    if(resetbasis) {
      j = lp->var_basic[i];
      if(j > lp->rows)
        lp->is_basic[j] = FALSE;
      lp->var_basic[i] = i;
      lp->is_basic[i] = TRUE;
    }
  }
 /* If we don't reset the basis, optionally sort it by index.
    (experimental code at one time used for debugging) */
#if 0
  if(!resetbasis) {
    j = lp->rows;
    for(i = 1; i <= lp->rows; i++) {
      if(lp->is_basic[i])
        lp->var_basic[i] = i;
      else {
        j++;
        while((j <= lp->sum) && !lp->is_basic[j])
          j++;
        lp->var_basic[i] = j;
      }
    }
  }
#endif

 /* Now do the refactorization */
  singularities = lp->bfp_factorize(lp, usercolB, k, usedpos, final);

 /* Do user reporting */
  if(userabort(lp, MSG_INVERT))
    goto Cleanup;

#ifdef Paranoia
  if(lp->spx_trace) {
    k = lp->bfp_nonzeros(lp, FALSE);
    report(lp, DETAILED, "invert: inv-length %4d, inv-NZ %4d\n",
                         lp->bfp_colcount(lp), k);
  }
#endif

 /* Finalize factorization/inversion */
  lp->bfp_finishfactorization(lp);

  /* Recompute the RHS ( Ref. lp_solve inverse logic and Chvatal p. 121 ) */
#ifdef DebugInv
  blockWriteLREAL(stdout, "RHS-values pre invert", lp->rhs, 0, lp->rows);
#endif
  recompute_solution(lp, shiftbounds);
  restartPricer(lp, AUTOMATIC);
#ifdef DebugInv
  blockWriteLREAL(stdout, "RHS-values post invert", lp->rhs, 0, lp->rows);
#endif

Cleanup:
  /* Check for numerical instability indicated by frequent refactorizations */
  test = get_refactfrequency(lp, FALSE);
  if(test < MIN_REFACTFREQUENCY) {
    test = get_refactfrequency(lp, TRUE);
    report(lp, NORMAL, "invert: Refactorization frequency %.1g indicates numeric instability\n",
                       test);
    lp->spx_status = NUMFAILURE;
  }

  FREE(usedpos);
  return((MYBOOL) (singularities <= 0));
} /* invert */


STATIC MYBOOL fimprove(lprec *lp, REAL *pcol, int *nzidx, REAL roundzero)
{
  REAL   *errors, sdp;
  int    j;
  MYBOOL Ok = TRUE;

  allocREAL(lp, &errors, lp->rows + 1, FALSE);
  if(errors == NULL) {
    lp->spx_status = NOMEMORY;
    Ok = FALSE;
    return(Ok);
  }
  MEMCOPY(errors, pcol, lp->rows + 1);
  lp->bfp_ftran_normal(lp, pcol, nzidx);
  prod_Ax(lp, SCAN_ALLVARS, pcol, NULL, 0, 0.0, -1, errors, NULL);
  lp->bfp_ftran_normal(lp, errors, NULL);

  sdp = 0;
  for(j = 1; j <= lp->rows; j++)
    if(fabs(errors[j])>sdp)
      sdp = fabs(errors[j]);
  if(sdp > lp->epsmachine) {
    report(lp, DETAILED, "Iterative FTRAN correction metric %g", sdp);
    for(j = 1; j <= lp->rows; j++) {
      pcol[j] += errors[j];
      my_roundzero(pcol[j], roundzero);
    }
  }
  FREE(errors);
  return(Ok);
}

STATIC MYBOOL bimprove(lprec *lp, REAL *rhsvector, int *nzidx, REAL roundzero)
{
  int    j;
  REAL   *errors, err, maxerr;
  MYBOOL Ok = TRUE;

  allocREAL(lp, &errors, lp->sum + 1, FALSE);
  if(errors == NULL) {
    lp->spx_status = NOMEMORY;
    Ok = FALSE;
    return(Ok);
  }
  MEMCOPY(errors, rhsvector, lp->sum + 1);

  /* Solve Ax=b for x, compute b back */
  lp->bfp_btran_normal(lp, errors, nzidx);
  prod_xA(lp, SCAN_ALLVARS+USE_BASICVARS, errors, NULL, XRESULT_FREE, 0.0, 1.0,
                                          errors, NULL);

  /* Take difference with ingoing values, while shifting the column values
     to the rows section and zeroing the columns again */
  for(j = 1; j <= lp->rows; j++)
    errors[j] = errors[lp->rows+lp->var_basic[j]] - rhsvector[j];
  for(j = lp->rows; j <= lp->sum; j++)
    errors[j] = 0;

  /* Solve the b errors for the iterative x adjustment */
  lp->bfp_btran_normal(lp, errors, NULL);

  /* Generate the adjustments and compute statistic */
  maxerr = 0;
  for(j = 1; j <= lp->rows; j++) {
    if(lp->var_basic[j]<=lp->rows) continue;
    err = errors[lp->rows+lp->var_basic[j]];
    if(fabs(err)>maxerr)
      maxerr = fabs(err);
  }
  if(maxerr > lp->epsmachine) {
    report(lp, DETAILED, "Iterative BTRAN correction metric %g", maxerr);
    for(j = 1; j <= lp->rows; j++) {
      if(lp->var_basic[j]<=lp->rows) continue;
      rhsvector[j] += errors[lp->rows+lp->var_basic[j]];
      my_roundzero(rhsvector[j], roundzero);
    }
  }
  FREE(errors);
  return(Ok);
}

STATIC void ftran(lprec *lp, REAL *rhsvector, int *nzidx, REAL roundzero)
{
  if((lp->improve & IMPROVE_FTRAN) && lp->bfp_pivotcount(lp))
    fimprove(lp, rhsvector, nzidx, roundzero);
  else {
    lp->bfp_ftran_normal(lp, rhsvector, nzidx);
  }
}

STATIC void btran(lprec *lp, REAL *rhsvector, int *nzidx, REAL roundzero)
{
  if((lp->improve & IMPROVE_BTRAN) && lp->bfp_pivotcount(lp))
    bimprove(lp, rhsvector, nzidx, roundzero);
  else {
    lp->bfp_btran_normal(lp, rhsvector, nzidx);
  }
}

STATIC MYBOOL fsolve(lprec *lp, int varin, REAL *pcol, int *nzidx, REAL roundzero, REAL ofscalar, MYBOOL prepareupdate)
/* Was setpivcol in versions earlier than 4.0.1.8 - KE */
{
  MYBOOL ok = TRUE;

  if(varin > 0)
    obtain_column(lp, varin, pcol, nzidx, NULL);

 /* Solve, adjusted for objective function scalar */
  pcol[0] *= ofscalar;
  if(prepareupdate)
    lp->bfp_ftran_prepare(lp, pcol, nzidx);
  else
    ftran(lp, pcol, nzidx, roundzero);

  return(ok);

} /* fsolve */


STATIC MYBOOL bsolve(lprec *lp, int row_nr, REAL *rhsvector, int *nzidx, REAL roundzero, REAL ofscalar)
{
  MYBOOL ok = TRUE;

  if(row_nr >= 0)
    row_nr = obtain_column(lp, row_nr, rhsvector, nzidx, NULL);

  /* Solve, adjusted for objective function scalar */
  rhsvector[0] *= ofscalar;
  btran(lp, rhsvector, nzidx, roundzero);

  return(ok);

} /* bsolve */


STATIC int prod_Ax(lprec *lp, int varset, REAL *input, int *nzinput,
                              int range, REAL roundzero, REAL multfactor,
                                          REAL *output, int *nzoutput)
/* prod_Ax is only used in fimprove; note that it is NOT VALIDATED/verified as of 20030801 - KE */
{
  int      j, colnr, Extrap,
           ib, ie, vb, ve;
  MYBOOL   omitfixed, omitnonfixed;
  MATrec   *mat = lp->matA;
  REAL     sdp;
  REAL     *value;
  int      *rownr;


  /* Find what variable range to scan (default is all) */
  /* First determine the starting position; add from the top, going down */
  Extrap = abs(lp->Extrap);
  vb = 1;
  if(varset & SCAN_ARTIFICIALVARS)
   vb = lp->sum - Extrap + 1;
  if(varset & SCAN_USERVARS)
    vb = lp->rows + 1;
  if(varset & SCAN_SLACKVARS)
    vb = 1;

  /* Then determine the ending position, add from the bottom, going up */
  ve = lp->sum;
  if(varset & SCAN_SLACKVARS)
    ve = lp->rows;
  if(varset & SCAN_USERVARS)
    ve = lp->sum - Extrap;
  if(varset & SCAN_ARTIFICIALVARS)
    ve = lp->sum;

  /* Determine exclusion columns */
  omitfixed = (MYBOOL) ((varset & OMIT_FIXED) != 0);
  omitnonfixed = (MYBOOL) ((varset & OMIT_NONFIXED) != 0);
  if(omitfixed && omitnonfixed)
    return(FALSE);

  /* Scan the basic columns */
  for(j = 1; j <= lp->rows; j++) {
    colnr = lp->var_basic[j];

    /* Do exclusions */
    if(colnr < vb || colnr > ve)
      continue;

    sdp = lp->upbo[colnr];
    if((omitfixed && (sdp == 0)) ||
       (omitnonfixed && (sdp != 0)))
      continue;

    /* Exclusions done, perform the multiplication */
    sdp = multfactor*input[j];
    if(colnr <= lp->rows)               /* A slack variable is in the basis */
      output[colnr] += sdp;
    else {                              /* A normal variable is in the basis */
      colnr -= lp->rows;
      ib = mat->col_end[colnr - 1];
      ie = mat->col_end[colnr];
      rownr = &COL_MAT_ROWNR(ib);
      value = &COL_MAT_VALUE(ib);
      for(; ib < ie;
          ib++, rownr += matRowColStep, value += matValueStep) {
        output[*rownr] += (*value)*sdp;
      }
    }
  }
  roundVector(output+1, lp->rows-1, roundzero);

  return(TRUE);
}

STATIC int prod_xA(lprec *lp, int varset, REAL *input, int *nzinput,
                              int range,  REAL roundzero, REAL ofscalar,
                                          REAL *output, int *nzoutput)
/* Note that the dot product xa is stored at the active column index of A, i.e. of a.
   This means that if the basis only contains non-slack variables, output may point to
   the same vector as input, without overwriting the [0..rows] elements. */
{
  int      i, rownr, varnr, Extrap, vb, ve, ib, ie;
  MYBOOL   omitfixed, omitnonfixed;
  RREAL    vmax;
  register RREAL v;
  REAL     value, *rhsvector;
  int      inz, *rowin, countNZ = 0;
  MATrec   *mat = lp->matA;
  register REAL     *matValue;
  register int      *matRownr;

  /* Clean output area (only necessary if we are returning the full vector) */
  if(nzoutput == NULL) {
    if(input == output)
      MEMCLEAR(output+lp->rows+1, lp->columns);
    else
      MEMCLEAR(output, lp->sum+1);
  }

  /* Find what variable range to scan - default is {SCAN_USERVARS} */
  /* First determine the starting position; add from the top, going down */
  Extrap = abs(lp->Extrap);
  vb = lp->rows + 1;
  if(varset & SCAN_ARTIFICIALVARS)
    vb = lp->sum - Extrap + 1;
  if(varset & SCAN_USERVARS)
    vb = lp->rows + 1;
  if(varset & SCAN_SLACKVARS)
    vb = 1;

  /* Then determine the ending position, add from the bottom, going up */
  ve = lp->sum;
  if(varset & SCAN_SLACKVARS)
    ve = lp->rows;
  if(varset & SCAN_USERVARS)
    ve = lp->sum - Extrap;
  if(varset & SCAN_ARTIFICIALVARS)
    ve = lp->sum;

  /* Adjust for partial pricing */
  if(varset & SCAN_PARTIALBLOCK) {
    if(vb > lp->rows)
      vb = MAX(vb, partial_blockStart(lp, FALSE));
    ve = MIN(ve, partial_blockEnd(lp, FALSE));
  }

  /* Determine exclusion columns */
  omitfixed = (MYBOOL) ((varset & OMIT_FIXED) != 0);
  omitnonfixed = (MYBOOL) ((varset & OMIT_NONFIXED) != 0);
  if(omitfixed && omitnonfixed)
    return(FALSE);

  /* Scan the target colums */
  vmax = 0;
  for(varnr = vb, rhsvector = output + vb;
      varnr <= ve; varnr++, rhsvector++) {

    /* Skip gap in the specified column scan range (possibly user variables) */
    if((varnr > lp->rows) && (varnr <= lp->sum-Extrap) && !(varset & SCAN_USERVARS))
      continue;

    /* Find if the variable is in the scope - default is {Ø} */
    i = lp->is_basic[varnr];
    if((varset & USE_BASICVARS) > 0 && (i))
      ;
    else if((varset & USE_NONBASICVARS) > 0 && (!i))
      ;
    else
      continue;

    v = lp->upbo[varnr];
    if((omitfixed && (v == 0)) ||
       (omitnonfixed && (v != 0)))
      continue;

    if(varnr <= lp->rows) {
      v = input[varnr];
    }
    else {
      i = varnr - lp->rows;
      v = 0;
      ie = mat->col_end[i];
      ib = mat->col_end[i - 1];
      if(ib < ie) {

        /* Do dense input vector version */
        if(nzinput == NULL) {

          /* Handle the objective row specially */
          matRownr = &COL_MAT_ROWNR(ib);
          matValue = &COL_MAT_VALUE(ib);
          if(*matRownr == 0) {
            value = *matValue;
            ib++;
            matValue += matValueStep;
            matRownr += matRowColStep;
          }
          else {
            value = 0;
          }
          if(modifyOF1(lp, varnr, &value, ofscalar))
            v += input[0] * value;

#define NoLoopUnroll
#ifdef NoLoopUnroll
          /* Then loop over all regular rows */
          if(ib < ie)
          for(; ib < ie; ib++) {
            v += input[*matRownr] * (*matValue);
            matValue += matValueStep;
            matRownr += matRowColStep;
          }
#else
          /* Prepare for simple loop unrolling */
          if(((ie-ib) % 2) == 1) {
            v += input[*matRownr] * (*matValue);
            ib++;
            matValue += matValueStep;
            matRownr += matRowColStep;
          }

          /* Then loop over remaining pairs of regular rows */
          while(ib < ie) {
            v += input[*matRownr] * (*matValue) +
                 input[*(matRownr+matRowColStep)] * (*(matValue+matValueStep));
            ib += 2;
            matValue += 2*matValueStep;
            matRownr += 2*matRowColStep;
          }
#endif
        }
        /* Do sparse input vector version */
        else {

          /* Initialize pointers */
          inz = 1;
          rowin = nzinput+inz;
          matRownr = &COL_MAT_ROWNR(ib);
          matValue = &COL_MAT_VALUE(ib);
          ie--;

          /* Handle OF row separately (since it can be overridden) */
          if(*rowin == 0) {
            if(*matRownr == 0) {
              ib++;
              /* Step forward at right */
              if(ib <= ie) {
                value = *matValue;
                matValue += matValueStep;
                matRownr += matRowColStep;
              }
            }
            else
              value = 0;
            if(modifyOF1(lp, varnr, &value, ofscalar))
              v += input[0]*value;
            /* Step forward at left */
            inz++;
            rowin++;
          }

          /* Then loop over all non-OF rows */
          while((inz <= *nzinput) && (ib <= ie)) {

           /* Try to synchronize at right */
            while((*rowin > *matRownr) && (ib < ie)) {
              ib++;
              matValue += matValueStep;
              matRownr += matRowColStep;
            }
            /* Try to synchronize at left */
            while((*rowin < *matRownr) && (inz < *nzinput)) {
              inz++;
              rowin++;
            }
            /* Perform dot product operation if there was a match */
            if(*rowin == *matRownr) {
              v += input[*rowin] * (*matValue);
              /* Step forward at left */
              inz++;
              rowin++;
            }
          }
        }
      }
#ifndef Prod_xA_RoundRelative
      my_roundzero(v, roundzero);
#endif
    }
    vmax = MAX(vmax, fabs((REAL) v));
    if((v != 0) && ((range == XRESULT_FREE) || (v*range > 0))) {
      countNZ++;
      if(nzoutput != NULL)
        nzoutput[countNZ] = varnr;
    }
    *rhsvector = (REAL) v;
  }

  /* Check if we should do relative rounding */
#ifdef Prod_xA_RoundRelative
  if((roundzero > 0) && (nzoutput != NULL)) {
    ie = 0;
    for(i = 1; i <= countNZ;  i++) {
      rownr = nzoutput[i];
      if(fabs(output[rownr])/vmax < roundzero)
        output[rownr] = 0;
      else if((range == XRESULT_FREE) || (output[rownr]*range > 0)) {
        ie++;
        nzoutput[ie] = rownr;
      }
    }
    countNZ = ie;
  }
#endif

  if(nzoutput != NULL)
    *nzoutput = countNZ;
  return(countNZ);
}

STATIC void prod_xA2(lprec *lp, REAL *prow, int prange, REAL proundzero, int *nzprow,
                                REAL *drow, int drange, REAL droundzero, int *nzdrow, REAL ofscalar)
{
#ifdef UseDouble_prod_xA
  prod_xA(lp, SCAN_ALLVARS+
              SCAN_PARTIALBLOCK+
              USE_NONBASICVARS+OMIT_FIXED,
              prow, NULL, XRESULT_RC, proundzero, ofscalar,
              prow, nzprow);
  prod_xA(lp, SCAN_ALLVARS+
              SCAN_PARTIALBLOCK+
              USE_NONBASICVARS+OMIT_FIXED,
              drow, NULL, XRESULT_FREE, droundzero, ofscalar,
              drow, nzdrow);
#else
  int      i, ii, varnr, rownr, ib, ie;
  RREAL    dmax, pmax;
  register RREAL d, p;
  MATrec   *mat = lp->matA;
  REAL     value;
  register REAL     *matValue;
  register int      *matRownr;

  /* Loop over slack variables only if we do sparse mapping; since slack variable
     coefficients are all 1, values in the incoming vectors are preserved */
  pmax = 0;
  dmax = 0;
  i  = partial_blockStart(lp, FALSE);
  ii = partial_blockEnd(lp, FALSE);
  if(ii == lp->sum)
    ii -= abs(lp->Extrap);
  if((nzprow != NULL) || (nzdrow != NULL)) {
    if(nzprow != NULL)
      *nzprow = 0;
    if(nzdrow != NULL)
      *nzdrow = 0;
    for(; i <= lp->rows; i++)
      if(!lp->is_basic[i] && (lp->upbo[i] > 0)) {
        pmax = MAX(pmax, fabs(prow[i]));
        dmax = MAX(dmax, fabs(drow[i]));
        if((nzprow != NULL) && (prow[i] != 0) &&
           ((prange == XRESULT_FREE) || (prow[i]*my_chsign(lp->is_lower[i], -1) > 0))) {
          (*nzprow)++;
          nzprow[*nzprow] = i;
        }
        if((nzdrow != NULL) && (drow[i] != 0) &&
           ((drange == XRESULT_FREE) || (drow[i]*my_chsign(lp->is_lower[i], -1) > 0))) {
          (*nzdrow)++;
          nzdrow[*nzdrow] = i;
        }
      }
  }

  /* Loop over user variables or any active subset */
  i = MAX(i, lp->rows + 1);
  for(; i <= ii; i++) {

    varnr = i - lp->rows;

    /* Only non-basic user variables */
    if(lp->is_basic[i])
      continue;

    /* Use only non-fixed variables */
    if(lp->upbo[i] == 0)
      continue;

    p = 0;
    d = 0;
    ie = mat->col_end[varnr];
    ib = mat->col_end[varnr - 1];

    if(ib < ie) {

      /* Handle the objective row specially */
      matRownr = &COL_MAT_ROWNR(ib);
      matValue = &COL_MAT_VALUE(ib);
      if(*matRownr == 0) {
        value = *matValue;
        ib++;
        matValue += matValueStep;
        matRownr += matRowColStep;
      }
      else {
        value = 0;
      }
      if(modifyOF1(lp, i, &value, ofscalar)) {
        p += prow[0] * value;
        d += drow[0] * value;
      }
#ifdef NoLoopUnroll
      /* Then loop over all regular rows */
      if(ib < ie)
      for( ; ib < ie; ib++) {
        p += prow[*matRownr] * (*matValue);
        d += drow[*matRownr] * (*matValue);
        matValue += matValueStep;
        matRownr += matRowColStep;
      }
#else
      /* Prepare for simple loop unrolling */
      if(((ie-ib) % 2) == 1) {
        p += prow[*matRownr] * (*matValue);
        d += drow[*matRownr] * (*matValue);
        ib++;
        matValue += matValueStep;
        matRownr += matRowColStep;
      }

      /* Then loop over remaining pairs of regular rows */
      while(ib < ie) {
        p += prow[*matRownr] * (*matValue) +
             prow[*(matRownr+matRowColStep)] * (*(matValue+matValueStep));
        d += drow[*matRownr] * (*matValue) +
             drow[*(matRownr+matRowColStep)] * (*(matValue+matValueStep));
        ib += 2;
        matValue += 2*matValueStep;
        matRownr += 2*matRowColStep;
      }
#endif

    }

    my_roundzero(p, proundzero);
    pmax = MAX(pmax, fabs((REAL) p));
    prow[i] = (REAL) p;
    if((nzprow != NULL) && (p != 0) &&
       ((prange == XRESULT_FREE) || (p*my_chsign(lp->is_lower[i], -1) > 0))) {
      (*nzprow)++;
      nzprow[*nzprow] = i;
    }
    my_roundzero(d, droundzero);
    dmax = MAX(dmax, fabs((REAL) d));
    drow[i] = (REAL) d;
    if((nzdrow != NULL) && (d != 0) &&
       ((drange == XRESULT_FREE) || (d*my_chsign(lp->is_lower[i], -1) > 0))) {
      (*nzdrow)++;
      nzdrow[*nzdrow] = i;
    }
  }

  /* Check if we should do relative rounding */
#ifdef Prod_xA_RoundRelative
  if((proundzero > 0) && (nzprow != NULL)) {
    ie = 0;
    for(i = 1; i <= *nzprow;  i++) {
      rownr = nzprow[i];
      if(fabs(prow[rownr])/pmax < proundzero)
        prow[rownr] = 0;
      else if((prange == XRESULT_FREE) || (prow[rownr]*my_chsign(lp->is_lower[rownr], -1) > 0)) {
        ie++;
        nzprow[ie] = rownr;
      }
    }
    *nzprow = ie;
  }
  if((droundzero > 0) && (nzdrow != NULL)) {
    ie = 0;
    for(i = 1; i <= *nzdrow;  i++) {
      rownr = nzdrow[i];
      if(fabs(drow[rownr])/dmax < droundzero)
        drow[rownr] = 0;
      else if((drange == XRESULT_FREE) || (drow[rownr]*my_chsign(lp->is_lower[rownr], -1) > 0)) {
        ie++;
        nzdrow[ie] = rownr;
      }
    }
    *nzdrow = ie;
  }
#endif

#endif

}

STATIC void bsolve_xA2(lprec *lp, int row_nr1, REAL *vector1, REAL roundzero1, int *nzvector1,
                                  int row_nr2, REAL *vector2, REAL roundzero2, int *nzvector2)
{
  REAL ofscalar = 1.0;
  int  rc_range = XRESULT_RC;

 /* Clear and initialize first vector */
  if(nzvector1 == NULL)
    MEMCLEAR(vector1, lp->sum + 1);
  else
    MEMCLEAR(vector1, lp->rows + 1);
  vector1[row_nr1] = 1;
/*  workINT[0] = 1;
  workINT[1] = row_nr1; */

  if(vector2 == NULL) {
    lp->bfp_btran_normal(lp, vector1, NULL);
    prod_xA(lp, SCAN_USERVARS+
                SCAN_PARTIALBLOCK+
                USE_NONBASICVARS+OMIT_FIXED,
                vector1, NULL, rc_range, roundzero1, ofscalar*0,
                vector1, nzvector1);
  }
  else {

   /* Clear and initialize second vector */
    if(nzvector2 == NULL)
      MEMCLEAR(vector2, lp->sum + 1);
    else
      MEMCLEAR(vector2, lp->rows + 1);
    vector2[row_nr2] = 1;
/*    workINT[2] = 1;
    workINT[3] = row_nr2; */

   /* A double BTRAN equation solver process is implemented "in-line" below in
      order to save time and to implement different rounding for the two */
    lp->bfp_btran_double(lp, vector1, NULL, vector2, NULL);

   /* Multiply solution vectors with matrix values */
    prod_xA2(lp, vector1, rc_range,     roundzero1, nzvector1,
                 vector2, XRESULT_FREE, roundzero2, nzvector2, ofscalar);
  }
}

