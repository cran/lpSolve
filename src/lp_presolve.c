
/* -------------------------------------------------------------------------
   Presolve routines for lp_solve v5.0+
   -------------------------------------------------------------------------
    Author:        Kjell Eikland
    Contact:       kjell.eikland@broadpark.no
    License terms: LGPL.

    Requires:      lp_lib.h, lp_presolve, lp_crash.h, lp_scale.h

    Release notes:
    v5.0.0  1 January 2004      Significantly expanded and repackaged
                                presolve routines.
    v5.0.1  1 April   2004      Added reference to new crash module
    v5.1.0  20 August 2004      Reworked infeasibility detection.
                                Added encapsulation of presolve undo logic.
    v5.1.1  10 September 2004   Added variable bound tightening based on 
                                full-constraint information, as well as 
                                variable fixing by duality.

   ------------------------------------------------------------------------- */

#include <string.h>
#include "commonlib.h"
#include "lp_lib.h"
#include "lp_presolve.h"
#include "lp_crash.h"
#include "lp_scale.h"
#include "lp_report.h"

#ifdef FORTIFY
# include "lp_fortify.h"
#endif


#define DoPresolveRounding
/*#define DoPresolveRelativeTest*/

#define AggressiveRowPresolve        /* Extra row presolve (beware of accuracy) */
/*#define ApplyDualPresolve */           /* Use dual information to fix variables */
/*#define UseFullConstraintInfo */   /* Use bound sums for non-singular constraints */

#define DeferredCompactRowA          /* Defer matrix row compacting to speed presolve */
#define DeferredCompactColA          /* Defer matrix column compacting to speed presolve */

#define PRESOLVE_EPSVALUE  lp->epsprimal
/*#define PRESOLVE_EPSROUND  lp->epsprimal*/
#define PRESOLVE_EPSROUND  lp->epsint 

STATIC MYBOOL presolve_createUndo(lprec *lp)
{
  if(lp->presolve_undo != NULL)
    presolve_freeUndo(lp);
  lp->presolve_undo = (presolveundorec *) calloc(1, sizeof(presolveundorec));
  if(lp->presolve_undo == NULL)
    return( FALSE );
  return( TRUE );
}
STATIC MYBOOL inc_presolve_space(lprec *lp, int delta, MYBOOL isrows)
{
  int i, ii, oldrowcolalloc, rowcolsum, oldrowalloc;
  presolveundorec *psundo = lp->presolve_undo;
  
  if(psundo == NULL) {
    presolve_createUndo(lp);
    psundo = lp->presolve_undo;
  }

  /* Set constants */
  oldrowalloc = lp->rows_alloc-delta;
  oldrowcolalloc = lp->sum_alloc;
  lp->sum_alloc += delta;
  rowcolsum = lp->sum_alloc + 1;

  /* Reallocate lp memory */
  allocREAL(lp, &psundo->fixed_rhs, lp->rows_alloc+1, AUTOMATIC);
  allocINT(lp, &psundo->var_to_orig, rowcolsum, AUTOMATIC);
  allocINT(lp, &psundo->orig_to_var, rowcolsum, AUTOMATIC);

  /* Fill in default values, where appropriate */
  for(i = oldrowcolalloc+1, ii = oldrowalloc+1; i < rowcolsum; i++, ii++) {
    psundo->var_to_orig[i] = 0;
    psundo->orig_to_var[i] = 0;
    if(isrows)
      psundo->fixed_rhs[ii] = 0;
  }

  return(TRUE);
}
STATIC MYBOOL presolve_setOrig(lprec *lp, int orig_rows, int orig_cols)
{
  presolveundorec *psundo = lp->presolve_undo;
  
  if(psundo == NULL)
    return( FALSE );
  psundo->orig_rows = orig_rows;
  psundo->orig_columns = orig_cols;
  if(lp->wasPresolved)
    presolve_fillUndo(lp, orig_rows, orig_cols, FALSE);
  return( TRUE );
}
STATIC MYBOOL presolve_fillUndo(lprec *lp, int orig_rows, int orig_cols, MYBOOL setOrig)
{
  int i;
  presolveundorec *psundo = lp->presolve_undo;

  for(i = 0; i <= orig_rows; i++) {
    psundo->var_to_orig[i] = i;
    psundo->orig_to_var[i] = i;
    psundo->fixed_rhs[i]   = 0;
  }
  for(i = 1; i <= orig_cols; i++) {
    psundo->var_to_orig[orig_rows + i] = i;
    psundo->orig_to_var[orig_rows + i] = i;
  }
  if(setOrig)
    presolve_setOrig(lp, orig_rows, orig_cols);
  
  return( TRUE );
}
STATIC MYBOOL presolve_freeUndo(lprec *lp)
{
  presolveundorec *psundo = lp->presolve_undo; 
  
  if(psundo == NULL)
    return( FALSE );
  FREE(psundo->orig_to_var);
  FREE(psundo->var_to_orig);
  FREE(psundo->fixed_rhs);
  FREE(lp->presolve_undo);
  return( TRUE );
}

STATIC REAL sumplumin(lprec *lp, int item, REAL *plu, REAL *neg)
{
  if(fabs(plu[item]) >= lp->infinite)
    return( plu[item] );
  else if(fabs(neg[item]) >= lp->infinite)
    return( neg[item] );
  else
    return( plu[item]+neg[item] );
}

STATIC REAL presolve_round(lprec *lp, REAL value, MYBOOL isGE)
{
  REAL epsvalue = PRESOLVE_EPSROUND;

#ifdef DoPresolveRounding
  value += my_chsign(isGE, epsvalue/SCALEDINTFIXRANGE);
  value = restoreINT(value, epsvalue);
#endif
  return( value );
}
STATIC REAL presolve_precision(lprec *lp, REAL value)
{
  REAL epsvalue = PRESOLVE_EPSROUND;

#ifdef DoPresolveRounding
  value = restoreINT(value, epsvalue);
#endif
  return( value );
}

STATIC int presolve_validate(lprec *lp, LLrec *rowmap, LLrec *colmap)
{
  int i, j, errc = 0;

  /* Validate constraint bounds */
  for(i = 1; i < lp->rows; i++) {
    if((rowmap != NULL) && !isActiveLink(rowmap, i))
      continue;
    /* Check if we have a negative range */
    if(lp->orig_upbo[i] < 0) {
      errc++;
      report(lp, SEVERE, "presolve_validate: Detected negative range %g for row %d\n",
                         lp->orig_upbo[i], i);
    }
  }
  /* Validate variables */
  for(j = 1; j < lp->columns; j++) {
    if((colmap != NULL) && !isActiveLink(colmap, j))
      continue;
    i = lp->rows+j;
    /* Check if we have infeasible  bounds */
    if(lp->orig_lowbo[i] > lp->orig_upbo[i]) {
      errc++;
      report(lp, SEVERE, "presolve_validate: Detected UB < LB for column %d\n",
                         j);
    }
  }
  /* Return total number of errors */
  return( errc );
}


STATIC MYBOOL presolve_singletonbounds(lprec *lp, int rownr, int colnr, REAL *lobound, REAL *upbound, REAL *aval)
{
  REAL   coeff_a, epsvalue = PRESOLVE_EPSVALUE;
  MYBOOL isneg;

  /* Compute row singleton variable range */
  if(is_constr_type(lp, rownr, EQ) && (fabs(*lobound) < epsvalue))
    *lobound = *upbound = 0;
  else {
    if(aval == NULL)
      coeff_a = get_mat(lp, rownr, colnr);
    else
      coeff_a = *aval;
    isneg = (MYBOOL) (coeff_a < 0);
    if(*lobound > -lp->infinite)
      *lobound /= coeff_a;
    else if(isneg)
      *lobound = -(*lobound);
    if(*upbound < lp->infinite)
      *upbound /= coeff_a;
    else if(isneg)
      *upbound = -(*upbound);
    if(isneg)
      swapREAL(lobound, upbound);
  }
  
  /* Check against bound - handle SC variables specially */
  if(is_semicont(lp, colnr)) {
    coeff_a = get_lowbo(lp, colnr);
    if(coeff_a > 0) {
      *lobound = MAX(*lobound, 0.0);
      *upbound = MIN(*upbound, get_upbo(lp, colnr));
    }
    else {
      coeff_a = get_upbo(lp, colnr);
      if(coeff_a > 0) {
        *lobound = MAX(*lobound, get_lowbo(lp, colnr));
        *upbound = MIN(*upbound, 0.0);
      }
    }
  }
  else {
    *lobound = MAX(*lobound, get_lowbo(lp, colnr));
    *upbound = MIN(*upbound, get_upbo(lp, colnr));
    /* Attempt numeric error management / correction */
#if 0
    if(!is_infinite(lp, *lobound) && !is_infinite(lp, *upbound)) {
      coeff_a = ((*lobound) + (*upbound))*0.5;
      if((*upbound < coeff_a) &&
         (my_reldiff(*upbound, coeff_a) > -epsvalue)) {
        my_roundzero(coeff_a, epsvalue);
        *lobound = coeff_a;
        *upbound = coeff_a;
      }
    }
#endif
  }
  
  /* Return with consistency status */
#ifdef DoPresolveRelativeTest
  return( (MYBOOL) (my_reldiff(*upbound, *lobound) >= - epsvalue) );
#else  
  return( (MYBOOL) (*upbound >= *lobound - epsvalue) );
#endif
}

STATIC MYBOOL presolve_altsingletonvalid(lprec *lp, int rownr, int colnr, REAL reflotest, REAL refuptest)
{
  REAL coeff_bl, coeff_bu, epsvalue = PRESOLVE_EPSVALUE;

  coeff_bl = get_rh_lower(lp, rownr);
  coeff_bu = get_rh_upper(lp, rownr);

  /* Check base data validity */
#ifdef DoPresolveRelativeTest
  if((my_reldiff(refuptest, reflotest) < -epsvalue) ||
#else  
  if((reflotest > refuptest + epsvalue) ||
#endif
     !presolve_singletonbounds(lp, rownr, colnr, &coeff_bl, &coeff_bu, NULL))
    return( FALSE );
  
  /* Base data is Ok, now check against against each other */
  if((reflotest > coeff_bu + epsvalue) ||
     (refuptest < coeff_bl - epsvalue))
    return( FALSE );
  else
    return( TRUE );
}

STATIC void presolve_initbounds(lprec *lp, int rownr, int colnr, REAL *lobound, REAL *upbound, REAL *aval,
                                REAL *UPpluvalue, REAL *UPnegvalue, REAL *LOpluvalue, REAL *LOnegvalue, int count)
{
  REAL   coeff_a, LHS, RHS, netX, Xupper, Xlower;
  
  /* Initialize to constraint bounds */
  *lobound = get_rh_lower(lp, rownr);
  *upbound = get_rh_upper(lp, rownr);
  
  if(count == 1)
    return;

  /* Get variable bounds for netting */
  LHS = *lobound;
  RHS = *upbound;
  Xlower = get_lowbo(lp, colnr);
  Xupper = get_upbo(lp, colnr);
  
  /* Identify opportunity for bound tightening */
  if(aval == NULL)
    coeff_a = get_mat(lp, rownr, colnr);
  else
    coeff_a = *aval;
  
#if 0
  if(is_chsign(lp, rownr)) {
    LHS = -LHS;
    RHS = -RHS;
    swapREAL(&LHS, &RHS);
  }
#endif
  netX = sumplumin(lp, rownr, UPpluvalue, UPnegvalue);
  if(!is_infinite(lp, LHS) && !is_infinite(lp, netX)) {
    if(coeff_a < 0) {
      LHS -= netX-coeff_a*Xupper;
      Xupper = MIN(Xupper, LHS/coeff_a);
    }
    else {
      LHS -= netX-coeff_a*Xlower;
      Xlower = MAX(Xlower, LHS/coeff_a);
    }
  }
    
  netX = sumplumin(lp, rownr, LOpluvalue, LOnegvalue);
  if(!is_infinite(lp, RHS) && !is_infinite(lp, netX)) {
    if(coeff_a < 0) {
      RHS -= netX-coeff_a*Xupper;
      Xlower = MAX(Xlower, RHS/coeff_a);
    }
    else {
      RHS -= netX-coeff_a*Xlower;
      Xupper = MIN(Xupper, RHS/coeff_a);
    }
  }
  
  *lobound = Xlower;
  *upbound = Xupper;
  
}

STATIC int presolve_tighten(lprec *lp, int i, int j, LLrec *colmap,
                            int items,
                            REAL *UPpluvalue, REAL *UPnegvalue,
                            REAL *LOpluvalue, REAL *LOnegvalue, int *count)
{
  REAL    RHlow, RHup, LObound, UPbound, Value, margin;
  int     elmnr, elmend, k, oldcount;
  MATrec  *mat = lp->matA;
  REAL    *value;
  int     *rownr;

  margin = PRESOLVE_EPSVALUE;
  
#ifdef Paranoia
  if(!isActiveLink(colmap, j))
    report(lp, SEVERE, "presolve_tighten: The selected column %d was already eliminated\n",
                       j);
#endif

  Value = get_mat(lp,i,j);
  if(Value == 0)
    return( RUNNING );

  /* Initialize and identify semicontinuous variable */
  LObound = get_lowbo(lp, j);
  UPbound = get_upbo(lp, j);
  if(is_semicont(lp, j) && (UPbound > LObound)) {
    if(LObound > 0)
      LObound = 0;
    else if(UPbound < 0)
      UPbound = 0;
  }

  /* Get singleton variable bounds */
  RHlow = get_rh_lower(lp, i);
  RHup  = get_rh_upper(lp, i);
  if(items == 1) {
    if(!presolve_singletonbounds(lp, i,j, &RHlow, &RHup, &Value))
      return( INFEASIBLE );
  }
  else {
#ifdef UseFullConstraintInfo
    presolve_initbounds(lp, i,j,  &RHlow, &RHup, &Value,
                                   UPpluvalue, UPnegvalue,
                                   LOpluvalue, LOnegvalue, count);
    if(!presolve_singletonbounds(lp, i,j, &RHlow, &RHup, &Value))
      return( INFEASIBLE );
#else
    RHlow = -lp->infinite;
    RHup = lp->infinite;
#endif
  }

  /* Look for opportunity to tighten upper and lower variable and constraint bounds */
  oldcount = (*count);
  if((RHup < lp->infinite) && (RHup+margin < UPbound)) {
    if(is_int(lp, j)) {
      if(lp->columns_scaled && is_integerscaling(lp) && (ceil(RHup) - RHup > margin))
        RHup = floor(RHup) + margin;
      else
        RHup = floor(RHup);
    }
    if(UPbound < lp->infinite) {
      elmnr = mat->col_end[j-1];
      elmend = mat->col_end[j];
      rownr = &COL_MAT_ROWNR(elmnr);
      value = &COL_MAT_VALUE(elmnr);
      for(; elmnr < elmend; 
          elmnr++, rownr += matRowColStep, value += matValueStep) {
        k = *rownr;
        Value = my_chsign(is_chsign(lp, k), *value);
        Value = unscaled_mat(lp, Value, k, j);
        if((Value > 0) && (UPpluvalue[k] < lp->infinite))
          UPpluvalue[k] += (RHup-UPbound)*Value;
        else if((Value < 0) && (UPnegvalue[i] < lp->infinite))
          UPnegvalue[k] += (RHup-UPbound)*Value;
      }
    }
    if(RHup < UPbound) {
      UPbound = RHup;
      (*count)++;
    }
  }
  if((RHlow > -lp->infinite) && (RHlow-margin > LObound)) {
    if(is_int(lp, j)) {
      if(lp->columns_scaled && is_integerscaling(lp) && (RHlow - floor(RHlow) > margin))
        RHlow = ceil(RHlow)-margin;
      else
        RHlow = ceil(RHlow);
    }
    if(LObound > -lp->infinite) {
      elmnr = mat->col_end[j-1];
      elmend = mat->col_end[j];
      rownr = &COL_MAT_ROWNR(elmnr);
      value = &COL_MAT_VALUE(elmnr);
      for(; elmnr < elmend; 
          elmnr++, rownr += matRowColStep, value += matValueStep) {
        k = *rownr;
        Value = my_chsign(is_chsign(lp, k), *value);
        Value = unscaled_mat(lp, Value, k, j);
        if((Value > 0) && (LOpluvalue[k] > -lp->infinite))
          LOpluvalue[k] += (RHlow-LObound)*Value;
        else if((Value < 0) && (LOnegvalue[k] > -lp->infinite))
          LOnegvalue[k] += (RHlow-LObound)*Value;
      }
    }
    if(RHlow > LObound) {
      LObound = RHlow;
      (*count)++;
    }
  }

  /* Now set the new variable bounds, if they are tighter */
  if((*count) > oldcount) {
#if 0     /* Experimental version */
    UPbound = presolve_round(lp, UPbound-margin, FALSE);
    LObound = presolve_round(lp, LObound+margin, TRUE);
#else     /* Safe version */
    UPbound = presolve_precision(lp, UPbound);
    LObound = presolve_precision(lp, LObound);
#endif
    if(LObound > UPbound) {
      if(LObound-UPbound < margin) {
        LObound = UPbound;
      }
      else {
        report(lp, IMPORTANT, "presolve_tighten: Found LB %g > UB %g in row %d, column %d\n",
                              LObound, UPbound, i, j);
        return(INFEASIBLE);
      }
    }
    if(lp->spx_trace || (lp->verbose > DETAILED))
      report(lp, NORMAL, "presolve_tighten: Replaced bounds on column %d to [%g ... %g]\n",
                         j, LObound, UPbound);
    set_bounds(lp, j, LObound, UPbound);
  }

  return( RUNNING );
}

STATIC void presolve_init(lprec *lp, int *plucount, int *negcount, int *pluneg)
{
  int    k, ix, ixx, colnr;
  MYBOOL isMI;
  REAL   Value, lobound, upbound;
  MATrec *mat = lp->matA;
  REAL   *value;
  int    *rownr;

  /* First do tallies; loop over nonzeros by column */
  for(colnr = 1; colnr <= lp->columns; colnr++) {

    upbound = get_upbo(lp, colnr);
    lobound = get_lowbo(lp, colnr);

    /* Handle strictly negative variables as converted to positive range with negative sign */
    isMI = is_negative(lp, colnr);

    if(is_semicont(lp, colnr) && (upbound > lobound)) {
      if(lobound > 0)
        lobound = 0;
      else if(upbound < 0)
        upbound = 0;
    }

    ix = mat->col_end[colnr - 1];
    ixx = mat->col_end[colnr];
    rownr = &COL_MAT_ROWNR(ix);
    value = &COL_MAT_VALUE(ix);
    for(; ix < ixx; 
        ix++, rownr += matRowColStep, value += matValueStep) {

     /* Retrieve row data and prepare */
      k = *rownr;
      Value = my_chsign(isMI, *value);
      Value = my_chsign(is_chsign(lp, k), Value);

     /* Cumulate counts */
      if(Value > 0)
        plucount[k]++;
      else
        negcount[k]++;
      if((lobound < 0) && (upbound > 0))
        pluneg[k]++;
    }
  }
}

STATIC int presolve_nextcol(MATrec *mat, int rownr, int prevcol, LLrec *colmap)
/* Find the first active (non-eliminated) nonzero column in rownr after prevcol */
{
  int j, jj, jb = 0, je = mat->row_end[rownr];

  if(rownr > 0)
    jb = mat->row_end[rownr-1];
  for(; jb < je; jb++) {
#if 1
    /* Narrow search window */
    jj = (jb + je) / 2;
    if(jj > jb) {
      j = ROW_MAT_COLNR(jj);
      if(j <= prevcol) {
        jb = jj;
        continue;
      }
    }
#endif
    /* Test on the left hand side of the window */
    j = ROW_MAT_COLNR(jb);
    if((j > prevcol) && isActiveLink(colmap, j))
       return( jb );
  }
  return( mat->row_end[rownr] );
}
STATIC int presolve_nextrow(MATrec *mat, int colnr, int prevrow, LLrec *rowmap)
/* Find the first active (non-eliminated) nonzero row in colnr after prevrow */
{
  int *rownr, i, ii, ib = mat->col_end[colnr-1], ie = mat->col_end[colnr];
  
  rownr = &COL_MAT_ROWNR(ib);
  for(; ib < ie; ib++, rownr += matRowColStep) {
#if 1
    /* Narrow search window */
    ii = (ib + ie) / 2;
    if(ii > ib) {
      i = COL_MAT_ROWNR(ii);
      if(i <= prevrow) {
        rownr += (ii - ib)*matRowColStep;
        ib = ii;
        continue;
      }
    }
#endif    
    /* Test on the left hand side of the window */
    if((*rownr > prevrow) && isActiveLink(rowmap, *rownr))
       return( ib );
  }
  return( mat->col_end[colnr] );
}

STATIC void presolve_rowupdate(lprec *lp, int rownr, int *collength, MYBOOL remove)
{
  if(collength != NULL) {
    if(remove) {
      int     i, ie;
      MATrec  *mat = lp->matA;

      i = mat->row_end[rownr-1];
      ie = mat->row_end[rownr];
      for(; i < ie; i++)
        collength[ROW_MAT_COLNR(i)]--;
    }
  }
}

/* Function to find if a variable can be fixed based on considering the dual */
STATIC MYBOOL presolve_coldualfix(lprec *lp, int colnr, REAL *fixValue)
{
  MYBOOL  isMI, doUP, isDualFREE = TRUE;
  int     i, ix, ie, *rownr;
  REAL    ValueOF, ValueA, *value, loX, upX;
  MATrec  *mat = lp->matA;

  /* First check basic variable range */
  loX = get_lowbo(lp, colnr);
  upX = get_upbo(lp, colnr);
  if((loX < 0) && (upX > 0))
    return( FALSE );
  isMI = (MYBOOL) (upX <= 0);

  /* Retrieve OF (standardized for minimization) */
  ix = mat->col_end[colnr - 1];
  ie = mat->col_end[colnr];
  rownr = &COL_MAT_ROWNR(ix);
  value = &COL_MAT_VALUE(ix);
  if(*rownr == 0) {
    ValueOF = *value;
    ix++;
    rownr += matRowColStep;
    value += matValueStep;
  }
  else
    ValueOF = 0;
  doUP = AUTOMATIC;
    
  /* Loop over all constraints involving active variable */
  for(; (ix < ie) && isDualFREE;
      ix++, rownr += matRowColStep, value += matValueStep) {

    /* Retrieve row data */
    i = *rownr;
    ValueA = my_chsign(is_chsign(lp, i), *value);
    
    /* Check for bound insensitivity */
    if(is_infinite(lp, get_rh_upper(lp, i))) {
      if(doUP == AUTOMATIC)
        doUP = (MYBOOL) (my_chsign(isMI, ValueA) < 0);
      isDualFREE = (MYBOOL) (my_chsign(doUP, ValueA) > 0);
    }
    else if(is_infinite(lp, get_rh_lower(lp, i))) {
      if(doUP == AUTOMATIC)
        doUP = (MYBOOL) (my_chsign(isMI, ValueA) > 0);
      isDualFREE = (MYBOOL) (my_chsign(doUP, ValueA) < 0);
    }
    else
      isDualFREE = FALSE;

  }

  /* Set fix value if we are successful (minimization objective) */
  if(isDualFREE) {
    if(ValueOF == 0)
      *fixValue = 0;
    else if(my_chsign(doUP, ValueOF) < 0) {
      doUP = (MYBOOL) (my_chsign(isMI, ValueOF) < 0);
      if(doUP) {
        if(is_infinite(lp, upX))
          isDualFREE = FALSE;
        else
          *fixValue = upX;
      }
      else {
        if(is_infinite(lp, loX))
          isDualFREE = FALSE;
        else
          *fixValue = loX;
      }
    }
    else
      isDualFREE = FALSE;
  }

  return( isDualFREE );
}

STATIC MYBOOL presolve_colupdate(lprec *lp, int colnr, int *pluneg, int *rowlength, 
                                                       int *plucount, int *negcount,
                                                       REAL *pluupper, REAL *negupper,
                                                       REAL *plulower, REAL *neglower, MYBOOL remove)
{
  int     i, ix, ie;
  MYBOOL  isneg, isMI;
  REAL    lobound, upbound, lovalue, upvalue,
          Value, fixvalue, fixprod, mult, epsvalue;
  MATrec  *mat = lp->matA;
  REAL    *value;
  int     *rownr;

  epsvalue = PRESOLVE_EPSVALUE;

  if(remove)
    mult = -1;
  else
    mult = 1;
  upbound = get_upbo(lp, colnr);
  lobound = get_lowbo(lp, colnr);

  /* Handle strictly negative variables as converted to positive range with negative sign */
  isMI = is_negative(lp, colnr);

  /* Set "exiting" value in case we are deleting a variable */
  if(upbound-lobound < epsvalue)
/*  if(my_reldiff(upbound, lobound) < epsvalue) */
    fixvalue = lp->orig_lowbo[lp->rows+colnr];
  else
    fixvalue = 0;

  /* Adjust semi-continuous variable bounds to zero-base */
  if(is_semicont(lp, colnr) && (upbound > lobound)) {
    if(lobound > 0)
      lobound = 0;
    else if(upbound < 0)
      upbound = 0;
  }

  ix = mat->col_end[colnr - 1];
  ie = mat->col_end[colnr];
  rownr = &COL_MAT_ROWNR(ix);
  value = &COL_MAT_VALUE(ix);
  for(; ix < ie; 
      ix++, rownr += matRowColStep, value += matValueStep) {

   /* Retrieve row data and adjust RHS if we are deleting a variable */
    i = *rownr;
    Value = *value;

    if(remove && (fixvalue != 0)) {
      fixprod = Value * fixvalue;
      lp->orig_rhs[i] -= fixprod;
      lp->presolve_undo->fixed_rhs[i] += fixprod;
      my_roundzero(lp->orig_rhs[i], epsvalue);
    }

   /* Prepare for further processing */
    Value = unscaled_mat(lp, Value, i, colnr);
    Value = my_chsign(is_chsign(lp, i), Value);
    isneg = (MYBOOL) (Value < 0);
    if(isMI)
      isneg = !isneg;

   /* Reduce row variable counts */
    if(remove) {
      if(rowlength != NULL)
        rowlength[i]--;
      if(isneg)
        negcount[i]--;
      else
        plucount[i]--;
      if((lobound < 0) && (upbound > 0))
        pluneg[i]--;
    }

   /* Compute associated constraint contribution values */
    if(isneg) {
      chsign_bounds(&lobound, &upbound);
      Value = -Value;
    }
    upvalue = my_if(upbound < lp->infinite, Value*upbound, lp->infinite);
    lovalue = my_if(lobound > -lp->infinite, Value*lobound, -lp->infinite);

   /* Cumulate effective upper row bound (only bother with non-finite bound) */
    if(isneg) {
      if((negupper[i] < lp->infinite) && (upvalue < lp->infinite)) {
        negupper[i] += mult*upvalue;
        negupper[i] = presolve_round(lp, negupper[i], FALSE);
      }
      else if(!remove)
        negupper[i] = lp->infinite;
    }
    else {
      if((pluupper[i] < lp->infinite) && (upvalue < lp->infinite)) {
        pluupper[i] += mult*upvalue;
        pluupper[i] = presolve_round(lp, pluupper[i], FALSE);
      }
      else if(!remove)
        pluupper[i] = lp->infinite;
    }

   /* Cumulate effective lower row bound (only bother with non-finite bound) */
    if(isneg) {
      if((neglower[i] > -lp->infinite) && (lovalue > -lp->infinite)) {
        neglower[i] += mult*lovalue;
        neglower[i] = presolve_round(lp, neglower[i], TRUE);
      }
      else if(!remove)
        neglower[i] = -lp->infinite;

      /* Remember to reset reversed bounds */
      chsign_bounds(&lobound, &upbound);
    }
    else {
      if((plulower[i] > -lp->infinite) && (lovalue > -lp->infinite)) {
        plulower[i] += mult*lovalue;
        plulower[i] = presolve_round(lp, plulower[i], TRUE);
      }
      else if(!remove)
        plulower[i] = -lp->infinite;
    }
    
   /* Validate consistency of eliminated singleton */
    if(remove && (rowlength[i] == 0)) {
      lovalue = unscaled_value(lp, sumplumin(lp, i, plulower, neglower), i);
      upvalue = unscaled_value(lp, sumplumin(lp, i, pluupper, negupper), i);
      if((upvalue < get_rh_lower(lp, i)) ||
         (lovalue > get_rh_upper(lp, i)))
        return( FALSE );
    }
  }
  return( TRUE );
}


STATIC void presolve_finalize(lprec *lp, LLrec *rowmap, LLrec *colmap)
{
  int i, ke, kb, n;

  if(colmap != NULL) {
    ke = lastInactiveLink(colmap);
    n = countInactiveLink(colmap);
    if((n > 0) && (ke > 0)) {
      kb = lastActiveLink(colmap);
      while(kb > ke)
        kb = prevActiveLink(colmap, kb);
      ke++;
      while ((n > 0) && (ke > 0)) {
        for(i = ke-1; i > kb; i--) {
#ifdef DeferredCompactColA
          del_column(lp, -i);
#else
          del_column(lp, i);
#endif
          n--;
        }
        ke = kb;
        kb = prevActiveLink(colmap, kb);
      }
    }
    freeLink(&colmap);
#ifdef DeferredCompactColA
    mat_colcompact(lp->matA, lp->presolve_undo->orig_columns);
#endif
  }
  if(rowmap != NULL) {
    ke = lastInactiveLink(rowmap);
    n = countInactiveLink(rowmap);
    if((n > 0) && (ke > 0)) {
      kb = lastActiveLink(rowmap);
      while(kb > ke)
        kb = prevActiveLink(rowmap, kb);
      ke++;
      while ((n > 0) && (ke > 0)) {
        for(i = ke-1; i > kb; i--) {
#ifdef DeferredCompactRowA
          del_constraint(lp, -i);
#else
          del_constraint(lp, i);
#endif
          n--;
        }
        ke = kb;
        kb = prevActiveLink(rowmap, kb);
      }
    }
    freeLink(&rowmap);
#ifdef DeferredCompactRowA
    mat_rowcompact(lp->matA);
#endif
  }
  mat_validate(lp->matA);
}

STATIC int presolve_collength(MATrec *mat, int column, LLrec *rowmap, int *collength)
{
  if(collength != NULL)
    column = collength[column];
  else if(rowmap == NULL)
    column = mat_collength(mat, column);
  else {
    int ib = mat->col_end[column-1],
        ie = mat->col_end[column];
    int *rownr;
    column = 0;
    rownr = &COL_MAT_ROWNR(ib);
    for(; ib < ie; ib++, rownr += matRowColStep) {
      if((*rownr == 0) ||
         isActiveLink(rowmap, *rownr))
        column++;
    }
  }
  return( column );
}
STATIC int presolve_rowlength(MATrec *mat, int row, LLrec *colmap, int *rowlength)
{
  if(rowlength != NULL)
    row = rowlength[row];
  else if(colmap == NULL)
    row = mat_rowlength(mat, row);
  else {
    int ib = 0, ie = mat->row_end[row];
    if(row > 0)
      ib = mat->row_end[row-1],
    row = 0;
    for(; ib < ie; ib++) {
      if(isActiveLink(colmap, ROW_MAT_COLNR(ib)))
        row++;
    }
  }
  return( row );
}

STATIC int presolve(lprec *lp)
{
  MYBOOL candelete;
  int    i,j,ix,iix, jx,jjx, n,nm,nn=0, nc,nv,nb,nr,ns,nt,nl;
  int    status;
  int    *plucount, *negcount, *pluneg;
  int    *rowlength = NULL, *collength = NULL;
  REAL   *pluupper, *negupper,
         *plulower, *neglower,
         Value, bound, test, epsvalue = PRESOLVE_EPSVALUE;
  MATrec *mat = lp->matA;
  LLrec  *rowmap = NULL, *colmap = NULL;

 /* Lock the variable mapping arrays and counts ahead of any row/column
    deletion or creation in the course of presolve, solvelp or postsolve */
  if(!lp->varmap_locked)  
    varmap_lock(lp);

 /* Check if we have already done presolve */
  status = RUNNING;
  mat_validate(mat);
  if(lp->wasPresolved)
    return(status);

  /* Produce original model statistics */
  REPORT_modelinfo(lp, TRUE, "Submitted:");

 /* Finalize basis indicators; if no basis was created earlier via
    set_basis or crash_basis then simply set the default basis. */
  if(!lp->basis_valid)
    lp->var_basic[0] = AUTOMATIC; /* Flag that we are presolving */

 /* Do the scaling of the problem (can also be moved to the end of the
    presolve block (before "Finish:") to possibly reduce rounding errors */
#ifndef PostScale
  auto_scale(lp);
#endif

#if 0
write_lp(lp, "test_in.lp");    /* Write to lp-formatted file for debugging */
/*write_mps(lp, "test_in.mps");*/  /* Write to lp-formatted file for debugging */
#endif

 /* Do traditional simple presolve */
  yieldformessages(lp);
  nm = 0;
  nl = 0;
  nv = 0;
  nc = 0;
  nb = 0;
  nr = 0;
  ns = 0;
  nt = 0;

  if((lp->do_presolve & PRESOLVE_LASTMASKMODE) == PRESOLVE_NONE) 
    mat_checkcounts(mat, NULL, NULL, TRUE);

  else {

    if(lp->full_solution == NULL)
      allocREAL(lp, &lp->full_solution, lp->sum_alloc+1, TRUE);

   /* Identify infeasible SOS'es prior to any pruning */
    for(i = 1; i <= SOS_count(lp); i++) {
      nn = SOS_infeasible(lp->SOS, i);
      if(nn > 0) {
        report(lp, NORMAL, "presolve: Found SOS %d (type %d) to be range-infeasible on variable %d\n",
                            i, SOS_get_type(lp->SOS, i), nn);
        status = INFEASIBLE;
        ns++;
      }
    }
    if(ns > 0)
      goto Finish;

   /* Create row and column counts */
#if 1
    allocINT(lp,  &rowlength, lp->rows + 1, TRUE);
    for(i = 1; i <= lp->rows; i++)
      rowlength[i] = presolve_rowlength(mat, i, NULL, NULL);
    allocINT(lp,  &collength, lp->columns + 1, TRUE);
    for(j = 1; j <= lp->columns; j++)
      collength[j] = presolve_collength(mat, j, NULL, NULL);
#endif
      
   /* Create NZ count and sign arrays, and do general initialization of row bounds */
    allocINT(lp,  &plucount,  lp->rows + 1, TRUE);
    allocINT(lp,  &negcount,  lp->rows + 1, TRUE);
    allocINT(lp,  &pluneg,    lp->rows + 1, TRUE);
    allocREAL(lp, &pluupper,  lp->rows + 1, TRUE);
    allocREAL(lp, &negupper,  lp->rows + 1, TRUE);
    allocREAL(lp, &plulower,  lp->rows + 1, TRUE);
    allocREAL(lp, &neglower,  lp->rows + 1, TRUE);
    createLink(lp->rows, &rowmap, NULL);
      fillLink(rowmap);
    createLink(lp->columns, &colmap, NULL);
      fillLink(colmap);
    presolve_init(lp, plucount, negcount, pluneg);

   /* Accumulate constraint bounds based on bounds on individual variables */
Restart:
    nm++;
    for(j = firstActiveLink(colmap); j != 0; j = nextActiveLink(colmap, j)) {
      presolve_colupdate(lp, j, pluneg, NULL, 
                                plucount, negcount,
                                pluupper, negupper,
                                plulower, neglower, FALSE);
      if((status == RUNNING) && userabort(lp, -1))
        goto Complete;
    }

    /* Reentry point */
Redo:
    nl++;
    nn = 0;

   /* Eliminate empty or fixed columns (including trivial OF column singletons) */
    if(is_presolve(lp, PRESOLVE_COLS) && mat_validate(mat)) {
      if(userabort(lp, -1))
        goto Complete;
      n = 0;
      for(j = lastActiveLink(colmap); j > 0; ) {
        candelete = FALSE;
        ix = lp->rows + j;
        iix = presolve_collength(mat, j, rowmap, collength);
        if(iix == 0) {
          if(SOS_is_member(lp->SOS, 0, j))
            report(lp, NORMAL, "presolve: Found empty variable %d as member of a SOS\n",
                                get_col_name(lp,j));
          else {
            if(lp->orig_lowbo[ix] != 0)
              report(lp, DETAILED, "presolve: Found empty non-zero variable %s\n",
                                    get_col_name(lp,j));
            candelete = TRUE;
          }
        }
        else if(isOrigFixed(lp, ix)) {
          if(SOS_is_member(lp->SOS, 0, j))
            continue;
          report(lp, DETAILED, "presolve: Eliminated variable %s fixed at %g\n",
                                get_col_name(lp,j), get_lowbo(lp, j));
          candelete = TRUE;
        }
        else if((iix == 1) && (COL_MAT_ROWNR(mat->col_end[j-1]) == 0) &&
                !SOS_is_member(lp->SOS, 0, j)) {

          if(get_OF_raw(lp, ix) > 0)
            test = get_lowbo(lp, j);
          else
            test = get_upbo(lp, j);
          if(is_infinite(lp, test)) {
            report(lp, DETAILED, "presolve: Unbounded variable %s\n",
                                  get_col_name(lp,j));
            status = UNBOUNDED;
          }
          else {
            /* Fix the value at its best bound */
            set_bounds(lp, j, test, test);
            report(lp, DETAILED, "presolve: Eliminated trivial variable %s fixed at %g\n",
                                  get_col_name(lp,j), test);
            candelete = TRUE;
          }
        }
#ifdef ApplyDualPresolve  /* Look for opportunity to fix column based on the dual */
        else if(presolve_coldualfix(lp, j, &Value)) {
          if(is_infinite(lp, Value)) {
            report(lp, DETAILED, "presolve: Unbounded variable %s\n",
                                  get_col_name(lp,j));
            status = UNBOUNDED;
          }
          else {
            /* Fix the value at its best bound */
            set_bounds(lp, j, Value, Value);
            report(lp, DETAILED, "presolve: Eliminated dual-zero variable %s fixed at %g\n",
                                get_col_name(lp,j), Value);
            candelete = TRUE;
          }
        }
#endif

        ix = j;
        j = prevActiveLink(colmap, j);
        if(candelete) {
          Value = get_lowbo(lp, ix);
          lp->full_solution[lp->presolve_undo->orig_rows + 
                            lp->presolve_undo->var_to_orig[lp->rows + ix]] = Value;
          if(!presolve_colupdate(lp, ix, pluneg, rowlength,
                                         plucount, negcount,
                                         pluupper, negupper,
                                         plulower, neglower, TRUE)) {
            report(lp, NORMAL, "presolve: Found variable bound infeasibility for column %d\n", ix);
            status = INFEASIBLE;
            nn = 0;
            break;
          }
          removeLink(colmap, ix);
          n++;
          nv++;
          nn++;
        }
      }
    }

   /* Eliminate linearly dependent rows; loop backwards over every row */
    if(is_presolve(lp, PRESOLVE_LINDEP) && mat_validate(mat)) {
      int firstix, RT1, RT2;
      if(userabort(lp, -1))
        goto Complete;
      n = 0;
      for(i = lastActiveLink(rowmap); (i > 0) && (status == RUNNING); ) {

        /* First scan for rows with identical row lengths */
        ix = prevActiveLink(rowmap, i);
        if(ix == 0)
          break;

        /* Don't bother about empty rows or row singletons since they are
           handled by PRESOLVE_ROWS */
        j = presolve_rowlength(mat, i, colmap, rowlength);
        if(j <= 1) {
          i = ix;
          continue;
        }

#if 0
        /* Enable this to scan all rows back */
        RT2 = lp->rows;
#else
        RT2 = 1;
#endif
        firstix = ix;
        for(RT1 = 0; (ix > 0) && (RT1 < RT2) && (status == RUNNING);
            ix = prevActiveLink(rowmap, ix), RT1++)  {
          candelete = FALSE;
          if(presolve_rowlength(mat, ix, colmap, rowlength) != j)
            continue;

          /* Check if the beginning columns are identical; if not, continue */
          iix = presolve_nextcol(mat, ix, 0, colmap);
          jjx = presolve_nextcol(mat, i,  0, colmap);

          if(ROW_MAT_COLNR(iix) != ROW_MAT_COLNR(jjx))
            continue;

          /* We have a candidate row; check if the entries have a fixed non-zero ratio */
          test  = get_mat_byindex(lp, iix, TRUE, FALSE);
          Value = get_mat_byindex(lp, jjx, TRUE, FALSE);
          bound = test / Value;
          Value = bound;

          /* Loop over remaining entries */
          jx = mat->row_end[i];
          jjx = presolve_nextcol(mat, i, ROW_MAT_COLNR(jjx), colmap);
          for(; (jjx < jx) && (Value == bound);
              jjx = presolve_nextcol(mat, i, ROW_MAT_COLNR(jjx), colmap)) {
             iix = presolve_nextcol(mat, ix, ROW_MAT_COLNR(iix), colmap);
             if(ROW_MAT_COLNR(iix) != ROW_MAT_COLNR(jjx))
               break;
             test  = get_mat_byindex(lp, iix, TRUE, FALSE);
             Value = get_mat_byindex(lp, jjx, TRUE, FALSE);

             /* If the ratio is different from the reference value we have a mismatch */
             Value = test / Value;
             if(bound == lp->infinite)
               bound = Value;
             else if(fabs(Value - bound) > epsvalue)
               break;
          }

          /* Check if we found a match (we traversed all active columns without a break) */
          if(jjx >= jx) {

            /* Get main reference values */
            test  = lp->orig_rhs[ix];
            Value = lp->orig_rhs[i] * bound;

            /* First check for inconsistent equalities */
            if((fabs(test - Value) > epsvalue) &&
               ((get_constr_type(lp, ix) == EQ) && (get_constr_type(lp, i) == EQ))) {
              status = INFEASIBLE;
            }

            else {

              /* Update lower and upper bounds */
              if(is_chsign(lp, i) != is_chsign(lp, ix))
                bound = -bound;

              test = get_rh_lower(lp, i);
              if(test <= -lp->infinite)
                test *= my_sign(bound);
              else
                test *= bound;

              Value = get_rh_upper(lp, i);
              if(Value >= lp->infinite)
                Value *= my_sign(bound);
              else
                Value *= bound;

              if((bound < 0))
                swapREAL(&test, &Value);

              if(get_rh_lower(lp, ix) < test)
                set_rh_lower(lp, ix, test);
              if(get_rh_upper(lp, ix) > Value)
                set_rh_upper(lp, ix, Value);

              /* Check results */
              test  = get_rh_lower(lp, ix);
              Value = get_rh_upper(lp, ix);
              if(fabs(Value-test) < epsvalue)
                set_constr_type(lp, ix, EQ);
              else if(Value < test) {
                status = INFEASIBLE;
              }  

              /* Verify if we can continue */
              candelete = (MYBOOL) (status == RUNNING);
              if(!candelete) {
                report(lp, IMPORTANT, "presolve: Range infeasibility found involving rows %d and %d\n",
                                      ix, i);
              }
            }
          }
          /* Perform i-row deletion if authorized */
          if(candelete) {
            presolve_rowupdate(lp, i, collength, TRUE);
            removeLink(rowmap, i);
            n++;
            nc++;
          }
        }
        i = firstix;
      }
    }

   /* Aggregate and tighten bounds using 2-element EQs */
    if(FALSE &&
       (lp->equalities > 0) && is_presolve(lp, PRESOLVE_AGGREGATE) && mat_validate(mat)) {
      if(userabort(lp, -1))
        goto Complete;
      n = 0;
      for(i = lastActiveLink(rowmap); (i > 0) && (status == RUNNING); ) {
        /* Find an equality constraint with 2 elements; the pivot row */
        if(!is_constr_type(lp, i, EQ) || (presolve_rowlength(mat, i, colmap, rowlength) != 2)) {
          i = prevActiveLink(rowmap, i);
          continue;
        }
        /* Get the column indeces of NZ-values of the pivot row */
        jx = mat->row_end[i-1];
        j =  mat->row_end[i];
        for(; jx < j; jx++)
          if(isActiveLink(colmap, ROW_MAT_COLNR(jx)))
            break;
        jjx = jx+1;
        for(; jjx < j; jjx++)
          if(isActiveLink(colmap, ROW_MAT_COLNR(jjx)))
            break;
        jx  = ROW_MAT_COLNR(jx);
        jjx = ROW_MAT_COLNR(jjx);
        if(SOS_is_member(lp->SOS, 0, jx) && SOS_is_member(lp->SOS, 0, jjx)) {
          i = prevActiveLink(rowmap, i);
          continue;
        }
        /* Determine which column we should eliminate (index in jx) :
           1) the longest column
           2) the variable not being a SOS member
           3) an integer variable  */
        if(presolve_collength(mat, jx, rowmap, collength) < 
           presolve_collength(mat, jjx, rowmap, collength))
          swapINT(&jx, &jjx);
        if(SOS_is_member(lp->SOS, 0, jx))
          swapINT(&jx, &jjx);
        if(!is_int(lp, jx) && is_int(lp, jjx))
          swapINT(&jx, &jjx);
        /* Whatever the priority above, we must have bounds to work with;
           give priority to the variable with the smallest bound */
        test  = get_upbo(lp, jjx)-get_lowbo(lp, jjx);
        Value = get_upbo(lp, jx)-get_lowbo(lp, jx);
        if(test < Value)
          swapINT(&jx, &jjx);
        /* Try to set tighter bounds on the non-eliminated variable (jjx) */
        test  = get_mat(lp, i, jjx); /* Non-eliminated variable coefficient a */
        Value = get_mat(lp, i, jx);  /* Eliminated variable coefficient     b */
#if 1
        bound = get_lowbo(lp, jx);
        if((bound > -lp->infinite)) {
          bound = (get_rh(lp, i)-Value*bound) / test;
          if(bound < get_upbo(lp, jjx)-epsvalue)
            set_upbo(lp, jjx, presolve_round(lp, bound, FALSE));
        }
        bound = get_upbo(lp, jx);
        if((bound < lp->infinite)) {
          bound = (get_rh(lp, i)-Value*bound) / test;
          if(bound > get_lowbo(lp, jjx)+epsvalue)
            set_lowbo(lp, jjx, presolve_round(lp, bound, TRUE));
        }
        i = prevActiveLink(rowmap, i);
#else
        /* Loop over the non-zero rows of the column to be eliminated;
           substitute jx-variable by updating rhs and jjx coefficients */
        int iiix;
        ix = mat->col_end[jx-1];
        iiix = mat->col_end[jx];
        rownr = &COL_MAT_ROWNR(ix);
        value = &COL_MAT_VALUE(ix);
        for(; ix < iiix; 
            ix++, rownr += matRowColStep, value += matValueStep) {
          REAL newvalue;
          iix = *rownr;
          if((iix == i) ||
             ((iix > 0) && !isActiveLink(rowmap, iix)))
            continue;
          /* Do the update */
          bound = unscaled_mat(lp, *value, iix, jx)/Value;
          bound = my_chsign(is_chsign(lp, iix), bound);
          newvalue = get_mat(lp, iix, jjx) - bound*test;
            set_mat(lp, iix, jjx, presolve_precision(lp, newvalue));
          newvalue = get_rh(lp, iix) - bound*get_rh(lp, i);
            set_rh(lp, iix, presolve_precision(lp, newvalue));
        }
        /* Delete the column */
        removeLink(colmap, jx);
        nc++;
        n++;
        /* Delete the row */
        ix = i;
        i = prevActiveLink(rowmap, i);
        presolve_rowupdate(lp, ix, collength, TRUE);
        removeLink(rowmap, ix);
        nr++;
        n++;
        mat_validate(mat);
#endif
      }
    }

#if 0
   /* Increase A matrix sparsity by discovering common subsets using 2-element EQs */
    if((lp->equalities > 0) && is_presolve(lp, PRESOLVE_SPARSER) && mat_validate(mat)) {
      int iiix;
      if(userabort(lp, -1))
        goto Complete;
      n = 0;
      for(i = lastActiveLink(rowmap); (i > 0) && (status == RUNNING); ) {
        candelete = FALSE;
        /* Find an equality constraint with 2 elements; the pivot row */
        if(!is_constr_type(lp, i, EQ) || (mat_rowlength(mat, i) != 2)) {
          i = prevActiveLink(rowmap, i);
          continue;
        }
        /* Get the column indeces of NZ-values of the pivot row */
        jx = mat->row_end[i-1];
        jx = ROW_MAT_COLNR(jx);
        jjx = mat->row_end[i];
        jjx = ROW_MAT_COLNR(jjx);
        /* Scan to find a row with matching column entries */
        ix = lp->col_end[jx-1];
        iiix = lp->col_end[jx];
        rownr = &COL_MAT_ROWNR(ix);
        for(; ix < iiix; 
            ix++, rownr += matRowColStep) {
          if(*rownr == i)
            continue;
          /* We now have a single matching value, find the next */
          iix = lp->col_end[jjx-1];
          for(; iix < lp->col_end[jjx]; iix++)
            if(COL_MAT_ROWNR(iix) >= ix)
              break;
          /* Abort this row if there was no second column match */
          if((iix >= lp->col_end[jjx]) || (COL_MAT_ROWNR(iix) > ix) )
            break;
          /* Otherwise, do variable subsitution and mark pivot row for deletion */
          candelete = TRUE;
          nc++;
          /*
           ... Add remaining logic later!
          */
        }
        ix = i;
        i = prevActiveLink(rowmap, i);
        if(candelete) {
          presolve_rowupdate(lp, ix, collength, TRUE);
          removeLink(rowmap, ix);
          n++;
        }
      }
    }
#endif

   /* Eliminate empty rows, convert row singletons to bounds,
      tighten bounds, and remove always satisfied rows */
    if(is_presolve(lp, PRESOLVE_ROWS) && mat_validate(mat)) {
      if(userabort(lp, -1))
        goto Complete;
      n = 0;
      for(i = lastActiveLink(rowmap); i > 0; ) {
      
        candelete = FALSE;

       /* First identify any full row infeasibilities */
#if 0       
        test = epsvalue;
#else
        test = DEF_EPSSOLUTION;
#endif
        if(!is_constr_type(lp, i, LE)) {
          Value = MAX(sumplumin(lp, i, plulower,neglower),
                      sumplumin(lp, i, pluupper,negupper));
          if(Value < get_rh(lp, i)-test) {
            report(lp, NORMAL, "presolve: Found upper bound infeasibility in row %d\n", i);
            status = INFEASIBLE;
            break;
          }
        }
        if(!is_constr_type(lp, i, GE)) {
          Value = MIN(sumplumin(lp, i, plulower,neglower),
                      sumplumin(lp, i, pluupper,negupper));
          if(Value > get_rh(lp, i)+test) {
            report(lp, NORMAL, "presolve: Found lower bound infeasibility in row %d\n", i);
            status = INFEASIBLE;
            break;
          }
        }  

        j = plucount[i]+negcount[i];

       /* Delete non-zero rows and variables that are completely determined;
          note that this step can provoke infeasibility in some tight models */
        if((j > 0)                                       /* Only examine non-empty rows, */
           && (fabs(lp->orig_rhs[i]) < epsvalue)         /* .. and the current RHS is zero, */
           && ((plucount[i] == 0) || (negcount[i] == 0)) /* .. and the parameter signs are all equal, */
           && (pluneg[i] == 0)                           /* .. and no (quasi) free variables, */
           && (is_constr_type(lp, i, EQ)
#define PresolveFindImpliedEqualities
#ifdef PresolveFindImpliedEqualities
               || (fabs(get_rh_lower(lp, i)-sumplumin(lp, i, pluupper,negupper)) < epsvalue)  /* Convert to equalities */
               || (fabs(get_rh_upper(lp, i)-sumplumin(lp, i, plulower,neglower)) < epsvalue)  /* Convert to equalities */
#endif
              )
              ) {

          /* Delete the columns of this row, but make sure we don't delete SOS variables */
          for(ix = mat->row_end[i]-1; ix >= mat->row_end[i-1]; ix--) {
            jx = ROW_MAT_COLNR(ix);
            if(isActiveLink(colmap, jx) && !SOS_is_member(lp->SOS, 0, jx)) {
              if(!presolve_colupdate(lp, jx, pluneg, rowlength,
                                             plucount, negcount,
                                             pluupper, negupper,
                                             plulower, neglower, TRUE)) {
                report(lp, NORMAL, "presolve: Found variable bound infeasibility for column %d\n", jx);
                status = INFEASIBLE;                             
                nn = 0;
                break;
              }
              removeLink(colmap, jx);
              nv++;
            }
          }
          /* Then delete the row, which is redundant */
          if(status == RUNNING) {
            candelete = TRUE;
            nc++;
          }  
        }
        else

       /* Then delete any empty or always satisfied / redundant row that cannot at
          the same time guarantee that we can also delete associated variables */
        if((j == 0) ||                                   /* Always delete an empty row */
           ((j > 1) &&
            (pluneg[i] == 0) && ((plucount[i] == 0) ||
                                 (negcount[i] == 0)) &&  /* Consider removing if block above is ON! */
            (sumplumin(lp, i, pluupper,negupper)-sumplumin(lp, i, plulower,neglower) < epsvalue))      /* .. or if it is always satisfied (redundant) */
          ) {
          candelete = TRUE;
          nc++;
        }

       /* Convert row singletons to bounds (delete fixed columns in columns section) */
        else if((j == 1) &&
                (sumplumin(lp, i, pluupper,negupper)-
                 sumplumin(lp, i, plulower,neglower) >= epsvalue)) {
          j = presolve_nextcol(mat, i, 0, colmap);
          j = ROW_MAT_COLNR(j);

          /* Make sure we do not have conflicting other singleton rows with this variable */
          Value = lp->infinite;
          test = -Value;
          if(presolve_collength(mat, j, rowmap, collength) > 1) {
            test  = get_rh_lower(lp, i);
            Value = get_rh_upper(lp, i);
            if(presolve_singletonbounds(lp, i, j, &test, &Value, NULL)) {
              jx = mat->col_end[j];
              for(ix = presolve_nextrow(mat, j, 0, rowmap);
                  ix < jx; ix = presolve_nextrow(mat, j, iix, rowmap)) {
                iix = COL_MAT_ROWNR(ix);
                if((iix != i) && 
                   (presolve_rowlength(mat, iix, colmap, rowlength) == 1) &&
                   !presolve_altsingletonvalid(lp, iix, j, test, Value)) {
                  status = INFEASIBLE;
                  nn = 0;
                  break;
                }
              }
            }
            else {
              status = INFEASIBLE;
            }
          }

          /* Proceed to remove variable */
          if(status == RUNNING) {
            if((fabs(test-Value) < epsvalue) && (fabs(test) > epsvalue) &&
               !SOS_is_member(lp->SOS, 0, j)) {
              set_bounds(lp, j, test, test);
              if(!presolve_colupdate(lp, j, pluneg, rowlength,
                                            plucount, negcount,
                                            pluupper, negupper,
                                            plulower, neglower, TRUE)) {
                status = INFEASIBLE; 
              }
              removeLink(colmap, j);
              nv++;
            }
            else
              status = presolve_tighten(lp, i, j, colmap, plucount[i]+negcount[i],
                                                          pluupper, negupper,
                                                          plulower, neglower, &nt);
          }
          if(status == INFEASIBLE) {
            report(lp, NORMAL, "presolve: Found variable bound infeasibility for column %d\n", j);
            nn = 0;
            break;
          }
          candelete = TRUE;
          nb++;
        }

       /* Check if we have a constraint made redundant through bounds on individual variables */
        else if((sumplumin(lp, i, plulower,neglower) >= get_rh_lower(lp, i)-epsvalue) &&
                (sumplumin(lp, i, pluupper,negupper) <= get_rh_upper(lp, i)+epsvalue)) {
          candelete = TRUE;
          nc++;
        }

#ifdef AggressiveRowPresolve
       /* Look for opportunity to tighten constraint bounds;
          know to create problems with scaled ADLittle.mps */
        else if(j > 1) {
          test = sumplumin(lp, i, plulower,neglower);
          if(test > get_rh_lower(lp, i)+epsvalue) {
            set_rh_lower(lp, i, presolve_round(lp, test, TRUE));
            nr++;
          }
          test = sumplumin(lp, i, pluupper,negupper);
          if(test < get_rh_upper(lp, i)-epsvalue) {
            set_rh_upper(lp, i, presolve_round(lp, test, FALSE));
            nr++;
          }
        }
#endif

        /* Get next row and do the deletion of the previous, if indicated */
        ix = i;
        i = prevActiveLink(rowmap, i);
        if(candelete) {
          presolve_rowupdate(lp, ix, collength, TRUE);
          removeLink(rowmap, ix);
          n++;
          nn++;
        }
        /* Look for opportunity to convert ranged constraint to equality-type */
        else if(!is_constr_type(lp, ix, EQ) && (get_rh_range(lp, ix) < epsvalue))
          set_constr_type(lp, ix, EQ);
      }
    }
   
   /* Try again if we were successful in this presolve loop */
    if((status == RUNNING) && !userabort(lp, -1)) {
      if(nn > 0) goto Redo;

      /* Optionally do an extra loop from scratch */
      if((nm < 1) && (nc+nb+nt+nv+nr > 0)) {
        MEMCLEAR(plucount, lp->rows+1);
        MEMCLEAR(negcount, lp->rows+1);
        MEMCLEAR(pluneg  , lp->rows+1);
        MEMCLEAR(pluupper, lp->rows+1);
        MEMCLEAR(negupper, lp->rows+1);
        MEMCLEAR(plulower, lp->rows+1);
        MEMCLEAR(neglower, lp->rows+1);
        goto Restart;
      }
    }

Complete:
   /* See if we can convert some constraints to SOSes (only SOS1 handled) */
    if(is_presolve(lp, PRESOLVE_SOS) &&
       (MIP_count(lp) > 0) && mat_validate(mat)) {
      n = 0;
      for(i = lastActiveLink(rowmap); i > 0; ) {
        candelete = FALSE;
        test = get_rh(lp, i);
        jx = get_constr_type(lp, i);
#ifdef EnableBranchingOnGUB
        if((test == 1) && (jx != GE)) {
#else
        if((test == 1) && (jx == LE)) {
#endif
          jjx = mat->row_end[i-1];
          iix = mat->row_end[i];
          for(; jjx < iix; jjx++) {
            j = ROW_MAT_COLNR(jjx);
            if(!isActiveLink(colmap, j))
              continue;
            if(!is_binary(lp, j) || (get_mat(lp, i, j) != 1))
              break;
          }
          if(jjx >= iix) {
            char SOSname[16];

            /* Define a new SOS instance */
            sprintf(SOSname, "SOS_%d", SOS_count(lp) + 1);
            ix = add_SOS(lp, SOSname, 1, 1, 0, NULL, NULL);
            if(jx == EQ)
              SOS_set_GUB(lp->SOS, ix, TRUE);
            Value = 0;
            jjx = mat->row_end[i-1];
            for(; jjx < iix; jjx++) {
              j = ROW_MAT_COLNR(jjx);
              if(!isActiveLink(colmap, j))
                continue;
              Value += 1;
              append_SOSrec(lp->SOS->sos_list[ix-1], 1, &j, &Value);
            }
            candelete = TRUE;
            nc++;
          }
        }

        /* Get next row and do the deletion of the previous, if indicated */
        ix = i;
        i = prevActiveLink(rowmap, i);
        if(candelete) {
          presolve_rowupdate(lp, ix, collength, TRUE);
          removeLink(rowmap, ix);
          n++;
          nn++;
        }
      }
      if(n)
        report(lp, NORMAL, "presolve: Converted %5d constraints to SOS1.\n", n);
    }

   /* Finalize presolve */
#ifdef Paranoia
    i = presolve_validate(lp, rowmap, colmap);
    if(i > 0)
      report(lp, SEVERE, "presolve: %d internal consistency failure(s) detected\n", i);
#endif
    presolve_finalize(lp, rowmap, colmap);

   /* Tighten MIP bound if possible (should ideally use some kind of smart heuristic) */
#ifndef PostScale
    if((MIP_count(lp) > 0) || (get_Lrows(lp) > 0)) {
      if(is_maxim(lp))
        lp->bb_heuristicOF = MAX(lp->bb_heuristicOF, sumplumin(lp, 0, plulower, neglower));
      else
        lp->bb_heuristicOF = MIN(lp->bb_heuristicOF, sumplumin(lp, 0, pluupper, negupper));
    }
#endif

   /* Report summary information */
    if(nv)
      report(lp, NORMAL, "presolve: Removed   %6d empty or fixed variables.\n", nv);
    if(nb)
      report(lp, NORMAL, "presolve: Converted %6d row singletons to variable bounds.\n", nb);
    if(nt)
      report(lp, NORMAL, "presolve: Tightened %6d other variable bounds.\n", nt);
    if(nc)
      report(lp, NORMAL, "presolve: Removed   %6d empty or redundant constraints.\n", nc);
    if(nr)
      report(lp, NORMAL, "presolve: Tightened %6d constraint bounds.\n", nr);

    if(nc+nb+nt+nv+nr > 0)
      report(lp, NORMAL, " \n");

    /* Report optimality or infeasibility */
    if(sumplumin(lp, 0, pluupper,negupper)-sumplumin(lp, 0, plulower,neglower) < epsvalue) {
      report(lp, NORMAL, "presolve: Identified optimal OF value of %g\n",
                         sumplumin(lp, 0, pluupper,negupper));
#ifndef PostScale
      lp->bb_limitOF = sumplumin(lp, 0, pluupper,negupper);
#endif
#if 0
      status = OPTIMAL;
#endif
    }
    else if(status != RUNNING)
      report(lp, NORMAL, "presolve: Infeasibility or unboundedness detected.\n");

    /* Clean up */
    FREE(rowlength);
    FREE(collength);
    FREE(plucount);
    FREE(negcount);
    FREE(pluneg);
    FREE(plulower);
    FREE(neglower);
    FREE(pluupper);
    FREE(negupper);

  }

  /* Signal that we are done presolving */
  if((lp->usermessage != NULL) && 
     ((lp->do_presolve & PRESOLVE_LASTMASKMODE) != 0) && (lp->msgmask & MSG_PRESOLVE))
     lp->usermessage(lp, lp->msghandle, MSG_PRESOLVE);
     
  /* Clean out empty SOS records */
  if(SOS_count(lp) > 0) {
    clean_SOSgroup(lp->SOS);
    if(lp->SOS->sos_count == 0)
      free_SOSgroup(&(lp->SOS));
  }
  
  /* Create master SOS variable list */
  if(SOS_count(lp) > 0)
    make_SOSchain(lp, (MYBOOL) ((lp->do_presolve & PRESOLVE_LASTMASKMODE) != PRESOLVE_NONE));

  /* Resolve GUBs */
#ifdef EnableBranchingOnGUB
  if(is_bb_mode(lp, NODE_GUBMODE))
    identify_GUB(lp, TRUE);
#endif

  /* Crash the basis, if specified */
  crash_basis(lp);

#ifdef PostScale
  if(status == RUNNING)
    auto_scale(lp);
#endif

  /* Produce presolved model statistics */
  if(nc+nb+nt+nv+nr > 0)
    REPORT_modelinfo(lp, FALSE, "Presolved:");

Finish:
  lp->wasPresolved = TRUE;

#if 0
  write_lp(lp, "test_out.lp");   /* Must put here due to variable name mapping */
#endif
#if 0
  REPORT_debugdump(lp, "testint2.txt", FALSE);
#endif

  return( status );

}

STATIC MYBOOL postsolve(lprec *lp, int status)
{
  /* Verify solution */
  if(lp->lag_status != RUNNING) {
    int itemp;
    
    if((status == OPTIMAL) || (status == SUBOPTIMAL)) {
      itemp = check_solution(lp, lp->columns, lp->best_solution, 
                                 lp->orig_upbo, lp->orig_lowbo, DEF_EPSSOLUTION);
      if((itemp != OPTIMAL) && (lp->spx_status == OPTIMAL))
          lp->spx_status = itemp;
      else if((itemp == OPTIMAL) && (status == SUBOPTIMAL))
        lp->spx_status = status;
    }
    else {
      report(lp, NORMAL, "lp_solve unsuccessful after %d iterations and a last best value of %g\n",
             lp->total_iter, lp->best_solution[0]);
      if(lp->bb_totalnodes > 0)
        report(lp, NORMAL, "lp_solve explored %d nodes before termination\n",
               lp->bb_totalnodes);
    }
  }
  
  /* Check if we can clear the variable map */
  if(varmap_canunlock(lp))
    lp->varmap_locked = FALSE;
    
  return( TRUE );
}
