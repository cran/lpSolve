#ifndef HEADER_lp_presolve
#define HEADER_lp_presolve

#include "lp_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Put function headers here */

STATIC MYBOOL presolve_createUndo(lprec *lp);
STATIC MYBOOL inc_presolve_space(lprec *lp, int delta, MYBOOL isrows);
STATIC MYBOOL presolve_setOrig(lprec *lp, int orig_rows, int orig_cols);
STATIC MYBOOL presolve_fillUndo(lprec *lp, int orig_rows, int orig_cols, MYBOOL setOrig);
STATIC MYBOOL presolve_freeUndo(lprec *lp);

STATIC int presolve(lprec *lp);
STATIC MYBOOL postsolve(lprec *lp, int status);

#ifdef __cplusplus
 }
#endif

#endif /* HEADER_lp_presolve */
