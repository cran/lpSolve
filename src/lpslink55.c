/*
** lpslink.c: function to link lpsolve with Excel, S-Plus or R.
*/
/*
** In addition to standard "include" files we need lpkit.h, supplied by
** lpsolve. This gives definitions for (for example) the "lprec" structure.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "lp_lib.h"

/*
** In R "integers" get passed as longs, whereas in S-Plus they're ints.
** So the user should turn on this flag to compile for R.
*/

#ifdef BUILDING_FOR_R
#define LONG_OR_INT int
#else
#define LONG_OR_INT long
#endif

/*
** The declaration of the link function.
*/

void lpslink (LONG_OR_INT *direction,         /* 1 for max, 0 for min        */
              LONG_OR_INT *x_count,           /* Number of x's               */
              double *objective,              /* Objective function          */
              LONG_OR_INT *const_count,       /* Number of constraints       */
              double *constraints,            /* Has extra element on front  */
              LONG_OR_INT *int_count,         /* Number of integer variables */
              LONG_OR_INT *int_vec,           /* Indices of int. variables   */
              double *obj_val,                /* Objective function value    */
              double *solution,               /* Result of call              */
              LONG_OR_INT *presolve,          /* Value of presolve           */
              LONG_OR_INT *compute_sens,      /* Want sensitivity?           */
              double *sens_coef_from,         /* Sens. coef. lower limit     */
              double *sens_coef_to,           /* Sens. coef. upper limit     */
              double *duals,                  /* Dual values                 */
              double *duals_from,             /* Lower limit dual values     */
              double *duals_to,               /* Lower limit dual values     */
              LONG_OR_INT *status);           /* Holds return value          */

/*
** Some globals for calling from VBScript. The caller will call lpslink,
*/
static long    vb_direction;
static long    vb_x_count;
static double *vb_objective;
static long    vb_const_count;
static double *vb_constraints;
static long    vb_int_count;
static long   *vb_int_vec;
static double  vb_obj_val;
static double *vb_solution;

/*
**************************** lps_vb_setup *****************************
**
**
** lps_vb_setup: set up an lp problem (that is, allocate necessary space)
*/

long lps_vb_setup (long direction,    /* Optimization direction (1 = max) */
                   long x_count,      /* Number of variables              */
                   long const_count,  /* Number of constraints            */
                   long int_count)    /* Number of integer variables      */
{
long  i; /* iteration variable */

/*
** Set globals (which start "vb") from parameters (which do not).
*/

vb_direction   = direction;
vb_x_count     = x_count;
vb_const_count = const_count;
vb_int_count   = int_count;

/*
** Allocate objective (which holds the coefficients in the objective function).
** We need an extra element at the front. If malloc() failes, get out.
*/
vb_objective = (double *) malloc (1 + sizeof (double) * vb_x_count);
if (vb_objective == (double *) NULL)
    return (-1);

vb_objective[0] = 0.0;

/*
** Allocate constraints. This, too, gets an extra entry. If the allocation
** fails, exit gracefully.
*/

vb_constraints = (double *) malloc (sizeof (double) * 
                     (1 + vb_const_count * (vb_x_count + 2)));

if (vb_constraints == (double *) NULL)
{
    free (vb_objective);
    return (-1);
}

vb_constraints[0] = 0.0;

/*
** If any variables are constrained to be integers, allocate one long
** for each such variable. Quit if the allocation fails. If it succeeds,
** put zeros in there. We will insert the proper numbers later.
*/

if (vb_int_count > 0) {
    vb_int_vec = (long *) malloc (1 + sizeof (long) * vb_int_count);
    if (vb_int_vec == (long *) NULL)
    {
        free (vb_objective);
        free (vb_constraints);
        return (-1);
    }

    for (i = 0L; i <= vb_int_count; i++) /* there's one extra */
        vb_int_vec[i] = 0L;
}

/*
** Allocate space for the solution. This will hold the optimal values for each
** coefficient. If the allocation fails, quit.
*/

vb_solution = (double *) malloc (sizeof (double) * vb_x_count);
if (vb_solution == (double *) NULL)
{
    free (vb_objective);
    free (vb_constraints);
    if (vb_int_count > 0)
        free (vb_int_vec);
    return (-1);
}
/*
** Our work here is done.
*/ 
return (0);
} /* end lps_vb_setup */

/***************************** lps_vb_set_element ********************/

long lps_vb_set_element (long type,     /* Place to do the setting */
                         long row,      /* Row in which to set     */
                         long column,   /* Column in which to set  */
                         double value)  /* Value to set            */
{
/*
** This function allows us to set an element of the lp. If "type" = 1,
** we ignore column and set the "row-th" element of vb_objective to "value."
** If type is 2, set the row, column-th element of constraints (allowing for
** the funny layout) to value; if type = 3, set the rowth variable to integer.
*/

    if (type == 1)
        vb_objective[row] = value;
    if (type == 2) {
        vb_constraints[(row - 1) * (vb_x_count + 2) + column] = value;
    }
    if (type == 3 && vb_int_count > 0) 
        vb_int_vec[row] = floor (value + 0.5);
    return (1);
}

/***************************** lps_vb_get_element ********************/
double lps_vb_get_element (long type,     /* Place to get from           */
                           long row,      /* Row to get from             */
                           long column)   /* Column to get from (unused) */
{
/*
** Get an element. If type is 1, get the objective value; if type is 2,
** get the rowth element of the solution. "Column" is currently unused.
*/
    if (type == 1)
        return (vb_obj_val);
    if (type == 2)
        return (vb_solution[row]);
    return (0.0);
}

/*
********************************* lps_vb_cleanup *************************
**
** lps_vb_cleanup: free all the things allocated in lps_vb_setup.
*/
long lps_vb_cleanup (long unused)
{
    free (vb_objective);
    free (vb_constraints);
    free (vb_int_vec);
    free (vb_solution);
    return (0);
}

/******************************** lpslink ************************************/

void lpslink (LONG_OR_INT *direction,         /* 1 for max, 0 for min        */
              LONG_OR_INT *x_count,           /* Number of x's               */
              double *objective,              /* Objective function          */
              LONG_OR_INT *const_count,       /* Number of constraints       */
              double *constraints,            /* Has extra element on front  */
              LONG_OR_INT *int_count,         /* Number of integer variables */
              LONG_OR_INT *int_vec,           /* Indices of int. variables   */
              double *obj_val,                /* Objective function value    */
              double *solution,               /* Result of call              */
              LONG_OR_INT *presolve,          /* Value of presolve           */
              LONG_OR_INT *compute_sens,      /* Want sensitivity?           */
              double *sens_coef_from,         /* Sens. coef. lower limit     */
              double *sens_coef_to,           /* Sens. coef. upper limit     */
              double *duals,                  /* Dual values                 */
              double *duals_from,             /* Lower limit dual values     */
              double *duals_to,               /* Lower limit dual values     */
              LONG_OR_INT *status)            /* Holds return value          */
{
/*
** This is the function called from the outside.
*/
int i,               /* Iteration variable      */
    result;          /* Holds result of calls   */

double *const_ptr;   /* Points to a constraint   */

lprec *lp;           /* Structure to hold the lp */

/*
** Make an empty lp with x_count variables. If it fails, return.
*/
lp = make_lp ((int) 0, *x_count);

if (lp == (lprec *) NULL)
    return;

set_verbose (lp, 1); /* CRITICAL */

/*
** "Objective" is a vector. Set the objective function. Return on fail.
*/
result = set_obj_fn (lp, objective);
if (result == 0)
    return;

/* Set the direction. The default is minimize, but set it anyway. */
if (*direction == 1)
    set_maxim (lp);
else
    set_minim (lp);

/*
** If there are any constraints, point "constr_ptr" at the first one.
*/
if ((int) *const_count > 0) {
    const_ptr = constraints;
/*
** Add constraints, one at a time; then move constr_ptr up.
*/

    for (i = 0; i < (int) *const_count; i++)
    {
        add_constraint (lp, const_ptr,
            (short) const_ptr[(int) (*x_count) + 1], 
                    const_ptr[(int) (*x_count) + 2]);
        const_ptr += (int) *x_count + 2;
    }
} /* end "if there are any constraints. */

if (*int_count > 0) {
    for (i = 0; i < (int) *int_count; i++)
        set_int (lp, (int) (int_vec[i]), TRUE);
}

/*
** Presolve if needed (that is, if we are going to want sensitivity in
** an integer program) then solve the lp. If "status" is non-zero, return.
*/

if (*compute_sens > 0) {
    if (*int_count > 0)
        set_presolve (lp, PRESOLVE_SENSDUALS, get_presolveloops (lp));
    }

*status = (LONG_OR_INT) solve (lp);

if ((int) *status != 0) {
    delete_lp (lp);
    return;
}

/* Now get the sensitivities, if requested. */
if (*compute_sens > 0) {
    get_sensitivity_obj (lp, sens_coef_from, sens_coef_to);
    get_sensitivity_rhs (lp, duals, duals_from, duals_to);
}

/*
** We've succeeded. Extract the objective function's value and
** the values of the variables.
*/

*obj_val = get_objective (lp);

get_variables (lp, solution);

/*
** 
*/

/*
** Free up the memory and return.
*/

delete_lp (lp);

return;

} /* end "lpslink" */

/*
****************************** lp_transbig ************************
**
** This function handles "big" transportation problem. It takes in
** the number of rows and columns, the cost matrix, and the signs and
** right-hand sides of the constraints and builds everything else on the
** fly.
*/
void lp_transbig (LONG_OR_INT *direction,     /* 1 for max, 0 for min       */
              LONG_OR_INT *r_count,           /* Number of rows             */
              LONG_OR_INT *c_count,           /* Number of columns          */
              double *costs,                  /* Objective function         */
              LONG_OR_INT *r_signs,           /* Signs of row constraints   */
              double *r_rhs,                  /* RHS of row constraints     */
              LONG_OR_INT *c_signs,           /* Signs of col constraints   */
              double *c_rhs,                  /* RHS of col constraints     */
              double *obj_val,                /* Objective function value   */
              LONG_OR_INT *int_count,         /* How many vars are integers?*/
              LONG_OR_INT *integers,          /* Which vars. are integer?   */
              double *solution,               /* Result of call             */
              LONG_OR_INT *presolve,          /* Value of presolve          */
              LONG_OR_INT *compute_sens,      /* Want sensitivity?          */
              double *sens_coef_from,         /* Sens. coef. lower limit    */
              double *sens_coef_to,           /* Sens. coef. upper limit    */
              double *duals,                  /* Dual values                */
              double *duals_from,             /* Lower limit dual values    */
              double *duals_to,               /* Lower limit dual values    */
              LONG_OR_INT *status)            /* Holds return value         */
{
long i;              /* Iteration variable       */
long result;         /* Holds result of calls    */
long this_element;   /* Which are we looking at? */
lprec *lp;           /* Structure to hold the lp */
double *row_vals;    /* Holds the values for row-type constraints */
int *col_inds;       /* Holds locations for col-type constraints  */
double *col_vals;    /* Holds the values for col-type constraints */
int *row_inds;       /* Holds locations for row-type constraints  */

long col_ind_ctr, row_ind_ctr;
long rc = *r_count, cc = *c_count, num_vars = *r_count * *c_count;

/*
** Make an empty lp with r_count x c_count variables. If it fails, return.
*/
lp = make_lp ((int) 0, *r_count * *c_count);

if (lp == (lprec *) NULL)
    return;

set_verbose (lp, 1); /* CRITICAL */

set_add_rowmode (lp, TRUE);
/*
** "Costs" is already a vector. Set the objective function. Return on fail.
*/
result = set_obj_fn (lp, costs);
if (result == 0)
    return;

/* Set the direction. The default is minimize, but set it anyway. */
if (*direction == 1)
    set_maxim (lp);
else
    set_minim (lp);

/*
** Add constraints. There are r_count row-type constraints, plus c_count
** col_type constraints.
*/
row_vals = calloc (cc, sizeof (double));
col_inds = calloc (cc, sizeof (int));

for (row_ind_ctr = 0L; row_ind_ctr < rc; row_ind_ctr++)
{
    for (col_ind_ctr = 0; col_ind_ctr < cc; col_ind_ctr++) {
        row_vals[col_ind_ctr] = 1.0;
        this_element = 1 + (col_ind_ctr * rc) + row_ind_ctr;
        col_inds[col_ind_ctr] = this_element;
    }
    add_constraintex (lp, cc, row_vals, col_inds, r_signs[row_ind_ctr], r_rhs[row_ind_ctr]);
}

free (row_vals);
free (col_inds);

col_vals = calloc (rc, sizeof (double));
row_inds = calloc (rc, sizeof (int));

for (col_ind_ctr = 0L; col_ind_ctr < cc; col_ind_ctr++)
{
    for (row_ind_ctr = 0; row_ind_ctr < rc; row_ind_ctr++) {
        col_vals[row_ind_ctr] = 1.0;
        this_element = 1 + row_ind_ctr + col_ind_ctr * rc;
        row_inds[row_ind_ctr] = this_element;
    }
    add_constraintex (lp, rc, col_vals, row_inds, c_signs[col_ind_ctr], c_rhs[col_ind_ctr]);
}
free (col_vals);
free (row_inds);

set_add_rowmode (lp, FALSE);

/*
** Set integers. set_int starts counting at 1.
*/
if (*int_count  > 0)
    for (i = 1; i <= *int_count; i++)
        set_int (lp, integers[i], 1); /* Variable in ith element of integers */

if (*compute_sens > 0) {
    set_presolve (lp, PRESOLVE_SENSDUALS, 10);
}

*status = (LONG_OR_INT) solve (lp);

if ((int) *status != 0) {
    return;
}

/* Now get the sensitivities, if requested. */
if (*compute_sens > 0) {
    get_sensitivity_obj (lp, sens_coef_from, sens_coef_to);
    get_sensitivity_rhs (lp, duals, duals_from, duals_to);
}

/*
** We've succeeded. Extract the objective function's value and
** the values of the variables.
*/

*obj_val = get_objective (lp);

get_variables (lp, solution);

/*
** 
*/

/*
** Free up the memory and return.
*/

delete_lp (lp);
}
