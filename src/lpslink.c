/*
** This is callable from S-Plus.
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include "lpkit.h"

#ifdef BUILDING_FOR_R
#define LONG_OR_INT int
#else
#define LONG_OR_INT long
#endif

char *status_file_name = "c:/lpstatus.txt";
char *status_file_2 = "c:/lpstatus2.txt";
FILE *status_file = (FILE *) NULL;

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
              LONG_OR_INT *status)            /* Holds return value          */
{
int i, result;
double *const_ptr;

lprec *lp;

/** status_file = fopen (status_file_name, "w"); **/

lp = make_lp ((int) 0, *x_count);

if (lp == (lprec *) NULL) {
/**
    fprintf (status_file, "Unable to create lp!\n");
    fclose (status_file);
**/
    *status = -1;
    return;
}

result = set_obj_fn (lp, objective);
if (result == 0) {
/**
    fprintf (status_file, "Unable to set objective!\n");
    fclose (status_file);
**/
    *status = -2;
    return;
}

if (*direction == 1)
    set_maxim (lp);

if ((int) *const_count > 0) {
    const_ptr = constraints;
    for (i = 0; i < (int) *const_count; i++)
    {
        add_constraint (lp, const_ptr,
            (short) const_ptr[(int) (*x_count) + 1],
                    const_ptr[(int) (*x_count) + 2]);
        const_ptr += (int) *x_count + 2;
    }
}

if (*int_count > 0) {
    for (i = 0; i < (int) *int_count; i++)
        set_int (lp, (int) (int_vec[i]), TRUE);
}

/**
fprintf (status_file, "As best I can tell, all is well.\n");
**/
/*
**
** Uncomment this to print out the lp to yet another file.
**
** print_file ("c:/lpzap.txt");
** print_lp (lp);
** print_file ("stdout");
*/

*status = (LONG_OR_INT) solve (lp);

if ((int) *status != 0) {
/**
    fprintf (status_file, "Unable to solve, error %i\n", (int) *status);
    fclose (status_file);
**/
    return;
}

*obj_val = get_objective (lp);
/** fprintf (status_file, "Solved, obj is %f\n", get_objective (lp)); **/

get_variables (lp, solution);


delete_lp (lp);


/**
fprintf (status_file, "Deleted!\n");
fclose (status_file);
**/

return;

} /* end "lpslink" */
