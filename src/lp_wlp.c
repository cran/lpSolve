
#include "lp_lib.h"
#include "lp_utils.h"
#include "lp_report.h"
#include "lp_wlp.h"

#ifdef FORTIFY
# include "lp_fortify.h"
#endif


/* ------------------------------------------------------------------------- */
/* Input and output of lp format model files for lp_solve                    */
/* ------------------------------------------------------------------------- */

STATIC void write_lpcomment(FILE *output, char *string, MYBOOL newlinebefore)
{
  fprintf(output, "%s/* %s */\n", (newlinebefore) ? "\n" : "", string);
}

STATIC MYBOOL write_lprow(lprec *lp, int rowno, FILE *output)
{
  int     i, ie, j;
  REAL    a;
  MATrec  *mat = lp->matA;
  MYBOOL  first = TRUE, rowwritten;

  if(rowno == 0)
    i = 0;
  else
    i = mat->row_end[rowno-1];
  ie = mat->row_end[rowno];
  rowwritten = (MYBOOL)(i < ie);
  for(; i < ie; i++) {
    j = ROW_MAT_COLNR(i);
    if(is_splitvar(lp, j))
      continue;
#if 0
    a = get_mat_byindex(lp, i, TRUE);
#else
    a = ROW_MAT_VALUE(i);
    a = my_chsign(is_chsign(lp, rowno), a);
    a = unscaled_mat(lp, a, rowno, j);
#endif
    if(!first)
      fputc(' ', output);
    else
      first = FALSE;
    if(a == -1)
      fprintf(output, "-");
    else if(a == 1)
      fprintf(output, "+");
    else
      fprintf(output, "%+.12g ", (double)a);
    fprintf(output, "%s", get_col_name(lp, j));
  }
  return(rowwritten);
}

MYBOOL LP_writefile(lprec *lp, char *filename)
{
  int    i, j, b;
  MYBOOL ok;
  REAL   a;
  FILE   *output = stdout;
  char   *ptr;

#ifdef Paranoia
  if(lp->matA->is_roworder) {
    report(lp, IMPORTANT, "LP_writefile: Cannot write to LP file while in row entry mode.\n");
    return(FALSE);
  }
#endif
  if(!mat_validate(lp->matA)) {
    report(lp, IMPORTANT, "LP_writefile: Could not validate the data matrix.\n");
    return(FALSE);
  }

  ok = (MYBOOL) ((filename == NULL) || ((output = fopen(filename,"w")) != NULL));
  if(!ok)
    return(ok);
  if(filename == NULL && lp->outstream != NULL)
    output = lp->outstream;

  /* Write name of model */
  ptr = get_lp_name(lp);
  if(ptr != NULL) {
    if(*ptr)
      write_lpcomment(output, ptr, FALSE);
    else
      ptr = NULL;
  }

  /* Write the objective function */
  write_lpcomment(output, "Objective function", (MYBOOL) (ptr != NULL));
  if(is_maxim(lp))
    fprintf(output, "max: ");
  else
    fprintf(output, "min: ");

  write_lprow(lp, 0, output);
  a = get_rh(lp, 0);
  if(a)
    fprintf(output, " %+.12g", a);
  fprintf(output, ";\n");

  /* Write constraints */
  if(lp->rows > 0)
    write_lpcomment(output, "Constraints", TRUE);
  for(j = 1; j <= lp->rows; j++) {
    if(lp->names_used && (lp->row_name[j] != NULL))
      ptr = get_row_name(lp, j);
    else
      ptr = NULL;
    if((ptr != NULL) && (*ptr))
      fprintf(output, "%s: ", ptr);

#ifndef SingleBoundedRowInLP
    /* Write the ranged part of the constraint, if specified */
    if ((lp->orig_upbo[j]) && (lp->orig_upbo[j] < lp->infinite)) {
      if(my_chsign(is_chsign(lp, j), lp->orig_rhs[j]) == -lp->infinite)
        fprintf(output, "-Inf %s ", (is_chsign(lp, j)) ? ">=" : "<=");
      else if(my_chsign(is_chsign(lp, j), lp->orig_rhs[j]) == lp->infinite)
        fprintf(output, "+Inf %s ", (is_chsign(lp, j)) ? ">=" : "<=");
      else
        fprintf(output, "%+.12g %s ",
                (lp->orig_upbo[j]-lp->orig_rhs[j]) * (is_chsign(lp, j) ? 1.0 : -1.0) / (lp->scaling_used ? lp->scalars[j] : 1.0),
                (is_chsign(lp, j)) ? ">=" : "<=");
    }
#endif

    if((!write_lprow(lp, j, output)) && (get_Ncolumns(lp) >= 1))
      fprintf(output, "0 %s", get_col_name(lp, 1));

    if(lp->orig_upbo[j] == 0)
      fprintf(output, " =");
    else if(is_chsign(lp, j))
      fprintf(output, " >=");
    else
      fprintf(output, " <=");
    if(fabs(get_rh(lp, j) + lp->infinite) < 1)
      fprintf(output, " -Inf;\n");
    else if(fabs(get_rh(lp, j) - lp->infinite) < 1)
      fprintf(output, " +Inf;\n");
    else
      fprintf(output, " %.12g;\n", get_rh(lp, j));

#ifdef SingleBoundedRowInLP
    /* Write the ranged part of the constraint, if specified */
    if ((lp->orig_upbo[j]) && (lp->orig_upbo[j] < lp->infinite)) {
      if(lp->names_used && (lp->row_name[j] != NULL))
        ptr = get_row_name(lp, j);
      else
        ptr = NULL;
      if((ptr != NULL) && (*ptr))
        fprintf(output, "%s: ", ptr);
      if((!write_lprow(lp, j, output)) && (get_Ncolumns(lp) >= 1))
        fprintf(output, "0 %s", get_col_name(lp, 1));
      fprintf(output, " %s %g;\n",
                     (is_chsign(lp, j)) ? "<=" : ">=",
                     (lp->orig_upbo[j]-lp->orig_rhs[j]) * (is_chsign(lp, j) ? 1.0 : -1.0) / (lp->scaling_used ? lp->scalars[j] : 1.0));
    }
#endif
  }

  /* Write bounds on variables */
  ok = FALSE;
  for(i = lp->rows + 1; i <= lp->sum; i++)
    if(!is_splitvar(lp, i - lp->rows)) {
      if(lp->orig_lowbo[i] == lp->orig_upbo[i]) {
        if(!ok) {
	  write_lpcomment(output, "Variable bounds", TRUE);
	  ok = TRUE;
	}
        fprintf(output, "%s = %.12g;\n", get_col_name(lp, i - lp->rows), get_upbo(lp, i - lp->rows));
      }
      else {
#ifndef SingleBoundedRowInLP
        if((lp->orig_lowbo[i] != 0) && (lp->orig_upbo[i] < lp->infinite)) {
          if(!ok) {
	    write_lpcomment(output, "Variable bounds", TRUE);
	    ok = TRUE;
	  }
          if(lp->orig_lowbo[i] == -lp->infinite)
            fprintf(output, "-Inf");
          else
            fprintf(output, "%.12g", get_lowbo(lp, i - lp->rows));
          fprintf(output, " <= %s <= ", get_col_name(lp, i - lp->rows));
          if(lp->orig_lowbo[i] == lp->infinite)
            fprintf(output, "+Inf");
          else
            fprintf(output, "%.12g", get_upbo(lp, i - lp->rows));
          fprintf(output, ";\n");
	}
        else
#endif
        {
          if(lp->orig_lowbo[i] != 0) {
            if(!ok) {
	      write_lpcomment(output, "Variable bounds", TRUE);
	      ok = TRUE;
	    }
      	    if(lp->orig_lowbo[i] == -lp->infinite)
	      fprintf(output, "%s >= -Inf;\n", get_col_name(lp, i - lp->rows));
      	    else if(lp->orig_lowbo[i] == lp->infinite)
	      fprintf(output, "%s >= +Inf;\n", get_col_name(lp, i - lp->rows));
	    else
              fprintf(output, "%s >= %.12g;\n",
                              get_col_name(lp, i - lp->rows), get_lowbo(lp, i - lp->rows));
	  }
	  if(lp->orig_upbo[i] != lp->infinite) {
            if(!ok) {
	      write_lpcomment(output, "Variable bounds", TRUE);
	      ok = TRUE;
	    }
            fprintf(output, "%s <= %.12g;\n",
                            get_col_name(lp, i - lp->rows), get_upbo(lp, i - lp->rows));
	  }
        }
      }
    }

  /* Write optional integer section */
  if(lp->int_count > 0) {
    write_lpcomment(output, "Integer definitions", TRUE);
    i = 1;
    while(i <= lp->columns && !is_int(lp, i))
      i++;
    if(i <= lp->columns) {
      fprintf(output, "int %s", get_col_name(lp, i));
      i++;
      for(; i <= lp->columns; i++)
        if((!is_splitvar(lp, i)) && (is_int(lp, i)))
          fprintf(output, ",%s", get_col_name(lp, i));
      fprintf(output, ";\n");
    }
  }

  /* Write optional SEC section */
  if(lp->sc_count > 0) {
    write_lpcomment(output, "Semi-continuous variables", TRUE);
    i = 1;
    while(i <= lp->columns && !is_semicont(lp, i))
      i++;
    if(i <= lp->columns) {
      fprintf(output, "sec %s", get_col_name(lp, i));
      i++;
      for(; i <= lp->columns; i++)
        if((!is_splitvar(lp, i)) && (is_semicont(lp, i)))
		    fprintf(output, ",%s", get_col_name(lp, i));
      fprintf(output, ";\n");
    }
  }

  /* Write optional SOS section */
  if(SOS_count(lp) > 0) {
    SOSgroup *SOS = lp->SOS;
    write_lpcomment(output, "SOS definitions", TRUE);
    for(b = 0, i = 0; i < SOS->sos_count; b = SOS->sos_list[i]->priority, i++) {
      fprintf(output, "SOS\n%s: ",
              (SOS->sos_list[i]->name == NULL) ||
              (*SOS->sos_list[i]->name==0) ? "SOS" : SOS->sos_list[i]->name); /* formatnumber12((double) lp->sos_list[i]->priority) */

      for(a = 0.0, j = 1; j <= SOS->sos_list[i]->size; a = SOS->sos_list[i]->weights[j], j++)
        if(SOS->sos_list[i]->weights[j] == ++a)
          fprintf(output, "%s%s",
                  (j > 1) ? "," : "",
                  get_col_name(lp, SOS->sos_list[i]->members[j]));
        else
          fprintf(output, "%s%s:%.12g",
                  (j > 1) ? "," : "",
                  get_col_name(lp, SOS->sos_list[i]->members[j]),
		  SOS->sos_list[i]->weights[j]);
      if(SOS->sos_list[i]->priority == ++b)
        fprintf(output, " <= %d;\n", SOS->sos_list[i]->type);
      else
        fprintf(output, " <= %d:%d;\n", SOS->sos_list[i]->type, SOS->sos_list[i]->priority);
    }
  }

  ok = TRUE;

  if(filename != NULL)
    fclose(output);
  return(ok);
}

MYBOOL LP_writehandle(lprec *lp, FILE *output)
{
  set_outputstream(lp, output);
  return(LP_writefile(lp, NULL));
}
