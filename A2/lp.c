/****************************************************************************/
/*                                                                          */
/*              CO759: LP interface for HW2, Subtour TSP Solver             */
/*              Date:  January 27, 2015                                     */
/*                     February 9, 2015 (added CO759lp_setbnd ())           */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <cplex.h>
#include "lp.h"

int CO759lp_init (CO759lp *lp)
{
    int rval = 0;

    lp->cplex_lp = (CPXLPptr) NULL;
    lp->cplex_env = CPXopenCPLEX (&rval);
    if (rval) {
        fprintf (stderr, "CPXoptnCPLEX failed, return code %d\n", rval);
        goto CLEANUP;
    }

CLEANUP:
    return rval;
}

void CO759lp_free (CO759lp *lp)
{
    if (lp) {
        if (lp->cplex_env) {
            if (lp->cplex_lp) {
                CPXfreeprob (lp->cplex_env, &(lp->cplex_lp));
            } 
            CPXcloseCPLEX (&lp->cplex_env);
            lp->cplex_env = (CPXENVptr) NULL;
        }
    }
}

int CO759lp_create (CO759lp *lp, const char *name)
{
    int rval = 0;
    char nambuf[32];

    if (!lp->cplex_env) {
        fprintf (stderr, "no cplex_env in lp\n");
        rval = 1; goto CLEANUP;
    }

    /* We copy the name since CPXcreateprob doesn't declare name as const */
    strncpy (nambuf, name, sizeof (nambuf));
    nambuf[sizeof(nambuf)-1] = '\0';

    lp->cplex_lp = CPXcreateprob (lp->cplex_env, &rval, nambuf);
    if (!lp->cplex_lp || rval) {
       fprintf (stderr, "CPXcreateprob failed, return code %d\n", rval);
       goto CLEANUP;
    }

CLEANUP:
    return rval;
}

int CO759lp_new_row (CO759lp *lp, char sense, double rhs)
{
    int rval = 0;
    char asense[1];
    double arhs[1];

    asense[0] = sense;
    arhs[0] = rhs;

    rval = CPXnewrows (lp->cplex_env, lp->cplex_lp, 1, arhs, asense,
                       (double *) NULL, (char **) NULL);
    if (rval) { fprintf (stderr, "CPXnewrows failed\n"); goto CLEANUP; }

CLEANUP:
    return rval;
}

int CO759lp_addrows (CO759lp *lp, int newrows, int newnz, double *rhs,
        char *sense, int *rmatbeg, int *rmatind, double *rmatval)
{
    int rval = 0;

    rval = CPXaddrows (lp->cplex_env, lp->cplex_lp, 0, newrows, newnz, rhs,
             sense, rmatbeg, rmatind, rmatval, (char **) NULL, (char **) NULL);
    if (rval) { fprintf (stderr, "CPXaddrows failed\n"); goto CLEANUP; }

CLEANUP:
    return rval;
}

int CO759lp_addcols (CO759lp *lp, int newcols, int newnz, double *obj,
        int *cmatbeg, int *cmatind, double *cmatval, double *lb, double *ub)
{
    int rval = 0;

    rval = CPXaddcols (lp->cplex_env, lp->cplex_lp, newcols, newnz, obj,
                   cmatbeg, cmatind, cmatval, lb, ub, (char **) NULL);
    if (rval) { fprintf (stderr, "CPXaddcols failed\n"); goto CLEANUP; }

CLEANUP:
    return rval;
}

int CO759lp_setbnd (CO759lp *lp, int col, char lower_or_upper, double bnd)
{
    int rval = 0;
    int cindex[1];
    double bd[1];
    char lu[1];

    cindex[0] = col;
    lu[0] = lower_or_upper;
    bd[0] = bnd;

    rval = CPXchgbds (lp->cplex_env, lp->cplex_lp, 1, cindex, lu, bd);
    if (rval) {
        fprintf (stderr, "Couldn't set bnd on variable %d in cplex\n", col);
        goto CLEANUP;
    }

CLEANUP:
    return rval;
}


int CO759lp_opt (CO759lp *lp, int *infeasible)
{
    int rval = 0, solstat;

    *infeasible = 0;
 
    rval = CPXdualopt (lp->cplex_env, lp->cplex_lp);
    if (rval) {
        fprintf (stderr, "CPXdualopt failed\n"); goto CLEANUP;
    }

    solstat =  CPXgetstat (lp->cplex_env, lp->cplex_lp);
    if (solstat == CPX_STAT_INFEASIBLE) {
        if (infeasible) *infeasible = 1;
    } else if (solstat != CPX_STAT_OPTIMAL &&
               solstat != CPX_STAT_OPTIMAL_INFEAS) {
        fprintf (stderr, "Cplex optimization status %d\n", solstat);
        rval = 1; goto CLEANUP;
    }

CLEANUP:
    return rval;
}

int CO759lp_objval (CO759lp *lp, double *obj)
{
    int rval = 0;

    rval = CPXgetobjval (lp->cplex_env, lp->cplex_lp, obj);
    if (rval) {
        fprintf (stderr, "CPXgetobjval failed\n"); goto CLEANUP;
    }

CLEANUP:
    return rval;
}

int CO759lp_x (CO759lp *lp, double *x)
{
    int rval = 0, ncols;

    ncols = CPXgetnumcols (lp->cplex_env, lp->cplex_lp);
    if (ncols == 0) {
        fprintf (stderr, "No columns in LP\n");
        rval = 1; goto CLEANUP;
    }
    rval = CPXgetx (lp->cplex_env, lp->cplex_lp, x, 0, ncols - 1);
    if (rval) {
        fprintf (stderr, "CPXgetx failed\n"); goto CLEANUP;
    }

CLEANUP:
    return rval;
}


int CO759lp_write (CO759lp *lp, const char *fname)
{
    int rval = 0;
    char nambuf[32], lpbuf[4];

    strncpy (nambuf, fname, sizeof (nambuf));
    nambuf[sizeof(nambuf)-1] = '\0';
    strncpy (lpbuf, "RLP", sizeof ("RLP"));

    rval = CPXwriteprob (lp->cplex_env, lp->cplex_lp, nambuf, lpbuf);
    if (rval) {
        fprintf (stderr, "CPXwriteprob failed\n"); goto CLEANUP;
    }

CLEANUP:
    return rval;
}


