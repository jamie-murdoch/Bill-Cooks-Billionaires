/****************************************************************************/
/*                                                                          */
/*              CO759: LP interface for HW2, Subtour TSP Solver             */
/*              Date:  January 27, 2015                                     */
/*                                                                          */
/****************************************************************************/

#include <cplex.h>

typedef struct CO759lp {
    CPXENVptr cplex_env;
    CPXLPptr  cplex_lp;
} CO759lp;

int CO759lp_init (CO759lp *lp);
int CO759lp_create (CO759lp *lp, const char *name);
int CO759lp_write (CO759lp *lp, const char *fname);
void CO759lp_free (CO759lp *lp);

int CO759lp_new_row (CO759lp *lp, char sense, double rhs);
int CO759lp_addrows (CO759lp *lp, int newrows, int newnz, double *rhs,
    char *sense, int *rmatbeg, int *rmatind, double *rmatval);
int CO759lp_addcols (CO759lp *lp, int newcols, int newnz, double *obj,
    int *cmatbeg, int *cmatind, double *cmatval, double *lb, double *ub);
int CO759lp_setbnd (CO759lp *lp, int col, char lower_or_upper, double bnd);

int CO759lp_opt (CO759lp *lp, int *infeasible);
int CO759lp_objval (CO759lp *lp, double *obj);
int CO759lp_x (CO759lp *lp, double *x);
