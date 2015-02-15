/****************************************************************************/
/*                                                                          */
/*              CO759: Model for HW2, TSP via Subtour Cuts                  */
/*              Date:  January 27, 2015                                     */
/*                     February 9, 2015   added add_all_subtours()          */
/*                     February 11, 2015  added add_connect()               */
/*                                                                          */
/****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <cplex.h>
#include "lp.h"
#include "util.h"

typedef struct ajdobj {
    int n;   /* index of neighbor node */
    int e;   /* index of adj joining neighbor */
} adjobj;

typedef struct node {
    int deg;
    adjobj *adj;
    int mark;
} node;

typedef struct graph {
    int ncount;
    int ecount;
    node *nodelist;
    adjobj *adjspace;
} graph;

static void usage (char *f);
static int getprob (char *fname, int *p_ncount, int *p_ecount, int **p_elist,
    int **p_elen);
static int parseargs (int ac, char **av);
static int subtour (int ncount, int ecount, int *elist, int *elen, int *tlist);
static int euclid_edgelen (int i, int j, double *x, double *y);
static int add_all_subtours (int ncount, int ecount, int *elist, CO759lp *lp);
static void next_set (int sz, int *Set);
static void get_delta (int nsize, int *nlist, int ecount, int *elist,
   int *deltacount, int *delta, int *marks);
static int add_subtour (CO759lp *lp, int deltacount, int *delta);
static int add_connect (int ncount, int ecount, int *elist, CO759lp *lp);
static void init_graph (graph *G);
static void free_graph (graph *G);
static int build_graph (int ncount, int ecount, int *elist, graph *G);
static int connected (graph *G, double *x, int *icount, int *island);
static void dfs (int n, graph *G, double *x, int *icount, int *island);

static char *fname = (char *) NULL;
static int seed = 0;
static int geometric_data = 0;
static int ncount_rand = 0;
static int gridsize_rand = 100;
static int use_all_subtours = 0;

int main (int ac, char **av)
{
    int rval  = 0, ncount = 0, ecount = 0;
    int *elist = (int *) NULL, *elen = (int *) NULL, *tlist = (int *) NULL;
    double szeit;

    seed = (int) CO759_real_zeit ();

    rval = parseargs (ac, av);
    if (rval) goto CLEANUP;

    if (!fname && !ncount_rand) {
        printf ("Must specify a problem file or use -k for random prob\n");
        rval = 1; goto CLEANUP;
    }
    printf ("Seed = %d\n", seed);
    srandom (seed);

    if (fname) {
        printf ("Problem name: %s\n", fname);
        if (geometric_data) printf ("Geometric data\n");
    }

    rval = getprob (fname, &ncount, &ecount, &elist, &elen);
    if (rval) {
        fprintf (stderr, "getprob failed\n"); goto CLEANUP;
    }

    if (use_all_subtours && ncount > 20) {
        fprintf (stderr, "Too many nodes to add all subtours\n"); goto CLEANUP;
    }

    tlist = (int *) malloc ((ncount)*sizeof (int));
    if (!tlist) {
        fprintf (stderr, "out of memory for tlist\n");
        rval = 1;  goto CLEANUP; 
    }

    szeit = CO759_zeit ();
    rval = subtour (ncount, ecount, elist, elen, tlist);
    if (rval) {
        fprintf (stderr, "subtour failed\n");
        goto CLEANUP;
    }
    printf ("Running Time: %.2f seconds\n", CO759_zeit() - szeit);
    fflush (stdout);

CLEANUP:

    if (tlist) free (tlist);
    if (elist) free (elist);
    if (elen) free (elen);
    return rval;
}

#define LP_EPSILON 0.000001

static int subtour (int ncount, int ecount, int *elist, int *elen, int *tlist)
{
    int rval = 0, i, j, infeasible = 0;
    double  obj[1], lb[1], ub[1], objval, *x = (double *) NULL;
    int     cmatbeg[1], cmatind[2];
    double  cmatval[2];
    CO759lp lp;

    rval = CO759lp_init (&lp);
    if (rval) { fprintf (stderr, "CO759lp_init failed\n"); goto CLEANUP; }

    rval = CO759lp_create (&lp, "subtour");
    if (rval) { fprintf (stderr, "CO759lp_create failed\n"); goto CLEANUP; }

    /* Build a row for each degree equation */

    for (i = 0; i < ncount; i++) {
        rval = CO759lp_new_row (&lp, 'E', 2.0);
        if (rval) {
            fprintf (stderr, "CO759lp_new_row failed\n"); goto CLEANUP;
        }
    }

    /* Build a column for each edge of the graph */

    cmatbeg[0] = 0;
    cmatval[0] = 1.0;
    cmatval[1] = 1.0;
    for (j = 0; j < ecount; j++) {
        obj[0]     = (double) elen[j];
        lb[0]      = 0.0;
        ub[0]      = 1.0;
        cmatind[0] = elist[2*j];
        cmatind[1] = elist[2*j+1];
        rval = CO759lp_addcols (&lp, 1 /* # of new variables */,
           2 /* # of new nonzeros */, obj, cmatbeg, cmatind, cmatval, lb, ub);
        if (rval) {
            fprintf (stderr, "CClp_addcols failed\n"); goto CLEANUP;
        }
    }

    rval = CO759lp_write (&lp, "subtour.lp");
    if (rval) {
        fprintf (stderr, "CO759lp_write failed\n"); goto CLEANUP;
    }

    rval = CO759lp_opt (&lp, &infeasible);
    if (rval) {
        fprintf (stderr, "CO759lp_opt failed\n"); goto CLEANUP;
    }
    if (infeasible) {
        fprintf (stderr, "LP is infeasible\n"); 
        rval = 1; goto CLEANUP;
    }

    rval = CO759lp_objval (&lp, &objval);
    if (rval) {
        fprintf (stderr, "CO759lp_objval failed\n"); goto CLEANUP;
    }

    printf ("Degree-Equation LP Value: %f\n", objval);
    fflush (stdout);

    x = (double *) malloc (ecount * sizeof (double));
    if (!x) {
        fprintf (stderr, "out of memory for x\n");
        rval = 1; goto CLEANUP;
    }

    rval = CO759lp_x (&lp, x);
    if (rval) {
        fprintf (stderr, "CO759lp_x failed\n"); goto CLEANUP;
    }
 
    for (i = 0, j = 0; j < ecount; j++) {
        if (x[j] > LP_EPSILON) i++;
    }

    printf ("LP graph has %d edges\n", i);
    for (j = 0; j < ecount; j++) {
        if (x[j] > LP_EPSILON) {
            printf ("%d %d %f\n", elist[2*j], elist[2*j+1], x[j]);
        }
    }
    fflush (stdout);

    if (use_all_subtours) {
        rval = add_all_subtours (ncount, ecount, elist, &lp);
        if (rval) {
            fprintf (stderr, "add_all_subtours failed\n"); goto CLEANUP;
        }
    } else {
        rval = add_connect (ncount, ecount, elist, &lp);
        if (rval) {
            fprintf (stderr, "add_connect failed\n"); goto CLEANUP;
        }
    }

    rval = CO759lp_opt (&lp, &infeasible);
    if (rval) {
        fprintf (stderr, "CO759lp_opt failed\n"); goto CLEANUP;
    }
    if (infeasible) {
        fprintf (stderr, "LP is infeasible\n"); 
        rval = 1; goto CLEANUP;
    }

    rval = CO759lp_objval (&lp, &objval);
    if (rval) {
        fprintf (stderr, "CO759lp_objval failed\n"); goto CLEANUP;
    }

    printf ("Current LP Value: %f\n", objval);
    fflush (stdout);

    rval = CO759lp_x (&lp, x);
    if (rval) {
        fprintf (stderr, "CO759lp_x failed\n"); goto CLEANUP;
    }
 
    for (i = 0, j = 0; j < ecount; j++) {
        if (x[j] > LP_EPSILON) i++;
    }

    printf ("Current LP graph has %d edges\n", i);
    for (j = 0; j < ecount; j++) {
        if (x[j] > LP_EPSILON) {
            printf ("%d %d %f\n", elist[2*j], elist[2*j+1], x[j]);
        }
    }
    fflush (stdout);

    for (i = 0, j = 0; j < ecount; j++) {
        if (x[j] > LP_EPSILON && x[j] < 1.0 - LP_EPSILON) break;
    }

    if (j == ecount) {
        printf ("LP solution is an optimal TSP tour\n");
    } else {
        printf ("Can use edge %d as branching variable\n", j);
    }


CLEANUP:
    CO759lp_free (&lp);
    if (x) free (x);
    return rval;
}

static int add_all_subtours (int ncount, int ecount, int *elist, CO759lp *lp)
{
    int rval = 0, i, sz, maxsz = ncount / 2, deltacount = 0;
    int *Set = (int *) NULL, *delta = (int *) NULL, *marks = (int *) NULL;
    int cnt = 0;

    /* To demonstrate the code for adding new constraints, this */
    /* routine will add ALL subtour inequalities.  This must be */
    /* replaced with the routine to add only violated subtours. */

    delta = (int *) malloc (ecount * sizeof(int));
    marks = (int *) malloc (ncount * sizeof(int));
    Set   = (int *) malloc (ncount * sizeof(int));
    if (!delta || !marks || !Set) {
        fprintf (stderr, "out of memory for delta, marks, or Set\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) marks[i] = 0;

    /* Run through all subsets of size sz from ncount elements */

    for (sz = 3; sz <= maxsz; sz++) {
        for (i = 0; i < sz; i++) Set[i] = i;
        for (; Set[sz-1] < ncount; next_set(sz,Set)) {
            get_delta (sz, Set, ecount, elist, &deltacount, delta, marks);
            rval = add_subtour (lp, deltacount, delta);
            if (rval) {
                fprintf (stderr, "add_subtour failed"); goto CLEANUP;
            }
            cnt++;
        }
    }

    printf ("Added all %d subtour constraints\n", cnt);
    fflush (stdout);

CLEANUP:
    if (delta) free (delta);
    if (marks) free (marks);
    if (Set) free (Set);
    return rval;
}

static void next_set (int sz, int *Set)
{
   int i;
   for (i=0; i < sz-1 && Set[i]+1 == Set[i+1]; i++) Set[i] = i;
   Set[i] = Set[i]+1;
}

static int add_connect (int ncount, int ecount, int *elist, CO759lp *lp)
{
    int rval = 0, icount, *island = (int *) NULL, *delta  = (int *) NULL;
    int round = 0, deltacount = 0, *marks = (int *) NULL;
    int infeasible, i;
    double *x = (double *) NULL, objval;
    graph G;

    init_graph (&G);

    rval = CO759lp_opt (lp, &infeasible);
    if (rval) { fprintf (stderr, "CO759lp_opt failed\n"); goto CLEANUP; }
    if (infeasible) {
        fprintf (stderr, "LP is infeasible\n"); rval = 1; goto CLEANUP;
    }

    rval = build_graph (ncount, ecount, elist, &G);
    if (rval) { fprintf (stderr, "build_graph failed\n"); goto CLEANUP; }

    x = (double *) malloc (ecount * sizeof (double));
    island = (int *) malloc (ncount * sizeof (int));
    delta  = (int *) malloc (ecount * sizeof(int));
    marks  = (int *) malloc (ncount * sizeof(int));
    if (!x || !island || !delta || !marks) {
        fprintf (stderr, "out of memory for x, island, delta, or marks\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < ncount; i++) marks[i] = 0;

    rval = CO759lp_x (lp, x);
    if (rval) { fprintf (stderr, "CO759lp_x failed\n"); goto CLEANUP; }

    while (connected (&G, x, &icount, island) == 0) {

        /*  just add one subtour; better to add one for each component */

        get_delta (icount, island, ecount, elist, &deltacount, delta, marks);

        rval = add_subtour (lp, deltacount, delta);
        if (rval) { fprintf (stderr, "add_subtour failed"); goto CLEANUP; }

        rval = CO759lp_opt (lp, &infeasible);
        if (rval) { fprintf (stderr, "CO759lp_opt failed\n"); goto CLEANUP; }
        if (infeasible) {
            fprintf (stderr, "LP is infeasible\n"); rval = 1; goto CLEANUP;
        }

        rval = CO759lp_objval (lp, &objval);
        if (rval) { fprintf (stderr, "CO759lp_objval failed\n"); goto CLEANUP; }

        printf ("Round %d LP: %f  (added subtour of size %d)\n",
                 round++, objval, icount); 
        fflush (stdout);

        rval = CO759lp_x (lp, x);
        if (rval) { fprintf (stderr, "CO759lp_x failed\n"); goto CLEANUP; }
    }

CLEANUP:
    free_graph (&G);
    if (x) free (x);
    if (island) free (island);
    if (delta) free (delta);
    if (marks) free (marks);
    return rval;
}

static int connected (graph *G, double *x, int *icount, int *island)
{
    int i;

    *icount = 0;
    for (i = 0; i < G->ncount; i++) G->nodelist[i].mark = 0;

    dfs (0, G, x, icount, island);

    if (*icount == G->ncount) return 1;
    else return 0;
}

static void dfs (int n, graph *G, double *x, int *icount, int *island)
{
    int i, neighbor;
    node *pn;

    island[*icount] = n;
    (*icount)++;

    pn = &G->nodelist[n];
    pn->mark = 1;

    for (i = 0; i < pn->deg; i++) {
        if (x[pn->adj[i].e] > LP_EPSILON) {
            neighbor = pn->adj[i].n;
            if (G->nodelist[neighbor].mark == 0) {
                dfs (neighbor, G, x, icount, island);
            }
        }
    }
}

static void init_graph (graph *G)
{
    if (G) {
        G->nodelist = (node *) NULL;
        G->adjspace = (adjobj *) NULL;
        G->ncount = 0;
        G->ecount = 0;
    }
}

static void free_graph (graph *G)
{
    if (G) {
        if (G->nodelist) free (G->nodelist);
        if (G->adjspace) free (G->adjspace);
    }
}

static int build_graph (int ncount, int ecount, int *elist, graph *G)
{
    int rval = 0, i, a, b;
    node *n;
    adjobj *p;

    G->nodelist = (node *) malloc (ncount * sizeof (node));
    G->adjspace = (adjobj *) malloc (2 * ecount * sizeof (node));
    if (!G->nodelist || !G->adjspace) {
        fprintf (stderr, "out of memory for nodelist or adjspace\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < ncount; i++) G->nodelist[i].deg = 0;
    for (i = 0; i < ecount; i++) {
        a = elist[2*i];  b = elist[2*i+1];
        G->nodelist[a].deg++;
        G->nodelist[b].deg++;
    }

    p = G->adjspace;
    for (i = 0; i < ncount; i++) {
        G->nodelist[i].adj = p;
        p += G->nodelist[i].deg;
        G->nodelist[i].deg = 0;
    }

    for (i = 0; i < ecount; i++) {
        a = elist[2*i];  b = elist[2*i+1];
        n = &G->nodelist[a];
        n->adj[n->deg].n = b;
        n->adj[n->deg].e = i;
        n->deg++;
        n = &G->nodelist[b];
        n->adj[n->deg].n = a;
        n->adj[n->deg].e = i;
        n->deg++;
    }

    G->ncount = ncount;
    G->ecount = ecount;

CLEANUP:
    return rval;
}

static void get_delta (int nsize, int *nlist, int ecount, int *elist,
       int *deltacount, int *delta, int *marks)
{
    int i, k = 0;

    for (i = 0; i < nsize; i++) marks[nlist[i]] = 1;

    for (i = 0; i < ecount; i++) {
        if (marks[elist[2*i]] + marks[elist[2*i+1]] == 1) {
            delta[k++] = i;
        }
    }
    *deltacount = k;

    for (i = 0; i < nsize; i++) marks[nlist[i]] = 0;
}

static int add_subtour (CO759lp *lp, int deltacount, int *delta)
{
    int rval = 0, i, newrows = 1, newnz = deltacount, *rmatind = delta;
    int rmatbeg[1];
    char sense[1];
    double rhs[1], *rmatval = (double *) NULL;

    rmatbeg[0] = 0; /* info for row starts at position 0 */
    rhs[0] = 2.0;   /* right-hand-side of subtour is 2.0 */
    sense[0] = 'G';    /* it is a >= inequality */

    rmatval = (double *) malloc (deltacount * sizeof (double));
    if (!rmatval) {
        fprintf (stderr, "out of memory for rmatval\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < deltacount; i++) rmatval[i] = 1.0;

    rval = CO759lp_addrows (lp, newrows, newnz, rhs, sense, rmatbeg,
                            rmatind, rmatval);
    if (rval) {
        fprintf (stderr, "CO759lp_addrows failed: %d\n", rval);
        goto CLEANUP;
    }

CLEANUP:
    if (rmatval) free (rmatval);
    return rval;
}

static int getprob (char *filename, int *p_ncount, int *p_ecount, int **p_elist,
    int **p_elen)
{
    FILE *f = (FILE *) NULL;
    int i, j, end1, end2, w, rval = 0, ncount, ecount;
    int *elist = (int *) NULL, *elen = (int *) NULL;
    double *x = (double *) NULL, *y = (double *) NULL;

    if (filename) {
        if ((f = fopen (filename, "r")) == NULL) {
    	    fprintf (stderr, "Unable to open %s for input\n",filename);
            rval = 1;  goto CLEANUP;
        }
    }

    if (filename && geometric_data == 0) {
        if (fscanf (f, "%d %d", &ncount, &ecount) != 2) {
       	    fprintf (stderr, "Input file %s has invalid format\n",filename);
            rval = 1;  goto CLEANUP;
        }

        printf ("Nodes: %d  Edges: %d\n", ncount, ecount);
        fflush (stdout);

        elist = (int *) malloc (2 * ecount * sizeof (int));
        if (!elist) {
            fprintf (stderr, "out of memory for elist\n");
            rval = 1;  goto CLEANUP;
        }

        elen = (int *) malloc (ecount * sizeof (int));
        if (!elen) {
            fprintf (stderr, "out of memory for elen\n");
            rval = 1;  goto CLEANUP;
        }

        for (i = 0; i < ecount; i++) {
    	    if (fscanf(f,"%d %d %d",&end1, &end2, &w) != 3) {
	        fprintf (stderr, "%s has invalid input format\n", filename);
                rval = 1;  goto CLEANUP;
	    }
	    elist[2*i] = end1;
	    elist[2*i+1] = end2;
	    elen[i] = w;
        }
    } else {
        if (filename) {
            if (fscanf (f, "%d", &ncount) != 1) {
       	        fprintf (stderr, "Input file %s has invalid format\n",filename);
                rval = 1;  goto CLEANUP;
            }
        } else {
            ncount = ncount_rand;
        }

        x = (double *) malloc (ncount * sizeof (double));
        y = (double *) malloc (ncount * sizeof (double));
        if (!x || !y) {
            fprintf (stdout, "out of memory for x or y\n");
            rval = 1; goto CLEANUP;
        }

        if (filename) {
            for (i = 0; i < ncount; i++) {
    	        if (fscanf(f,"%lf %lf",&x[i], &y[i]) != 2) {
	            fprintf (stderr, "%s has invalid input format\n", filename);
                    rval = 1;  goto CLEANUP;
	        }
            }
        } else {
            rval = CO759_build_xy (ncount, x, y, gridsize_rand);
            if (rval) {
                fprintf (stderr, "CO759_build_xy failed\n");
                goto CLEANUP;
            }
    
            printf ("%d\n", ncount);
            for (i = 0; i < ncount; i++) {
                printf ("%.0f %.0f\n", x[i], y[i]);
            }
            printf ("\n");
        }

        ecount = (ncount * (ncount - 1)) / 2;
        printf ("Complete graph: %d nodes, %d edges\n", ncount, ecount);

        elist = (int *) malloc (2 * ecount * sizeof (int));
        if (!elist) {
            fprintf (stderr, "out of memory for elist\n");
            rval = 1;  goto CLEANUP;
        }

        elen = (int *) malloc (ecount * sizeof (int));
        if (!elen) {
            fprintf (stderr, "out of memory for elen\n");
            rval = 1;  goto CLEANUP;
        }

        ecount = 0;
        for (i = 0; i < ncount; i++) {
            for (j = i+1; j < ncount; j++) {
                elist[2*ecount] = i;
                elist[2*ecount+1] = j;
                elen[ecount] = euclid_edgelen (i, j, x, y);
                ecount++;
            }
        }
    }

    *p_ncount = ncount;
    *p_ecount = ecount;
    *p_elist = elist;
    *p_elen = elen;

CLEANUP:
    if (f) fclose (f);
    if (x) free (x);
    if (y) free (y);
    return rval;
}

static int euclid_edgelen (int i, int j, double *x, double *y)
{
    double t1 = x[i] - x[j], t2 = y[i] - y[j];
    return (int) (sqrt (t1 * t1 + t2 * t2) + 0.5);
}

static int parseargs (int ac, char **av)
{
    int c;

    if (ac == 1) {
        usage (av[0]);
        return 1;
    }

    while ((c = getopt (ac, av, "ab:gk:s:")) != EOF) {
        switch (c) {
        case 'a':
            use_all_subtours = 1;
            break;
        case 'b':
            gridsize_rand = atoi (optarg); 
            break;
        case 'g':
            geometric_data = 1;
            break;
        case 'k':
            ncount_rand = atoi (optarg);
            break;
        case 's':
            seed = atoi (optarg);
            break;
        case '?':
        default:
            usage (av[0]);
            return 1;
        }
    }

    if (optind < ac) fname = av[optind++];

    if (optind != ac) {
        usage (av[0]);
        return 1;
    }

    return 0;
}

static void usage (char *f)
{
    fprintf (stderr, "Usage: %s [-see below-] [prob_file]\n", f);
    fprintf (stderr, "   -a    add all subtours cuts at once\n");
    fprintf (stderr, "   -b d  gridsize d for random problems\n");
    fprintf (stderr, "   -g    prob_file has x-y coordinates\n");
    fprintf (stderr, "   -k d  generate problem with d cities\n");
    fprintf (stderr, "   -s d  random seed\n");
}

