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

typedef struct adjdobj {
    int n;   /* index of neighbor Node */
    int e;   /* index of adj joining neighbor */
} adjobj;

typedef struct Node {
    int degree;
    adjobj *adj_objs;
    int mark;
} Node;

typedef struct Graph {
    int node_count;
    int edge_count;
    Node *nodelist;
    adjobj *adjspace;
} Graph;

static void usage (char *f);
static int getprob(char *fname, int *p_ncount, int *p_ecount, int **p_elist, int **p_elen);
static int parseargs (int ac, char **av);
static int subtour_init (int node_count, int edge_count, int *elist, int *elen, int *tlist);
static int subtour (CO759lp *lp, int edge_count, int node_count, int *elist, int *elen, int *tlist);
static int euclid_edgelen (int i, int j, double *x, double *y);
static void get_delta (int nsize, int *nlist, int edge_count, int *elist, int *deltacount, int *delta, int *marks);
static int add_subtour (CO759lp *lp, int deltacount, int *delta);
static int add_connect (int node_count, int edge_count, int *elist, CO759lp *lp);
static void init_graph (Graph *G);
static void free_graph (Graph *G);
static int build_graph (int node_count, int edge_count, int *elist, Graph *G);
static int connected (Graph *G, double *x, int *icount, int *island);
static void dfs (int n, Graph *G, double *x, int *icount, int *island);

static char *fname = (char *) NULL;
static int seed = 0;
static int geometric_data = 1;
static int ncount_rand = 0;
static int gridsize_rand = 100;
static int use_all_subtours = 0;
double min_tour_value = INFINITY;
int *min_tour = (int *) NULL;

int old_main (int ac, char **av)
{
    int rval  = 0, node_count = 0, edge_count = 0, j;
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

    rval = getprob (fname, &node_count, &edge_count, &elist, &elen);
    if (rval) {
        fprintf (stderr, "getprob failed\n"); goto CLEANUP;
    }

    if (use_all_subtours && node_count > 20) {
        fprintf (stderr, "Too many nodes to add all subtours\n"); goto CLEANUP;
    }

    min_tour = (int *) malloc((node_count-1) * sizeof(int));
    tlist = (int *) malloc ((node_count)*sizeof (int));
    if (!tlist) {
        fprintf (stderr, "out of memory for tlist\n");
        rval = 1;  goto CLEANUP; 
    }

    szeit = CO759_zeit ();
    rval = subtour_init (node_count, edge_count, elist, elen, tlist);
    if (rval) {
        fprintf (stderr, "subtour failed\n");
        goto CLEANUP;
    }

    printf("Optimal tour:\n");
    for (j = 0; j < node_count; j++) {
      printf ("%d %d %f\n", elist[2*min_tour[j]], elist[2*min_tour[j]+1], 1.0);
    }
    printf("Optimal tour value: %f\n",min_tour_value);
    fflush (stdout);
    
    printf ("Running Time: %.2f seconds\n", CO759_zeit() - szeit);
    fflush (stdout);

CLEANUP:

    if (tlist) free (tlist);
    if (elist) free (elist);
    if (elen) free (elen);
    return rval;
}

static int subtour_init (int node_count, int edge_count, int *elist, int *elen, int *tlist)
{
    int rval = 0, i, j, infeasible = 0;
    double  objective_val[1], lower_bound[1], upper_bound[1], objval;
    int     cmatbeg[1], cmatind[2];
    double  cmatval[2];
    CO759lp lp;

    rval = CO759lp_init (&lp);
    if (rval) { fprintf (stderr, "CO759lp_init failed\n"); goto CLEANUP; }

    rval = CO759lp_create (&lp, "subtour");
    if (rval) { fprintf (stderr, "CO759lp_create failed\n"); goto CLEANUP; }

    /* Build a row for each degree equation */

    for (i = 0; i < node_count; i++) {
        rval = CO759lp_new_row (&lp, 'E', 2.0);
        if (rval) {
            fprintf (stderr, "CO759lp_new_row failed\n"); goto CLEANUP;
        }
    }

    /* Build a column for each edge of the Graph */

    cmatbeg[0] = 0;
    cmatval[0] = 1.0;
    cmatval[1] = 1.0;
    for (j = 0; j < edge_count; j++) {
        objective_val[0]     = (double) elen[j];
        lower_bound[0]      = 0.0;
        upper_bound[0]      = 1.0;
        cmatind[0] = elist[2*j];
        cmatind[1] = elist[2*j+1];
        rval = CO759lp_addcols (&lp, 1 /* # of new variables */,
           2 /* # of new nonzeros */, objective_val, cmatbeg, cmatind, cmatval, lower_bound, upper_bound);
        if (rval) {
            fprintf (stderr, "CO759lp_addcols failed\n"); goto CLEANUP;
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
    rval = subtour(&lp, edge_count, node_count, elist, elen, tlist);
CLEANUP:
    CO759lp_free (&lp);
    return rval;

}

static int subtour (CO759lp *lp, int edge_count, int node_count, int *elist, int *elen, int *tlist){
    int rval = 0, i, j, infeasible = 0;
    double objval, *x = (double *) NULL, mindist = 0;
    x = (double *) malloc (edge_count * sizeof (double));
    if (!x) {
        fprintf (stderr, "out of memory for x\n");
        rval = 1; goto CLEANUP;
    }

    rval = CO759lp_opt (lp, &infeasible); //TODO - Is this being run twice right at the start?  
    if (rval) {
        fprintf (stderr, "CO759lp_opt failed\n"); goto CLEANUP;
    }
    if (infeasible) {
      printf("LP is infeasible, exitting\n");
      goto CLEANUP;
    }
    
    rval = CO759lp_x (lp, x);
    if (rval) {
        fprintf (stderr, "CO759lp_x failed\n"); goto CLEANUP;
    }
 
    for (i = 0, j = 0; j < edge_count; j++) {
        if (x[j] > LP_EPSILON) i++;
    }

    printf ("LP Graph has %d edges\n", i);
    /*    for (j = 0; j < edge_count; j++) {
        if (x[j] > LP_EPSILON) {
            printf ("%d %d %f\n", elist[2*j], elist[2*j+1], x[j]);
        }
    }*/
    fflush (stdout);

    rval = add_connect (node_count, edge_count, elist, lp);
    if (rval) {
        fprintf (stderr, "add_connect failed\n"); goto CLEANUP;
    }

    rval = CO759lp_opt (lp, &infeasible);
    if (rval) {
        fprintf (stderr, "CO759lp_opt failed\n"); goto CLEANUP;
    }
    if (infeasible) {
        printf ("LP is infeasible, exitting\n"); 
        goto CLEANUP;
    }

    rval = CO759lp_objval (lp, &objval);
    if (rval) {
        fprintf (stderr, "CO759lp_objval failed\n"); goto CLEANUP;
    }

    printf ("Current LP Value: %f\n", objval);
    fflush (stdout);

    if (objval > min_tour_value){
      printf ("Current LP value is higher than min tour value, exitting\n"); fflush(stdout); goto CLEANUP;    
    }
    
    rval = CO759lp_x (lp, x);
    if (rval) {
        fprintf (stderr, "CO759lp_x failed\n"); goto CLEANUP;
    }
 
    for (i = 0, j = 0; j < edge_count; j++) {
        if (x[j] > LP_EPSILON) i++;
    }

    printf ("Current LP Graph has %d edges\n", i);
    /*    for (j = 0; j < edge_count; j++) {
        if (x[j] > LP_EPSILON) {
            printf ("%d %d %f\n", elist[2*j], elist[2*j+1], x[j]);
        }
	}*/
    fflush (stdout);

    for (i = 0, j = 0; j < edge_count; j++) {
        double m = fmin(x[j], 1 - x[j]);
        if (m > mindist){
	       mindist = m;
	       i = j;
        }
    }

    rval = CO759lp_objval (lp, &objval);
    if (rval) {
      fprintf (stderr, "CO759lp_objval failed\n"); goto CLEANUP;
    }

    if(objval > min_tour_value){
      printf ("Current LP value is higher than min tour value, exitting\n"); fflush(stdout); goto CLEANUP;    
    }
    
    if (mindist < LP_EPSILON) {
        printf ("LP solution is an optimal TSP tour\n");
    	if (objval < min_tour_value){
    	  printf("NEW OPTIMAL TOUR VALUE: %f\n", objval);
    	  min_tour_value = objval;
    	  for (i = 0, j = 0; j < edge_count; j++){
    	    if (x[j] > LP_EPSILON){
    	      min_tour[i++] = j;
    	    }
    	  }
    	  if(i != node_count){
    	    printf ("Computed tour isn't a tour\n"); goto CLEANUP;
    	  }
    	}
    } else {
        printf ("Branching on edge %d\n", i);
        rval = CO759lp_setbnd(lp, i, 'U', 0.0);
    	if(rval){
    	  fprintf(stderr, "CO759lp_setbnd failed\n"); goto CLEANUP;
    	}
    	rval = subtour(lp, edge_count, node_count, elist, elen, tlist);
    	if(rval) goto CLEANUP;
    	rval = CO759lp_setbnd(lp, i, 'U', 1.0);
    	if(rval){
    	  fprintf(stderr, "CO759lp_setbnd failed\n"); goto CLEANUP;
    	}
    	rval = CO759lp_setbnd(lp, i, 'L', 1.0);
    	if(rval){
    	  fprintf(stderr, "CO759lp_setbnd failed\n"); goto CLEANUP;
    	}
    	rval = subtour(lp, edge_count, node_count, elist, elen, tlist);
    	if(rval) goto CLEANUP;
    	rval = CO759lp_setbnd(lp, i, 'L', 0.0);
    	if(rval){
    	  fprintf(stderr, "CO759lp_setbnd failed\n"); goto CLEANUP;
    	}
    }


CLEANUP:
    if (x) free (x);
    return rval;
}

static int add_connect (int node_count, int edge_count, int *elist, CO759lp *lp)
{
    int rval = 0, icount, *island = (int *) NULL, *delta  = (int *) NULL;
    int round = 0, deltacount = 0, *marks = (int *) NULL;
    int infeasible = 0, i;
    double *x = (double *) NULL, objval;
    Graph G;

    init_graph (&G);

    rval = CO759lp_opt (lp, &infeasible);
    if (rval) { fprintf (stderr, "CO759lp_opt failed\n"); goto CLEANUP; }
    if (infeasible) {
        fprintf (stderr, "LP is infeasible\n"); rval = 1; goto CLEANUP;
    }

    rval = build_graph (node_count, edge_count, elist, &G);
    if (rval) { fprintf (stderr, "build_graph failed\n"); goto CLEANUP; }

    x = (double *) malloc (edge_count * sizeof (double));
    island = (int *) malloc (node_count * sizeof (int));
    delta  = (int *) malloc (edge_count * sizeof(int));
    marks  = (int *) malloc (node_count * sizeof(int));
    if (!x || !island || !delta || !marks) {
        fprintf (stderr, "out of memory for x, island, delta, or marks\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < node_count; i++) marks[i] = 0;

    rval = CO759lp_x (lp, x);
    if (rval) { fprintf (stderr, "CO759lp_x failed\n"); goto CLEANUP; }

    while (connected (&G, x, &icount, island) == 0) {

        /*  just add one subtour; better to add one for each component */

        get_delta (icount, island, edge_count, elist, &deltacount, delta, marks);

        rval = add_subtour (lp, deltacount, delta);
        if (rval) { fprintf (stderr, "add_subtour failed"); goto CLEANUP; }

        rval = CO759lp_opt (lp, &infeasible);
        if (rval) { fprintf (stderr, "CO759lp_opt failed\n"); goto CLEANUP; }
        if (infeasible) {
            fprintf (stderr, "LP is infeasible, exitting\n"); goto CLEANUP;
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

static int connected (Graph *G, double *x, int *icount, int *island)
{
    int i;

    *icount = 0;
    for (i = 0; i < G->node_count; i++) G->nodelist[i].mark = 0;

    dfs (0, G, x, icount, island);

    if (*icount == G->node_count) return 1;
    else return 0;
}

static void dfs (int n, Graph *G, double *x, int *icount, int *island)
{
    int i, neighbor;
    Node *pn;

    island[*icount] = n;
    (*icount)++;

    pn = &G->nodelist[n];
    pn->mark = 1;

    for (i = 0; i < pn->degree; i++) {
        if (x[pn->adj_objs[i].e] > LP_EPSILON) {
            neighbor = pn->adj_objs[i].n;
            if (G->nodelist[neighbor].mark == 0) {
                dfs (neighbor, G, x, icount, island);
            }
        }
    }
}

static void init_graph (Graph *G)
{
    if (G) {
        G->nodelist = (Node *) NULL;
        G->adjspace = (adjobj *) NULL;
        G->node_count = 0;
        G->edge_count = 0;
    }
}

static void free_graph (Graph *G)
{
    if (G) {
        if (G->nodelist) free (G->nodelist);
        if (G->adjspace) free (G->adjspace);
    }
}

static int build_graph (int node_count, int edge_count, int *elist, Graph *G)
{
    int rval = 0, i, a, b;
    Node *n;
    adjobj *p;

    G->nodelist = (Node *) malloc (node_count * sizeof (Node));
    G->adjspace = (adjobj *) malloc (2 * edge_count * sizeof (Node));
    if (!G->nodelist || !G->adjspace) {
        fprintf (stderr, "out of memory for nodelist or adjspace\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < node_count; i++) G->nodelist[i].degree = 0;
    for (i = 0; i < edge_count; i++) {
        a = elist[2*i];  b = elist[2*i+1];
        G->nodelist[a].degree++;
        G->nodelist[b].degree++;
    }

    p = G->adjspace;
    for (i = 0; i < node_count; i++) {
        G->nodelist[i].adj_objs = p;
        p += G->nodelist[i].degree;
        G->nodelist[i].degree = 0;
    }

    for (i = 0; i < edge_count; i++) {
        a = elist[2*i];  b = elist[2*i+1];
        n = &G->nodelist[a];
        n->adj_objs[n->degree].n = b;
        n->adj_objs[n->degree].e = i;
        n->degree++;
        n = &G->nodelist[b];
        n->adj_objs[n->degree].n = a;
        n->adj_objs[n->degree].e = i;
        n->degree++;
    }

    G->node_count = node_count;
    G->edge_count = edge_count;

CLEANUP:
    return rval;
}

static void get_delta (int nsize, int *nlist, int edge_count, int *elist,
       int *deltacount, int *delta, int *marks)
{
    int i, k = 0;

    for (i = 0; i < nsize; i++) marks[nlist[i]] = 1;

    for (i = 0; i < edge_count; i++) {
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
    int i, j, end1, end2, w, rval = 0, node_count, edge_count;
    int *elist = (int *) NULL, *elen = (int *) NULL;
    double *x = (double *) NULL, *y = (double *) NULL;

    if (filename) {
        if ((f = fopen (filename, "r")) == NULL) {
    	    fprintf (stderr, "Unable to open %s for input\n",filename);
            rval = 1;  goto CLEANUP;
        }
    }

    if (filename && geometric_data == 0) {
        if (fscanf (f, "%d %d", &node_count, &edge_count) != 2) {
       	    fprintf (stderr, "Input file %s has invalid format\n",filename);
            rval = 1;  goto CLEANUP;
        }

        printf ("Nodes: %d  Edges: %d\n", node_count, edge_count);
        fflush (stdout);

        elist = (int *) malloc (2 * edge_count * sizeof (int));
        if (!elist) {
            fprintf (stderr, "out of memory for elist\n");
            rval = 1;  goto CLEANUP;
        }

        elen = (int *) malloc (edge_count * sizeof (int));
        if (!elen) {
            fprintf (stderr, "out of memory for elen\n");
            rval = 1;  goto CLEANUP;
        }

        for (i = 0; i < edge_count; i++) {
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
            if (fscanf (f, "%d", &node_count) != 1) {
       	        fprintf (stderr, "Input file %s has invalid format\n",filename);
                rval = 1;  goto CLEANUP;
            }
        } else {
            node_count = ncount_rand;
        }

        x = (double *) malloc (node_count * sizeof (double));
        y = (double *) malloc (node_count * sizeof (double));
        if (!x || !y) {
            fprintf (stdout, "out of memory for x or y\n");
            rval = 1; goto CLEANUP;
        }

        if (filename) {
            for (i = 0; i < node_count; i++) {
    	        if (fscanf(f,"%lf %lf",&x[i], &y[i]) != 2) {
	            fprintf (stderr, "%s has invalid input format\n", filename);
                    rval = 1;  goto CLEANUP;
	        }
            }
        } else {
            rval = CO759_build_xy (node_count, x, y, gridsize_rand);
            if (rval) {
                fprintf (stderr, "CO759_build_xy failed\n");
                goto CLEANUP;
            }
    
            printf ("%d\n", node_count);
            for (i = 0; i < node_count; i++) {
                printf ("%.0f %.0f\n", x[i], y[i]);
            }
            printf ("\n");
        }

        edge_count = (node_count * (node_count - 1)) / 2;
        printf ("Complete Graph: %d nodes, %d edges\n", node_count, edge_count);

        elist = (int *) malloc (2 * edge_count * sizeof (int));
        if (!elist) {
            fprintf (stderr, "out of memory for elist\n");
            rval = 1;  goto CLEANUP;
        }

        elen = (int *) malloc (edge_count * sizeof (int));
        if (!elen) {
            fprintf (stderr, "out of memory for elen\n");
            rval = 1;  goto CLEANUP;
        }

        edge_count = 0;
        for (i = 0; i < node_count; i++) {
            for (j = i+1; j < node_count; j++) {
                elist[2*edge_count] = i;
                elist[2*edge_count+1] = j;
                elen[edge_count] = euclid_edgelen (i, j, x, y);
                edge_count++;
            }
        }
    }

    *p_ncount = node_count;
    *p_ecount = edge_count;
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

