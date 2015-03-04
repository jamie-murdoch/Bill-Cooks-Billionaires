#include <iostream>
#include <getopt.h>
#include <math.h>
#include <cplex.h>

#include "find_tour.h"

//valgrind --dsymutil=yes ./subtour r20.dat

bool TSP_Solver::find_min_tour(Graph &graph, vector<int> &tour_indices) {
	bool success = true;

	tour_indices.resize(graph.node_count, 0);
	CO759lp lp;

	double *x;
	try {
		//Build the initial lp
		build_lp(&lp, graph);

		//Run it the first time to check if it is feasible
		run_lp(&lp);
		cout << "Degree-Equation LP Value: " << get_obj_val(&lp) << endl;

	    subtour(&lp, graph, tour_indices);

	}
	catch(const char* error) {
		cerr << error << endl;
		success = false;
	}


	//Cleanup
	CO759lp_free (&lp);
	return success;
}

void TSP_Solver::build_lp(CO759lp *lp, const Graph &graph) {
	CO759lp_init (lp);
	CO759lp_create (lp, "subtour");

	/* Build a row for each degree equation */
    for(int i = 0; i < graph.node_count; i++) {
    	CO759lp_new_row (lp, 'E', 2.0);
    }

    /* Build a column for each edge of the Graph */
    int cmatbeg = 0, num_vars = 1, num_non_zero = 2;
    double coefficients[2] = {1.0, 1.0};
    double lower_bound = 0.0;
	double upper_bound = 1.0;
    for(int j = 0; j < graph.edge_count; j++) {
    	int *nodes = (int*)graph.edges[j].end;
        double objective_val = (double)graph.edges[j].len;
        CO759lp_addcols (lp, num_vars, num_non_zero, &objective_val,
        	             &cmatbeg, nodes, coefficients, &lower_bound, &upper_bound);
    }

    CO759lp_write (lp, "subtour.lp");

    //TODO necessary?
    //Try runnint subtour outside of subtour_init in original
    run_lp(lp);
    double obj = get_obj_val(lp);

    //cout << "Done building LP" << endl;
}

int TSP_Solver::subtour (CO759lp *lp, Graph &graph, vector<int> &tour_indices) {
	double objval;
	int rval = 0, infeasible = 0;

	run_lp(lp); //TODO - Is this being run twice right at the start?  
	print_num_edges(lp, graph.edge_count);

	add_subtour_inequalities(lp, graph);

	infeasible = run_lp(lp);
	if(infeasible) {
		cout << "LP is infeasible, exitting." << endl;
		return rval;
	}

	objval = get_obj_val(lp);
	cout << "Current LP Value: " << objval << endl;

    if (objval > m_min_tour_value){
      cout << "Current LP value is higher than min tour value, exitting" << endl;
      return rval;
    }
    
	print_num_edges(lp, graph.edge_count);

	//Branch on edge with x[edge] closest to 0.5
	double *x = get_edges(lp, graph.edge_count);
	double max_dist = 0.0;
	int branching_edge = 0;
    for (int edge = 0; edge < graph.edge_count; edge++) {
        double dist_from_int = min(x[edge], 1.0 - x[edge]);
        if (dist_from_int > max_dist){
	       max_dist = dist_from_int;
	       branching_edge = edge;
        }
    }
    
    if (max_dist < LP_EPSILON) { //TODO should the indenting on this if be different?
        cout << "LP solution is integral." << endl;

    	if (objval < m_min_tour_value){
    	  cout << "NEW OPTIMAL TOUR VALUE: " << objval << endl;
    	  m_min_tour_value = objval;

    	  int index = 0;
    	  for (int j = 0; j < graph.edge_count; j++){
    	    if (x[j] > LP_EPSILON){
                if(index > graph.node_count - 1) {
                        cout << "wow" << endl;
                }
    	      tour_indices[index++] = j;
    	    }
    	  }

    	  if(index != graph.node_count){
    	    cout << "Computed tour isn't a tour" << endl;
    	    return rval;
    	  }
    	}

    } else {
        cout << "Branching on edge " << branching_edge << endl;
		
		CO759lp_setbnd(lp, branching_edge, 'U', 0.0);

    	rval = subtour(lp, graph, tour_indices);
    	if(rval) return rval;

		CO759lp_setbnd(lp, branching_edge, 'U', 1.0);
		CO759lp_setbnd(lp, branching_edge, 'L', 1.0);

    	subtour(lp, graph, tour_indices);
    	if(rval) return rval;

		CO759lp_setbnd(lp, branching_edge, 'L', 0.0);
    }

    return rval;
}

void TSP_Solver::add_subtour_inequalities(CO759lp *lp, Graph &graph)
{
    int rval = 0, icount;
    int round = 0, deltacount = 0;
    int infeasible = 0, i;
    double objval;
    ComponentGraph G;

    cout << "Adding subtour inequalities" << endl;

    init_graph (&G);

    infeasible = run_lp(lp); //TODO not necessary?
    if(infeasible) throw "Adding subtour inequalities failed.";
    cout <<"df" << endl;
    //exit(0);

    rval = build_graph (graph.node_count, graph.edge_count, graph.edges, &G);
    if (rval) { fprintf (stderr, "build_graph failed\n"); exit(1); }

    double *x = (double *) malloc (graph.edge_count * sizeof (double));
    int *island = (int *) malloc (graph.node_count * sizeof (int));
    int *delta  = (int *) malloc (graph.edge_count * sizeof(int));
    int *marks  = (int *) malloc (graph.node_count * sizeof(int));
    if (!x || !island || !delta || !marks) {
        fprintf (stderr, "out of memory for x, island, delta, or marks\n");
        rval = 1; goto CLEANUP;
    }
    for (i = 0; i < graph.node_count; i++) marks[i] = 0;

    rval = CO759lp_x (lp, x);
    if (rval) { fprintf (stderr, "CO759lp_x failed\n"); goto CLEANUP; }

    while (connected (&G, x, &icount, island) == 0) {
        /*  just add one subtour; better to add one for each component */

        get_delta (icount, island, graph.edge_count, graph.edges, &deltacount, delta, marks);

        rval = add_subtour (lp, deltacount, delta);
        if (rval) { fprintf (stderr, "add_subtour failed"); goto CLEANUP; }

        rval = CO759lp_opt (lp, &infeasible);
        if (rval) { fprintf (stderr, "CO759lp_opt failed\n"); goto CLEANUP; }
        if (infeasible) {
            fprintf (stderr, "LP is infeasible, exitting\n"); goto CLEANUP;
        }

        rval = CO759lp_objval (lp, &objval);
        if (rval) { fprintf (stderr, "CO759lp_objval failed\n"); goto CLEANUP; }

        // printf ("Round %d LP: %fL  (added subtour of size %d)\n",
        //          round++, objval, icount); 
        // fflush (stdout);
        cout << "Round " << round++ << " LP: " << objval << "  (added subtour of size " << icount << ")" << endl;

        rval = CO759lp_x (lp, x);
        if (rval) { fprintf (stderr, "CO759lp_x failed\n"); goto CLEANUP; }

    }

CLEANUP:
    free_graph (&G);
    if (x) free (x);
    if (island) free (island);
    if (delta) free (delta);
    if (marks) free (marks);
}

int TSP_Solver::connected (ComponentGraph *G, double *x, int *icount, int *island)
{
    int i;

    *icount = 0;
    for (i = 0; i < G->node_count; i++) G->nodelist[i].mark = 0;

    dfs (0, G, x, icount, island);

    if (*icount == G->node_count) return 1;
    else return 0;
}

void TSP_Solver::dfs (int n, ComponentGraph *G, double *x, int *icount, int *island)
{
    int i, neighbor;
    GNode *pn;

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

void TSP_Solver::init_graph (ComponentGraph *G)
{
    if (G) {
        G->nodelist = (GNode *) NULL;
        G->adjspace = (adjobj *) NULL;
        G->node_count = 0;
        G->edge_count = 0;
    }
}

void TSP_Solver::free_graph (ComponentGraph *G)
{
    if (G) {
        if (G->nodelist) free (G->nodelist);
        if (G->adjspace) free (G->adjspace);
    }
}

int TSP_Solver::build_graph (int node_count, int edge_count, vector<Edge> &elist, ComponentGraph *G)
{
    int rval = 0, i, a, b;
    GNode *n;
    adjobj *p;

    G->nodelist = (GNode *) malloc (node_count * sizeof (GNode));
    G->adjspace = (adjobj *) malloc (2 * edge_count * sizeof (GNode));
    if (!G->nodelist || !G->adjspace) {
        fprintf (stderr, "out of memory for nodelist or adjspace\n");
        rval = 1; goto CLEANUP;
    }

    for (i = 0; i < node_count; i++) G->nodelist[i].degree = 0;
    for (i = 0; i < edge_count; i++) {
        a = elist[i].end[0];  b = elist[i].end[1];
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
        a = elist[i].end[0];  b = elist[i].end[1];
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

void TSP_Solver::get_delta (int nsize, int *nlist, int edge_count, vector<Edge> &elist,
       int *deltacount, int *delta, int *marks)
{
    int i, k = 0;

    for (i = 0; i < nsize; i++) marks[nlist[i]] = 1;

    for (i = 0; i < edge_count; i++) {
        if (marks[elist[i].end[0]] + marks[elist[i].end[1]] == 1) {
            delta[k++] = i;
        }
    }
    *deltacount = k;

    for (i = 0; i < nsize; i++) marks[nlist[i]] = 0;
}

int TSP_Solver::add_subtour (CO759lp *lp, int deltacount, int *delta)
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


//Helpers//
int TSP_Solver::run_lp(CO759lp *lp) {
	int infeasible = 0;
	int rval = CO759lp_opt (lp, &infeasible);

    if (rval) {
        throw "run_lp failed"; 
    }

    return infeasible;
}

double TSP_Solver::get_obj_val(CO759lp *lp) {
	double objval;
	CO759lp_objval(lp, &objval);
	return objval;
}

double* TSP_Solver::get_edges(CO759lp *lp, int max_edges) {
    if(x.size() < max_edges) x.resize(max_edges);

    CO759lp_x (lp, &x[0]);
    
    return &x[0];
}

int TSP_Solver::get_num_edges(CO759lp *lp, int max_edges) {
    double *x = get_edges(lp, max_edges);

    int i = 0;
    for (int j = 0; j < max_edges; j++) {
        if (x[j] > LP_EPSILON) i++;
    }

    return i;
}

void TSP_Solver::print_num_edges(CO759lp *lp, int max_edges) {
	int num_tour_edges = get_num_edges(lp, max_edges);
	cout << "LP Graph has " << num_tour_edges << " edges" << endl;
}
