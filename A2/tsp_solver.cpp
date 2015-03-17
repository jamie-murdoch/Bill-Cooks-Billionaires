#include <iostream>
#include <getopt.h>
#include <math.h>
#include <cplex.h>

#include "tsp_solver.h"

//valgrind --dsymutil=yes ./subtour r20.dat

TSP_Solver::TSP_Solver(Graph &graph) : m_graph(graph), m_min_tour_value(INFINITY) {
    //Build the basic LP
    CO759lp_init (&m_lp);
    CO759lp_create (&m_lp, "subtour");

    /* Build a row for each degree equation */
    for(int i = 0; i < graph.node_count; i++) {
        CO759lp_new_row (&m_lp, 'E', 2.0);
    }

    /* Build a column for each edge of the Graph */
    int cmatbeg = 0, num_vars = 1, num_non_zero = 2;
    double coefficients[2] = {1.0, 1.0};
    double lower_bound = 0.0;
    double upper_bound = 1.0;
    for(int j = 0; j < graph.edge_count; j++) {
        int *nodes = (int*)graph.edges[j].end;
        double objective_val = (double)graph.edges[j].len;
        CO759lp_addcols (&m_lp, num_vars, num_non_zero, &objective_val,
                         &cmatbeg, nodes, coefficients, &lower_bound, &upper_bound);
    }

    CO759lp_write (&m_lp, "subtour.lp");

    //Moving some stuff for the DFS up here to save malloc calls
    init_graph (&G);
    islands = (int *) malloc (m_graph.node_count * sizeof (int));
    island_sizes = (int *) malloc (m_graph.node_count * sizeof (int));
    delta  = (int *) malloc (m_graph.edge_count * sizeof(int));
    edge_marks  = (int *) malloc (m_graph.node_count * sizeof(int));
    if (!islands || !delta || !edge_marks) {
        fprintf (stderr, "out of memory for x, islands, delta, or edge_marks\n");
    }
    build_graph (m_graph.node_count, m_graph.edge_count, m_graph.edges, &G);
}

TSP_Solver::~TSP_Solver() {
    CO759lp_free (&m_lp);

    free(islands);
    free(island_sizes);
    free(delta);
    free(edge_marks);

    free_graph(&G);
}

bool TSP_Solver::find_min_tour(vector<int> &tour_indices) {
    //Make sure we have room for the tour
    tour_indices.resize(m_graph.node_count, 0);

    //Use try to catch any errors from CPLEX
	try {
		//Run it the first time to check if it is feasible
		int infeasible = run_lp();
        if(infeasible) {
            cout << "Initial LP is infeasible." << endl;
            return false;
        } else {
            cout << "Degree-Equation LP Value: " << get_obj_val() << endl;    
        }
		
        //Run the main solver code
	    tsp_branch_and_bound(tour_indices);
	}
	catch(const char* error) {
		cerr << error << endl;
        return false;
	}

    return true;
}

int TSP_Solver::tsp_branch_and_bound(vector<int> &tour_indices) {
	double objval;
	int rval = 0, infeasible = 0;

    infeasible = run_lp();
    if(infeasible) {
        cout << "LP is infeasible, exitting." << endl;
        return rval;
    }
    print_num_edges();

    //Add constraints for disconnected components
	rval = add_subtour_inequalities();
    if(rval) return rval;

    //Run and check if the new LP is feasible
	infeasible = run_lp();
	if(infeasible) {
		cout << "LP is infeasible, exitting." << endl;
		return rval;
	}

    //Quit if LP val got worse
	objval = get_obj_val();
	cout << "Current LP Value: " << objval << endl;
    if (objval > m_min_tour_value){ //Can I mike this >=? Testing says yes...
      cout << "Current LP value is higher than min tour value, exitting" << endl;
      return rval;
    }
    
	print_num_edges();

	//Branch on edge with m_lp_edges[edge] closest to 0.5
	int branching_edge = compute_branch_edge();
    double branch_edge_val = m_lp_edges[branching_edge];

    if (is_almost_integral(branch_edge_val)) {
        cout << "LP solution is integral." << endl;

    	if (objval < m_min_tour_value){ //TODO do we really need this check? Only if objval == m_min_tour_val
    	    cout << "NEW OPTIMAL TOUR VALUE: " << objval << endl;
    	    m_min_tour_value = objval;

            //Update the tour
            if(!update_current_tour_indices(tour_indices)) {
                cout << "Computed tour isn't a tour" << endl;
                return rval;
            } 
    	}
    } else {
        cout << "Branching on edge " << branching_edge << endl;

        //Clamp it to 0  //TODO Try setting to 1.0 first
		CO759lp_setbnd(&m_lp, branching_edge, 'U', 0.0);

    	rval = tsp_branch_and_bound(tour_indices);
    	if(rval) return rval;

        //Clamp it to 1
		CO759lp_setbnd(&m_lp, branching_edge, 'U', 1.0);
		CO759lp_setbnd(&m_lp, branching_edge, 'L', 1.0);

    	rval = tsp_branch_and_bound(tour_indices);
    	if(rval) return rval;

        //Relax again
		CO759lp_setbnd(&m_lp, branching_edge, 'L', 0.0);
    }

    return rval;
}

int TSP_Solver::add_subtour_inequalities() {
    int rval = 0, icount;
    int round = 0, deltacount = 0;
    int infeasible = 0, i;
    double objval;


    get_edges();

    cout << "Looking for islands" << endl;

    int num_islands = 0;
    find_islands(islands, island_sizes, &num_islands);

    cout << "num islands: " << num_islands <<endl;
    // while(num_islands != 1) {

    //     cout << "Adding inequalities for " << num_islands << " islands" << endl;
    //     int island_index = 0;
    //     for(int i = 0; i < num_islands; i++) {
    //         int island_size = island_sizes[i];

    //         int delta_count = 0;
    //         get_delta(&islands[island_index], island_size, &delta_count, delta);

    //         add_subtour(deltacount, delta);     

    //         island_index += island_size;
    //     }
    //     cout << "Done adding subtour inequalities" << endl;

    //     int infeasible = run_lp();
    //     if(infeasible) return 1;

    //     get_edges();

    //     find_islands(islands, island_sizes, &num_islands);
    // }

    
    // while (!connected (&G, x, &icount, island, 0)) {

    //     get_delta (icount, island, m_graph.edge_count, m_graph.edges, &deltacount, delta, edge_marks);
    //     rval = add_subtour(deltacount, delta);

        
    //     for(int i = 1; i < m_graph.node_count; i++) {
    //         if(G.nodelist[i].mark == 0) { //Not yet traversed by 
    //             if(!connected(&G, x, &icount, island, i)) {
    //                 get_delta (icount, island, m_graph.edge_count, m_graph.edges, &deltacount, delta, edge_marks);
    //                 rval = add_subtour(deltacount, delta);
    //             }
    //             else {
    //                 break; //Graph is now connected
    //             }
    //         }
    //     }
        
    //     if (rval) { fprintf (stderr, "add_subtour failed"); goto CLEANUP; }

    //     rval = CO759lp_opt (&m_lp, &infeasible);
    //     if (rval) { fprintf (stderr, "CO759lp_opt failed\n"); goto CLEANUP; }
    //     if (infeasible) {
    //         fprintf (stderr, "LP is infeasible, exitting\n"); goto CLEANUP;
    //     }

    //     rval = CO759lp_objval (&m_lp, &objval);
    //     if (rval) { fprintf (stderr, "CO759lp_objval failed\n"); goto CLEANUP; }

    //     cout << "Round " << round++ << " LP: " << objval << "  (added subtour of size " << icount << ")" << endl;

    //     rval = CO759lp_x (&m_lp, x);
    //     if (rval) { fprintf (stderr, "CO759lp_x failed\n"); goto CLEANUP; }

    // }

CLEANUP:
    return rval;
}

int TSP_Solver::compute_branch_edge() {
    double *x = get_edges();
    double max_dist = 0.0;
    int branching_edge = 0;
    for (int edge = 0; edge < m_graph.edge_count; edge++) {
        double dist_from_int = min(x[edge], 1.0 - x[edge]);
        if (dist_from_int > max_dist){
           max_dist = dist_from_int;
           branching_edge = edge;
        }
    }

    return branching_edge;
}

//Return false if not a tour
bool TSP_Solver::update_current_tour_indices(vector<int> &tour_indices) {
    int index = 0;
    for (int j = 0; j < m_graph.edge_count; j++){
        if (m_lp_edges[j] > LP_EPSILON){
            //if(is_almost_integral(m_lp_edges[j])){ //TODO do I need to check this?
                tour_indices[index++] = j;    
            //}
        }
    }

    return index == m_graph.node_count;
}

// int TSP_Solver::connected (ComponentGraph *G, double *x, int *icount, int *island, int starting_node)
// {
//     int i;

//     *icount = 0;
//     for (i = 0; i < G->node_count; i++) G->nodelist[i].mark = 0;

//     dfs (starting_node, G, x, icount, island);

//     if (*icount == G->node_count) return 1;
//     else return 0;
// }


//new -- try using different marks for diffferent
//THIS IS BROKEN
void TSP_Solver::find_islands(int *islands, int *island_sizes, int *num_islands)
{
    for (int i = 0; i < G.node_count; i++) G.nodelist[i].mark = 0;

    int mark_colour = 1;
    int total_explored = 0;
    *num_islands = 0;

    for(int i = 0; i < G.node_count; i++) {
        int island_size = 0;
        dfs (i, &G, &m_lp_edges[0], &island_size, &islands[total_explored], mark_colour);
        
        if(island_size != 0) { //If we found a new island
            *num_islands++;
            total_explored += island_size;
            mark_colour++;

            if (total_explored == G.node_count) { //If we traversed all nodes quit
                cout << *num_islands <<endl;
                return;
            }
        } else {
            cout << "iszero" << endl;
        }
    }


}

void TSP_Solver::dfs (int n, ComponentGraph *G, double *x, int *icount, int *island, int mark_colour)
{
    int i, neighbor;
    GNode *pn;
    
    if(*icount >= m_graph.node_count) {
        cout << "interesting " << endl;
    }

    island[*icount] = n;
    (*icount)++;

    pn = &G->nodelist[n];
    pn->mark = mark_colour;

    for (i = 0; i < pn->degree; i++) {
        if (x[pn->adj_objs[i].e] > LP_EPSILON) {
            //cout <<
            neighbor = pn->adj_objs[i].n;
            if (G->nodelist[neighbor].mark == 0) {
                dfs (neighbor, G, x, icount, island, mark_colour);
            }
        }
    }
}

void TSP_Solver::get_delta(int *island, int island_size, int *delta_count, int *delta) {

    *delta_count = 0;

    for(int i = 0; i < island_size; i++) { //For every node in the island
        int node_id = island[i];
        GNode *node = &G.nodelist[node_id];
        int mark_colour = node->mark;

        for(int j = 0; j < node->degree; j++) { //Look through all neighboring nodes
            int neighb_node_id = node->adj_objs[j].n;
            int neighb_edge = node->adj_objs[j].e;

            GNode *neighb_node = &G.nodelist[neighb_node_id];
            int neighb_mark_colour = neighb_node->mark;

            //Neighbor in a different isalnd
            if(mark_colour != neighb_mark_colour) { //If they are marked differently, they aren't in the island
                delta[*delta_count] = neighb_edge;
                delta_count++;
            }
        }
    }
}

// void TSP_Solver::get_delta (int nsize, int *nlist, int edge_count, vector<Edge> &elist,
//        int *deltacount, int *delta, int *edge_marks)
// {
//     int i, k = 0;

//     //TODO Better way?
//     for(i = 0; i < m_graph.node_count; i++) edge_marks[i] = 0;

//     for (i = 0; i < nsize; i++) edge_marks[nlist[i]] = 1;

//     for (i = 0; i < edge_count; i++) {
//         if (edge_marks[elist[i].end[0]] + edge_marks[elist[i].end[1]] == 1) {
//             delta[k++] = i;
//         }
//     }
//     *deltacount = k;

//     for (i = 0; i < nsize; i++) edge_marks[nlist[i]] = 0;
// }

// void TSP_Solver::dfs (int n, ComponentGraph *G, double *x, int *icount, int *island)
// {
//     int i, neighbor;
//     GNode *pn;

//     island[*icount] = n;
//     (*icount)++;

//     pn = &G->nodelist[n];
//     pn->mark = 1;

//     for (i = 0; i < pn->degree; i++) {
//         if (x[pn->adj_objs[i].e] > LP_EPSILON) {
//             neighbor = pn->adj_objs[i].n;
//             if (G->nodelist[neighbor].mark == 0) {
//                 dfs (neighbor, G, x, icount, island);
//             }
//         }
//     }
// }

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



int TSP_Solver::add_subtour (int deltacount, int *delta)
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

    rval = CO759lp_addrows (&m_lp, newrows, newnz, rhs, sense, rmatbeg,
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
int TSP_Solver::run_lp() {
	int infeasible = 0;
	int rval = CO759lp_opt (&m_lp, &infeasible);

    if (rval) {
        throw "run_lp failed"; 
    }

    return infeasible;
}

double TSP_Solver::get_obj_val() {
	double objval;
	CO759lp_objval(&m_lp, &objval);
	return objval;
}

double* TSP_Solver::get_edges() {
    if(m_lp_edges.size() < m_graph.edge_count) m_lp_edges.resize(m_graph.edge_count);

    CO759lp_x (&m_lp, &m_lp_edges[0]);
    
    return &m_lp_edges[0];
}

int TSP_Solver::get_num_edges() {
    double *x = get_edges();

    int i = 0;
    for (int j = 0; j < m_graph.edge_count; j++) {
        if (x[j] > LP_EPSILON) i++;
    }

    return i;
}

void TSP_Solver::print_num_edges() {
	int num_tour_edges = get_num_edges();
	cout << "LP Graph has " << num_tour_edges << " edges" << endl;
}


  //   //for (i = 0; i < m_graph.node_count; i++) edge_marks[i] = 0;

  //   double *x = get_edges();

  // //  while (!connected (&G, x, &icount, island, 0)) {
  //       /*  just add one subtour; better to add one for each component */

  //     //  get_delta (icount, island, m_graph.edge_count, m_graph.edges, &deltacount, delta, edge_marks);
  //       //rval = add_subtour(deltacount, delta);

  //       //
  //       for(int i = 0; i < m_graph.node_count; i++) {
  //           if(G.nodelist[i].mark == 0) { //Not yet traversed by dfs
  //               if(!connected(&G, x, &icount, island, i)) {
  //                   get_delta (icount, island, m_graph.edge_count, m_graph.edges, &deltacount, delta, edge_marks);
  //                   rval = add_subtour(deltacount, delta);
  //               }
  //           }
  //       }
        
  //   //     if (rval) { fprintf (stderr, "add_subtour failed"); goto CLEANUP; }

  //   //     rval = CO759lp_opt (&m_lp, &infeasible);
  //   //     if (rval) { fprintf (stderr, "CO759lp_opt failed\n"); goto CLEANUP; }
  //   //     if (infeasible) {
  //   //         fprintf (stderr, "LP is infeasible, exitting\n"); goto CLEANUP;
  //   //     }

  //   //     rval = CO759lp_objval (&m_lp, &objval);
  //   //     if (rval) { fprintf (stderr, "CO759lp_objval failed\n"); goto CLEANUP; }

  //   //     cout << "Round " << round++ << " LP: " << objval << "  (added subtour of size " << icount << ")" << endl;

  //   //     rval = CO759lp_x (&m_lp, x);
  //   //     if (rval) { fprintf (stderr, "CO759lp_x failed\n"); goto CLEANUP; }

  //   // }