#ifndef FIND_TOUR_H
#define FIND_TOUR_H

#include <vector>
#include "Graph.h"

#include "lp.h"
#include "util.h"


using namespace std;


typedef struct adjdobj {
    int n;   /* index of neighbor Node */
    int e;   /* index of adj joining neighbor */
} adjobj;

typedef struct GNode {
    int degree;
    adjobj *adj_objs;
    int mark;
} GNode;

typedef struct ComponentGraph {
    int node_count;
    int edge_count;
    GNode *nodelist;
    adjobj *adjspace;
} ComponentGraph;


class TSP_Solver {
public:
	TSP_Solver(Graph &graph);
	~TSP_Solver();

	bool find_min_tour(vector<int> &tour_indices);	

private:
	int tsp_branch_and_bound(vector<int> &tour_indices);
	int add_subtour_inequalities();
	int compute_branch_edge();
	bool update_current_tour_indices(vector<int> &tour_indices);


	void get_delta (int nsize, int *nlist, int edge_count, vector<Edge> &elist, int *deltacount, int *delta, int *marks);
	int add_subtour (int deltacount, int *delta);
	int add_connect (int node_count, int edge_count, vector<Edge> &elist);
	void init_graph (ComponentGraph *G);
	void free_graph (ComponentGraph *G);
	int build_graph (int node_count, int edge_count, vector<Edge> &elist, ComponentGraph *G);
	int connected (ComponentGraph *G, double *x, int *icount, int *island, int starting_node);
	void dfs (int n, ComponentGraph *G, double *x, int *icount, int *island);

	//Helper Functions
	int run_lp();
	double get_obj_val();
	double* get_edges();
	int get_num_edges();
	void print_num_edges();

	Graph &m_graph;
	CO759lp m_lp;
	vector<double> m_lp_edges;
	double m_min_tour_value;

	int *island;
	int *delta;
	int *edge_marks;
	ComponentGraph G;
};

#endif


