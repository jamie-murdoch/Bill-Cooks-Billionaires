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
	TSP_Solver() : m_min_tour_value(INFINITY) {}

	bool find_min_tour(Graph &graph, vector<int> &tour_indices);	

private:
	void build_lp(CO759lp *lp, const Graph &graph);
	int subtour (CO759lp *lp, Graph &graph, vector<int> &tour_indices);
	void add_subtour_inequalities(CO759lp *lp, Graph &graph);
	void get_delta (int nsize, int *nlist, int edge_count, vector<Edge> &elist, int *deltacount, int *delta, int *marks);
	int add_subtour (CO759lp *lp, int deltacount, int *delta);
	int add_connect (int node_count, int edge_count, vector<Edge> &elist, CO759lp *lp);
	void init_graph (ComponentGraph *G);
	void free_graph (ComponentGraph *G);
	int build_graph (int node_count, int edge_count, vector<Edge> &elist, ComponentGraph *G);
	int connected (ComponentGraph *G, double *x, int *icount, int *island);
	void dfs (int n, ComponentGraph *G, double *x, int *icount, int *island);

	//Helper Functions
	int run_lp(CO759lp *lp);
	double get_obj_val(CO759lp *lp);
	double* get_edges(CO759lp *lp, int max_edges);
	int get_num_edges(CO759lp *lp, int max_edges);
	void print_num_edges(CO759lp *lp, int max_edges);

	vector<double> x;
	double m_min_tour_value;
};

#endif


