#ifndef GRAPH_H
#define GRAPH_H

#include <vector>

using namespace std;

struct Edge {
    Edge() {}
    Edge(int e0, int e1, int _len);

    int end[2];
    int rounded_len;
    double len;

    bool operator<(const Edge& val) const {
    return len < val.len;
    }

    static bool ptr_compare(Edge *e0, Edge *e1);
};

struct Graph {
    Graph() : node_count(0), edge_count(0) {}
    Graph (const char *filename);

    unsigned long sum_edge_weights(vector<int> &edge_indices);
    void print_edges();

    //Returns false if graph is connected
    //bool find_island(double *x_edges, vector<int> &components)

    int node_count;
    int edge_count;
    double *x; 
    double *y;
    vector<Edge> edges;
    ~Graph(){
      if(x) delete x;
      if(y) delete y;
    }
};

/*
 * Header file for node struct used to build sets of nodes
 * In Kruskal, a connected component is a set of nodes stored as a tree
 *  - each component has an arbitrary ``canonical node'', one that pts to itself
 *  - two nodes are in the same component if they have the same parent
 *  Operations:
 *  - making a singleton, finding parent, swapping, and linking
 */

// struct Node {
//     Node();

//     // finds the canonical element of Node  performs path-compression en route
//     Node* find_canonical();
//     Node* find_canonical_with_compression();

//     // x and y are the canonical nodes of two disjoint trees.
//     // assumes (or swaps, or increments so that) y has higher rank than x
//     // merges the trees by making y the parent of x
//     // returns a pointer to the merged tree
//     static Node* link(Node* x, Node* y);

//     // swaps the nodes x, y
//     static void swap(Node*& x, Node*& y);

//     int rank;
//     Node *parent;
// };

#endif
