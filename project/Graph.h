#ifndef GRAPH_H
#define GRAPH_H

#include "algebra.h"
#include "kdtree.h"
#include <vector>
#include <cmath>

using namespace std;

struct Edge {
    Edge() : int_len(0), len(0.0), useless(false) { end[0] = 0; end[1] = 0; }
    Edge(int e0, int e1, double length);

    int end[2]; //Index into points vector in graph
    int int_len;
    double len;
    bool useless;
};

struct Graph {
    Graph() : kd_tree(NULL) {}
    Graph(const char *tsp_file);

    ~Graph() { }

    unsigned long sum_edge_weights(vector<int> &edge_indices);
    void print_edges();
    int count_useless();

    int node_count() { return points.size(); }
    int edge_count() { return edges.size(); }

    vector<Point2D> points;
    vector<Edge> edges;
    vector<vector<int> > int_lengths;
    vector<vector<double> > lengths;

    KdTree *kd_tree;

};


#endif
