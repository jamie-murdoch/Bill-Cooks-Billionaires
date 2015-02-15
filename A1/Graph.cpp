#include <fstream>
#include <iostream>

#include "Graph.h"


Edge::Edge(int e0, int e1, int _len) {
    end[0] = e0;
    end[1] = e1;
    len = _len;
}

bool Edge::ptr_compare(Edge *e0, Edge *e1) {
    return *e0 < *e1;
}

//Load a graph from a file
Graph::Graph(const char *filename) {
    // fstream fin;
    // fin.open(filename);

    // fin >> node_count >> edge_count;

    // edges.reserve(edge_count);

    // Edge edge;
    // for(int i = 0; i < edge_count; i++) {
    //     fin >> edge.end[0] >> edge.end[1] >> edge.len;
    //     edges.push_back(edge);
    // }


    //-------------------
    FILE *f = (FILE *) NULL;
    int i, end1, end2, w, rval = 0, ncount, ecount;

    if ((f = fopen (filename, "r")) == NULL) {
    fprintf (stderr, "Unable to open %s for input\n",filename);
        rval = 1;  goto CLEANUP;
    }

    if (fscanf (f, "%d %d", &ncount, &ecount) != 2) {
    fprintf (stderr, "Input file %s has invalid format\n",filename);
        rval = 1;  goto CLEANUP;
    }

    printf ("Nodes: %d  Edges: %d\n", ncount, ecount);
    fflush (stdout);

    edge_count = ecount;
    node_count = ncount;
    edges.resize(ecount);
    for (i = 0; i < ecount; i++) {
    if (fscanf(f,"%d %d %d",&end1, &end2, &w) != 3) {
        fprintf (stderr, "%s has invalid input format\n",filename);
            rval = 1;  goto CLEANUP;
    }
    edges[i].end[0] = end1;
    edges[i].end[1] = end2;
    edges[i].len = w;
    }

CLEANUP:
    if (f) fclose (f);
}

unsigned long Graph::sum_edge_weights(vector<int> &edge_indices) {
    unsigned long tot = 0;
    for(int i = 0; i < (int)edge_indices.size(); i++) {
        tot += edges[edge_indices[i]].len;
    }

    return tot;
}

void Graph::print_edges() {
    for (int i = 0; i < edges.size(); ++i)
    {
        cout << edges[i].end[0] << ", " << edges[i].end[1] << endl;
    }
}


Node::Node(){
    parent = this;
    rank = 0;
}


Node* Node::find_canonical() {
    Node* p;
    for (p = this; p->parent != p; p = p->parent);

    return p;
}

Node* Node::find_canonical_with_compression() {
    if(parent == this) {
        return this;
    }
    else {
        Node* new_parent = parent->find_canonical_with_compression();
        parent = new_parent;
        return parent;
    }
}

Node* Node::link(Node* x, Node* y){
    if (x->rank > y->rank) {
        Node::swap(x, y);
    }
    else if (x->rank == y->rank) {
        y->rank += 1;  
    }

    x->parent = y;

    return y;
}

void Node::swap(Node*& x, Node*& y){
    Node* temp = x;
    x = y;
    y = temp;
}
