#include <fstream>
#include <iostream>

#include "Graph.h"

//Load a graph from a file
Graph::Graph(const char *filename) {
    int edge_count;

    fstream fin;
    fin.open(filename);

    fin >> node_count >> edge_count;

    edges.resize(edge_count);

    for(int i = 0; i < edge_count; i++) {
        fin >> edges[i].end[0] >> edges[i].end[1] >> edges[i].len;
    }
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
        x->parent = y;
    }

    return y;
}

void Node::swap(Node*& x, Node*& y){
    Node* temp = x;
    x = y;
    y = temp;
}
