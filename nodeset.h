/*
 * Header file for node struct used to build sets of nodes
 * In Kruskal, a connected component is a set of nodes stored as a tree
 *  - each component has an arbitrary ``canonical node'', one that pts to itself
 *  - two nodes are in the same component if they have the same parent
 *  Operations:
 *  - making a singleton, finding parent, swapping, and linking
 */

struct node {
  //int val;
  int rank;
  node *parent;
};

// creates a singleton out of x: a self-pointing, rank-0 node
void makeset(node* x);

// finds the canonical element of node x; performs path-compression en route
node* find(node* x);

// swaps the nodes x, y
void swap(node*& x, node*& y);

// x and y are the canonical nodes of two disjoint trees.
// assumes (or swaps, or increments so that) y has higher rank than x
// merges the trees by making y the parent of x
// returns a pointer to the merged tree
node* link(node* x, node* y);
