#include "nodeset.h"

void makeset(node* x){
  x->parent = x; x->rank = 0;
}

/*
node* find(node* x){
  node* p;
  for (p = x; p->parent != p; p = p->parent)
    ;
  return p;
}
*/
node* find(node* x){
  if (x != x->parent)
    x->parent = find(x->parent);
  return x->parent;
}

void swap(node*& x, node*& y){
  node* temp = x;
  x = y;
  y = temp;
}

node* link(node* x, node* y){
  if (x->rank > y->rank)
    swap(x, y);
  else if (x->rank == y->rank)
    y->rank += 1;
  x->parent = y;
  return y;
}
