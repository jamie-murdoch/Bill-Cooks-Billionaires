#include<fstream>
#include<iostream>
#include<assert.h>

#include "kdtree.h"

int findmaxspread(int l, int u){
  int x_min = perm[l][0];
  int x_max = perm[l][0];
  int y_min = perm[l][1];
  int y_max = perm[l][1];
  for (int i = l; i <= u; i++){
    if (perm[i][0] < x_min)
      x_min = perm[i][0];
    else if (perm[i][0] > x_max)
      x_max = perm[i][0];

    if (perm[i][1] < y_min)
      y_min = perm[i][1];
    else if (perm[i][1] > y_max)
      y_max = perm[i][1];
  }

  if ((x_max - x_min) > (y_max - y_min))
    return 0;
  else return 1;
}

inline bool x_less(vector<int> p1, vector<int> p2){
  return p1[0] < p2[0];
}

inline bool y_less(vector<int> p1, vector<int> p2){
  return p1[1] < p2[1];
}

void select(int l, int u, int dim){
  if (dim == 0)
    sort(perm.begin() + l, perm.begin() + u + 1, x_less);
  else
    sort(perm.begin() + l, perm.begin() + u + 1, y_less);
}

kdnode *build(int l, int u){
  p = new kdnode;
  int m = (l + u) / 2;
  if (u - l + 1 <= BUCKETSIZE){
    p->bucket = 1;
    p->lopt = l;
    p->hipt = u;
  } else {
    p->bucket = 0;
    p->cutdim = findmaxspread(l, u);
    select(l, u, p->cutdim);
    p->cutval = perm[m][p->cutdim];
    p->loson = build(l, m);
    p->hison = build(m+1, u);
  }
  return p;
}

