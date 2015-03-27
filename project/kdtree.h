#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include "algebra.h"

#define BUCKETSIZE 3

//used to initialize the nearest-neighbor search
//#define HUGE 9999

using namespace std;


struct KdNode {
  char bucket;
  int cutdim;
  float cutval;
  KdNode *loson, *hison;
  int lopt, hipt;
};

class KdTree {
 public:
  KdTree(const vector<Point2D> &points);

  float px(int i, int j); 
  int findmaxspread(int l, int u);
  bool pt_less(int i, int j, int dim);
  void select(int l, int u, int m, int dim);
  KdNode* build(int l, int u);
  void print_point(int i);
  void print_tree(KdNode *node);

 private:
  vector<int> mPerm;
  const vector<Point2D> *mPoints;
  KdNode *mRoot;
};

#endif
