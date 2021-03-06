#include<fstream>
#include<iostream>
#include<cmath>

#include "kdtree.h"



void PointSet::addpt(Point newpt){
  ptvec.push_back(newpt);
}

//returns the jth coord of ptvec[i]
float PointSet::x(int i, int j){
  return (j == 0) ? ptvec[i].x : ptvec[i].y;
}

float PointSet::dist(int i, int j){
  return sqrt(pow((x(i, 0) - x(j, 0)), 2) + pow((x(i, 1) - x(j, 1)), 2));
}

int PointSet::numpts(void){
  return ptvec.size();
}

// builds a kdtree from a given pointset
KdTree::KdTree(PointSet *ptset){
  int n = ptset->numpts();
  perm.resize(n);
  source_ptset = ptset;

  for(int i = 0; i < n; i++)
    perm[i] = i;

  root = build(0, n - 1);
}

//accesses jth coord of perm[i]
float KdTree::px(int i, int j){
  return source_ptset->x(perm[i], j);
}

//prints the point perm[i] in the format (x, y)-
void KdTree::print_point(int i){
  cout << "(" << px(i, 0) << "," << px(i, 1) << ")-";
}

//comparison operator to test if perm[i] has dimth coordinate less than perm[j]
bool KdTree::pt_less(const int i, const int j, const int dim){
  return px(i, dim) < px(j, dim);
}

//searches perm[l...u] and returns int corresponding to dimension with greatest
// range of vals
int KdTree::findmaxspread(int l, int u){
  float x_min = px(l, 0); float x_max = px(l, 1);
  float y_min = px(l, 1); float y_max = px(l, 1);
  for(int i = l; i <= u; i++){
    if (px(i, 0) < x_min){
      x_min = px(i, 0);
    } else if (px(i, 0) > x_max){
      x_max = px(i, 0);
    }
    if (px(i, 1) < y_min){
      y_min = px(i, 1);
    } else if (px(i, 1) > y_max){
      y_max = px(i, 1);
    }
  }

  return ((x_max - x_min) > (y_max - y_min)) ? 0 : 1;
}


// rearranges perm[l...u] (in practice m = (l+u)/2) such that perm[m] has dimth
// coordinate no greater than any point to its left, and no less than any point to
// its right...without sorting the array! yay!
void KdTree::select(int l, int u, int m, int dim){
  int again = 1;
  int temp;
  while(again){
    again = 0;
    // tests if perm[m] has dimth coordinate greater than any pt to its left
    for(int i = l; i <= m; i++){
      if(pt_less(i, m, dim)){
	  temp = perm[m];
	  perm[m] = perm[i];
	  perm[i] = temp;
	  again = 1;
	}
    }
    //tests if perm[m] has dimth coordinate less than any pt to its right
    for(int i = m; i <= u; i++){
      if(pt_less(m, i, dim)){
	  temp = perm[m];
	  perm[m] = perm[i];
	  perm[i] = temp;
	  again = 1;
      }
    }
  }
}

// recursive function used to partition the points
KdNode* KdTree::build(int l, int u){
  KdNode *p = new KdNode;
  if (u - l + 1 <= BUCKETSIZE){
    p->bucket = 1;
    p->lopt = l;
    p->hipt = u;
  } else {
    p->bucket = 0;
    p->cutdim = findmaxspread(l, u);
    int m = (l + u) / 2;
    select(l, u, m, p->cutdim);
    p->cutval = px(m, p->cutdim);
    p->loson = build(l, m);
    p->hison = build(m + 1, u);
  }
  return p;
}

void KdTree::print_tree(KdNode *node){
  int low, high;
  if (node == NULL) node = root;
  if (node->bucket){
    cout << "bucket node:" << endl;
    low = node->lopt; high = node->hipt;
    cout << "bucket low and high range " << low << " " << high << endl;
    for (int i = low; i <= high; i++)
      print_point(i);
    cout << endl;
  } else {
    cout << "non-bucket w cut-dim " << node->cutdim << "and val " << node->cutval
	 << endl;
    print_tree(node->loson);
    print_tree(node->hison);
  }
}

int main(void){
  PointSet *ptset = new PointSet;
  int x, y;
  Point newpt;
  int num;
  cout << "How many points to enter ";
  cin >> num;
  for(; num > 0; num--){
    cout << "Enter two points for an x-y coord " << num << "remaining ";
    cin >> x >> y;
    newpt.x = x; newpt.y = y;
    ptset->addpt(newpt);
  }

  KdTree *kdt;
  kdt = new KdTree(ptset);
  kdt->print_tree(NULL);
  

  return 0;
}

