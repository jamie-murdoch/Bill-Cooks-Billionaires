#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include "nodeset.h"
#include "util.h"

using namespace std;

struct edge{
  int end[2];
  int len;
  bool operator<(const edge& val) const {
    return len < val.len;
  }
};

static int kruskal_tree(int ncount, int ecount, vector<edge>& elist,
			vector<int>& tlist)
{
  int i, j = 0;
  vector<node*> block;
  block.reserve(ncount);
  for (i = 0; i < ncount; i++){
    block[i] = new node;
    makeset(block[i]);
  }

  tlist.reserve(ncount - 1);
  for (i = 0; i < ncount - 1; i++) tlist[i] = -1;

  sort(elist.begin(), elist.end());

  for (i = 0; i < ecount; i++){
    if (find(block[elist[i].end[0]]) != find(block[elist[i].end[1]])){
      link(find(block[elist[i].end[0]]), find(block[elist[i].end[1]]));
      tlist[j++] = i;
    }
  }

  return 0;
}



int main(int argc, char **argv){
  int ncount; int ecount;
  vector<int> tlist;
  vector<edge> edgelist;
  edge myedge; 

  int i;
  int e0, e1, cost;

  ifstream fin;

  if(argc > 1) {
    fin.open(argv[1]);
  }
  else {
    cerr << "Pass input file as first argument." << endl;
    return 1;
  }

  if (fin >> ncount >> ecount){
    while(fin >> e0 >> e1 >> cost){
      myedge.end[0] = e0; myedge.end[1] = e1;
      myedge.len = cost;
      edgelist.push_back(myedge);
    }
  }
  fin.close();

  double szeit = CO759_zeit();
  kruskal_tree(ncount, ecount, edgelist, tlist);
  cout << "Running time " << CO759_zeit() - szeit << endl;

  unsigned long weight = 0;
  for (i = 0; i < ncount - 1; i++){
    weight += edgelist[tlist[i]].len;
  }
  cout << "mst has weight " << weight << endl;
  
}
