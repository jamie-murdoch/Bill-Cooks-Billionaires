#include<iostream>
#include<vector>
#include<fstream>

#include "nodeset.h"

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
  int i;
  vector<node*> block;
  block.reserve(ncount);
  for (i = 0; i < ncount; i++){
    block[i] = new node;
    makeset(block[i]);
  }

  tlist.reserve(ncount - 1);
  for (i = 0; i < ncount - 1; i++) tlist[i] = 0;

  sort(elist.begin(), elist.end());

  for (i = 0; i < ecount; i++){
    if (find(block[elist[i].end[0]]) != find(block[elist[i].end[1]])){
      link(find(block[elist[i].end[0]]), find(block[elist[i].end[1]]));
      tlist[i] = 1;
    }
  }

  return 0;
}



int main(){
  int ncount; int ecount;
  vector<int> tlist;
  vector<edge> edgelist;
  edge myedge; 

  int i;
  int e0, e1, cost;

  ifstream fin;
  fin.open("g100.1012.txt");

  if (fin >> ncount >> ecount){
    while(fin >> e0 >> e1 >> cost){
      myedge.end[0] = e0; myedge.end[1] = e1;
      myedge.len = cost;
      edgelist.push_back(myedge);
    }
  }
  fin.close();
  
  kruskal_tree(ncount, ecount, edgelist, tlist);


  int weight = 0;
  for (i = 0; i < ecount - 1; i++){
    if (tlist[i]){
      weight += edgelist[i].len;
    }
  }
  cout << "mst has weight " << weight << endl;
  
}
