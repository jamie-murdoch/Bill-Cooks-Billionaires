#include <fstream>
#include <iostream>
#include <limits>

#include "Graph.h"
#include "util.h"


Edge::Edge(int e0, int e1, double length) {
    end[0] = e0;
    end[1] = e1;
    len = length;
    int_len = nint(len);
    useless = false;
}

//Load a graph from a file
Graph::Graph(const char *tsp_file) : kd_tree(NULL) {
    fstream fin;
    fin.open(tsp_file);

    string temp;
    getline(fin, temp); //Skip name
    getline(fin, temp); //Skip comment
    getline(fin, temp); //Skip type

    int num_points;
    fin >> temp >> temp >> num_points; //Read edge count
    getline(fin, temp); //Skip weight type
    getline(fin, temp); //Skip weight type
    getline(fin, temp); //Skip section type


    //test
    // num_points = 50;
    // vector<double> x(num_points);
    // vector<double> y(num_points);
    // CO759_build_xy(num_points, x, y, 100);
    // points.reserve(num_points);
    // for(int i = 0; i < num_points; i++) {
    //     int index;
    //     double x, y;
    //     fin >> index >> x >> y;
    //     points.push_back(Point2D(x,y));
    // }

    //
    points.reserve(num_points);
    for(int i = 0; i < num_points; i++) {
        int index;
        double x, y;
        fin >> index >> x >> y;
        points.push_back(Point2D(x,y));
    }

    fin.close();

    //Assume it is a complete graph
    edges.reserve(num_points * (num_points - 1) / 2);
    int_lengths.resize(num_points);
    lengths.resize(num_points);
    useless.resize(num_points);
    for(int i = 0; i < num_points; i++) {
        int_lengths[i].resize(num_points,  numeric_limits<int>::max());
        lengths[i].resize(num_points, numeric_limits<double>::infinity());
	useless[i].resize(num_points, false);
    }

    for(int i = 0; i < num_points; i++) {
        for(int j = i; j < num_points; j++) {
            if(i != j) {
                double len = (points[i] - points[j]).length();
                edges.push_back(Edge(i, j, len));
                int_lengths[i][j] = nint(len); int_lengths[j][i] = nint(len);
                lengths[i][j] = len; lengths[j][i] = len;
                // cout << num_points << endl;
                // cout << edge_count() << endl;
            }
        }
    }

    kd_tree = new KdTree(points, lengths);
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
        cout << "Edge " << i <<": " << "ends = (" << edges[i].end[0] << ", " << edges[i].end[1] << ") len = " << edges[i].int_len<< endl;
    }
}

int Graph::count_useless() {
    int sum = 0;
    for(int i = 0; i < (int)edges.size(); i++) {
        if(edges[i].useless) sum++;
    }
    return sum;
}

void Graph::get_bounding_box(double &minX, double &minY, double &maxX, double &maxY) const{
    minX = minY = numeric_limits<double>::infinity();
    maxX = maxY = -numeric_limits<double>::infinity();

    for(int i = 0; i < (int)points.size(); i++) {
        const Point2D &p = points[i];
        minX = min(p.x(), minX);
        minY = min(p.y(), minY);
        maxX = max(p.x(), maxX);
        maxY = max(p.y(), maxY);
    }
}


void Graph::save_edges(string fname, bool output_useless) {
    ofstream fout;
    fout.open(fname.c_str());

    fout << node_count() << " ";

    if(output_useless) {
        fout << edge_count() << endl;
    }
    else {
        fout << edge_count() - count_useless() << endl;
    }

    for(int i = 0; i < edge_count(); i++) {
        Edge &e = edges[i];
        if(!e.useless || output_useless) {
            fout << e.end[0] << " " << e.end[1] << " " << e.int_len << endl;
        }
    }

    fout.close();
}








