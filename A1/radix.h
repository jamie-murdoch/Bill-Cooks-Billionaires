#ifndef RADIX_H
#define RADIX_H

#include <iostream> /* cout */
#include <cstdlib> /* srand, rand */
#include <ctime> /* time() */
#include <list> /* list */
#include <vector>

#include "Graph.h"
 
#define BASE 10 // # of buckets to use
#define ARRAY_SIZE 50 // max # of elements in array
 
using namespace std;
 
// function prototypes
void radix(int* nums, int length, int max);

void radix(vector<Edge> &edges, int length, int max) {
	list<Edge> bucket[BASE];
	int i;
 
	// iterate through each radix until n>max
	for (int n=1; max >= n; n *= BASE) {
		// sort list of numbers into buckets
		for (i=0; i<length; i++)
			bucket[(edges[i].len/n)%BASE].push_back(edges[i]);
 
		// merge buckets back to list
		for (int k=i=0; i<BASE; bucket[i++].clear())
			for (list<Edge>::iterator j = bucket[i].begin(); j != bucket[i].end(); edges[k++] = *(j++));
	}
}

// void radix(int* nums, int length, int max) {
// 	list<int> bucket[BASE];
// 	int i;
 
// 	// iterate through each radix until n>max
// 	for (int n=1; max >= n; n *= BASE) {
// 		// sort list of numbers into buckets
// 		for (i=0; i<length; i++)
// 			bucket[(nums[i]/n)%BASE].push_back(nums[i]);
 
// 		// merge buckets back to list
// 		for (int k=i=0; i<BASE; bucket[i++].clear())
// 			for (list<int>::iterator j = bucket[i].begin(); j != bucket[i].end(); nums[k++] = *(j++));
// 	}
// }
#endif