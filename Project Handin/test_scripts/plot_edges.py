import matplotlib as m
import json
import os
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

#["pr20.tsp","a80.tsp", "a180.tsp","a280.tsp",
files = ["pr1002.tsp","u1060.tsp","vm1084.tsp","pcb1173.tsp","d1291.tsp","rl1304.tsp","rl1323.tsp","nrw1379.tsp","u1432.tsp","d1655.tsp","vm1748.tsp","u1817.tsp","rl1889.tsp","d2103.tsp","u2152.tsp","u2319.tsp","pr2392.tsp","pcb3038.tsp","fnl4461.tsp","brd14051.tsp","d15112.tsp","d18512.tsp"]

#Note: need to run this from test/results with nothing else present
def main():
	for fname in files:
		dat_fname = fname + "-results.txt"
		path = os.path.join("edges/", dat_fname)

		if os.path.isfile(path):
			tot_edges, step1_edges_rem, step2_edges_rem = get_edges_rem(path)
			print "%s, %d, %d, %d" % (fname , tot_edges, step1_edges_rem , step2_edges_rem)

def get_edges_rem(path):
	with open(path, 'r') as f:
		lbl, tot_edges = f.readline().split()
		tot_edges = int(tot_edges)

		lbl, step1_edges_rem = f.readline().split()
		step1_edges_rem = tot_edges - int(step1_edges_rem)

		lbl, step1_time = f.readline().split()
		step1_time = float(step1_time)

		lbl, step2_edges_rem = f.readline().split()
		step2_edges_rem = step1_edges_rem - int(step2_edges_rem)

		lbl, step2_time = f.readline().split()
		step2_time = float(step2_time)

	return tot_edges, step1_edges_rem, step2_edges_rem

if __name__ == "__main__":
	main()