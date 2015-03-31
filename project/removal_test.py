import subprocess
import sys, os
import ntpath

def main():
	problem_file = sys.argv[1]
	prob_name = ntpath.basename(problem_file)

	subprocess.call(["./edge_deletion", problem_file])

	pruned_fname = prob_name + "-pruned.edg"
	orig_fname = prob_name + "-orig.edg"
	useless_fname = prob_name + "-useless.edg"

	subprocess.call(["./test/prob2tsp", pruned_fname])
	subprocess.call(["./test/concorde", "out.tsp"])

	subprocess.call(["./test/prob2tsp", orig_fname])
	subprocess.call(["./test/concorde", "-f", "out.tsp"]) #-f save as edge file "out.sol"

	solution_edges = load_edges("out.sol")
	useless_edges = load_edges(useless_fname)

	print useless_edges

def load_edges(fname):
	edges = {}

	with open(fname, 'r') as sol_file:
		node_count, edge_count = sol_file.readline().split()
		
		for line in sol_file:
			e0, e1 = line.split()[:2]
			
			#Canonical storage
			if(int(e0) > int(e1)):
				e0, e1 = e1, e0

			key = "%s,%s" % (e0, e1)
			edges[key] = True

	return edges

if __name__ == "__main__":
	main()