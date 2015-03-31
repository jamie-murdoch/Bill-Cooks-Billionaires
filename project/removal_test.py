import subprocess
import sys, os
import ntpath

def main():
	problem_file = sys.argv[1]
	node_count = 0

	useless_fname = run_solvers(problem_file)
	
	solution_edges, node_count, node_count = load_edges("out.sol")
	useless_edges, node_count, useless_edge_count = load_edges(useless_fname)
	total_edge_count = node_count * node_count

	print
	conflict = check_for_conflicts(solution_edges, useless_edges)

	if not conflict:
		print "Test passed! No Conflicts found."
	else:
		print "Try again!"

	print
	print "Removed a total of %d edges out of %d. (%%%f)" % (useless_edge_count, total_edge_count, useless_edge_count/float(total_edge_count) * 100)
	print "There are %d edges remaining." % (total_edge_count - useless_edge_count, )

def run_solvers(problem_file):
	prob_name = ntpath.basename(problem_file)

	subprocess.call(["./edge_deletion", problem_file])

	useless_fname = prob_name + "-useless.edg"
	pruned_fname = prob_name + "-pruned.edg"
	orig_fname = prob_name + "-orig.edg"

	subprocess.call(["./test/prob2tsp", pruned_fname])
	subprocess.call(["./test/concorde", "out.tsp"])

	subprocess.call(["./test/prob2tsp", orig_fname])
	subprocess.call(["./test/concorde", "-f", "out.tsp"]) #-f save as edge file "out.sol"

	return useless_fname

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

	return edges, int(node_count), int(edge_count)

def check_for_conflicts(solution_edges, useless_edges):
	found_conflict = False

	for edge in useless_edges:
		if edge in solution_edges:
			print "You messed up!"
			print "Removed edge %s that was found in the optimal tour!" % edge

			found_conflict = True

	return found_conflict

if __name__ == "__main__":
	main()




