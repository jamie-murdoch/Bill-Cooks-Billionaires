import subprocess
from subprocess import Popen, PIPE
import multiprocessing
import json
import os
import re

files = ["a80.tsp", "a180.tsp","a280.tsp","pr1002.tsp","u1060.tsp","vm1084.tsp","pcb1173.tsp","d1291.tsp","rl1304.tsp","rl1323.tsp","nrw1379.tsp","u1432.tsp","d1655.tsp","vm1748.tsp","u1817.tsp","rl1889.tsp","d2103.tsp","u2152.tsp","u2319.tsp","pr2392.tsp","pcb3038.tsp","fnl4461.tsp","brd14051.tsp","d15112.tsp","d18512.tsp"]

N_SAMPLES = 1

#Note: need to run this from test/results with nothing else present
def main():
	pool = multiprocessing.Pool(4)

	pool.map(time_file, files)

	print "Done all problems."

def time_file(prob_name):
	print "Timing problem " + prob_name
	print

	output_fname = prob_name + "-times.txt"

	oldir = os.getcwd()
	working_dir = prob_name.split('.')[0]
	
	try:
		os.mkdir(working_dir)
	except Exception:
		if os.path.isfile(working_dir + "/" + prob_name + "-times.txt"):
			return None

	os.chdir(working_dir)

	#subprocess.call(["../../../edge_deletion", "../../../data/" + prob_name])

	pruned_fname = prob_name + "-pruned1+2.edg"
	orig_fname = prob_name + "-orig1+2.edg"

	pruned_time = 0.0
	orig_time = 0.0

	for i in xrange(N_SAMPLES):
		#Run on pruned graph
		subprocess.call(["../test/prob2tsp", pruned_fname])
		pruned_results = subprocess.check_output(["../test/concorde", "out.tsp"])

		#run on original graph
		subprocess.call(["../test/prob2tsp", orig_fname])
		orig_results = subprocess.check_output(["../test/concorde", "out.tsp"])

		pruned_time += get_runtime(pruned_results)
		orig_time += get_runtime(orig_results)

	#Save to disk
	data = {prob_name: {"pruned": pruned_time/N_SAMPLES, "orig": orig_time/N_SAMPLES}}
	with open(output_fname, 'w') as times_file:
		json.dump(data, times_file)

	os.chdir(oldir)
	return None

def get_runtime(result_str):
	lines = result_str.split('\n')
	time_str = lines[-2]
	try:
		return float(re.findall(r"[-+]?\d*\.\d+|\d+",time_str)[0])
	except:
		return time_str


if __name__ == "__main__":
	main()