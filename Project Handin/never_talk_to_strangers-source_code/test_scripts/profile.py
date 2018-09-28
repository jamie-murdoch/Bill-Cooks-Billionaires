import subprocess
import multiprocessing
#files = ["pr1002.tsp","u1060.tsp","vm1084.tsp","pcb1173.tsp","d1291.tsp","rl1304.tsp","rl1323.tsp","nrw1379.tsp","u1432.tsp","d1655.tsp","vm1748.tsp","u1817.tsp","rl1889.tsp","d2103.tsp","u2152.tsp","u2319.tsp","pr2392.tsp","pcb3038.tsp","fnl4461.tsp",
files = ["brd14051.tsp","d15112.tsp","d18512.tsp"]


def main():
	pool = multiprocessing.Pool(1)
	pool.map(run_deletion, files)

	print "Done all problems."

def run_deletion(filename):
	print
	print "Running " + filename
	subprocess.call(["../edge_deletion", "../data/" + filename])
	print "Done."
	print

	return None

if __name__ == "__main__":
	main()