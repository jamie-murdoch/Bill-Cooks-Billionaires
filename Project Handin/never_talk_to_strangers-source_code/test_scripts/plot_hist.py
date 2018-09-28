import matplotlib as m
import json
import os
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab


#files = ["pr20.tsp","a80.tsp", "a180.tsp","a280.tsp","pr1002.tsp","u1060.tsp","vm1084.tsp","pcb1173.tsp","d1291.tsp","rl1304.tsp","rl1323.tsp","nrw1379.tsp","u1432.tsp","d1655.tsp","vm1748.tsp","u1817.tsp","rl1889.tsp","d2103.tsp","u2152.tsp","u2319.tsp","pr2392.tsp","pcb3038.tsp","fnl4461.tsp","brd14051.tsp","d15112.tsp","d18512.tsp"]
files = ["pr1002.tsp","u1060.tsp","vm1084.tsp","pcb1173.tsp","d1291.tsp","rl1304.tsp","rl1323.tsp","nrw1379.tsp","u1432.tsp","d1655.tsp","vm1748.tsp","u1817.tsp","rl1889.tsp","d2103.tsp","u2152.tsp","u2319.tsp","pr2392.tsp","pcb3038.tsp","fnl4461.tsp","brd14051.tsp","d15112.tsp","d18512.tsp"]

#Note: need to run this from test/results with nothing else present
def main():
	time_data = load_time_data()

	speedups = []
	probs = []
	for prob_name in files:
		if prob_name in time_data:
			speedups.append(compute_speedup(time_data[prob_name]))
			probs.append(prob_name)

	print zip(probs, speedups)
	print len(speedups)

	num_bins = 5
	# the histogram of the data
	n, bins, patches = plt.hist(speedups, num_bins, normed=True, facecolor='green', alpha=0.5)

	# add a 'best fit' line
	#y = mlab.normpdf(bins, mu, sigma)
	#plt.plot(bins, y, 'r--')
	plt.xlabel('Speedup')
	plt.ylabel('Frequency')
	plt.title(r'Histogram of TSP Solve Speedup')

	# Tweak spacing to prevent clipping of ylabel
	#plt.subplots_adjust(left=0.15)
	plt.show()


def load_time_data():
	time_data = {}

	for prob_name in files:
		time_fname = prob_name + "-times.txt"
		time_path = os.path.join("results/", prob_name.split(".")[0], time_fname)

		if os.path.isfile(time_path):
			with open(time_path, 'r') as f:
				data = json.load(f)
				for key in data:
					time_data[key] = data[key]

	return time_data

def compute_speedup(data):
	t_old = data["orig"]
	t_new = data["pruned"]

	return t_old/t_new

if __name__ == "__main__":
	main()