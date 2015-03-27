def main():
	pass

class Point:
	def __init__(self, x, y):
		self.x = x
		self.y = y

class Edge:
	def __init__(self, e0, e1):
		self.end = [e0, e1]

class Graph:
	def __init__(self, filename):
		self.points = []
		self.edges = []
		self.maxP = Point(float("-inf"), float("-inf"))
		self.minP = Point(float("inf"), float("inf"))

		with open(filename, 'r') as tspfile:
			tspfile.readline()
			tspfile.readline()
			tspfile.readline()
			tspfile.readline()
			tspfile.readline()
			tspfile.readline()

			for line in tspfile:
				tokens = line.split()
				if tokens[0] != "EOF":
					#print tokens[1]
					self.points.append(Point(float(tokens[1]), float(tokens[2])))
					self.maxP.x = max(self.maxP.x, self.points[-1].x)
					self.maxP.y = max(self.maxP.y, self.points[-1].y)
					self.minP.x = min(self.minP.x, self.points[-1].x)
					self.minP.y = min(self.minP.y, self.points[-1].y)

		for i in xrange(0, len(self.points)):
			for j in xrange(i + 1, len(self.points)):
				self.edges.append(Edge(i,j))


if __name__ == "__main__":
	main()