from tsp_loader import Graph
from graphics import *
import sys

def main():
	graph = Graph("../data/pr20.tsp")

	win = GraphWin(width = 600, height = 600) # create a window

	width = max(graph.maxP.x - graph.minP.x, graph.maxP.y - graph.minP.y)
	offset = width / 20.0
	win.setCoords(graph.minP.x - offset, graph.minP.y - offset, graph.minP.x + width + offset, graph.minP.y + width + offset) # set the coordinates of the window; bottom left is (0, 0) and top right is (10, 10)

	#Draw points
	for i,point in enumerate(graph.points):
		offset = width/200
		pt = Point(point.x, point.y)
		cir = Circle(pt,offset)
		cir.setFill('blue')
		cir.draw(win)


		message = Text(Point(point.x + offset* 3, point.y + offset* 3), str(i))
		message.draw(win)


	# for edge in graph.edges:
	# 	points = [graph.points[i] for i in edge.end]
	# 	pts = [Point(p.x,p.y) for p in points]

	# 	line = Line(pts[0], pts[1])
	# 	line.draw(win)

	win.getMouse() # pause before closing



if __name__ == "__main__":
	main()