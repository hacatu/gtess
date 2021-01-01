from math import hypot, sqrt, log, acos, pi

class Point:
	def __init__(self, x, y):
		self.x = x
		self.y = y

class Polygon:
	def __init__(self, vertices, vertex_properties, edge_properties):
		self.vertices = vertices
		self.vertex_properties = vertex_properties
		self.edge_properties = edge_properties

class Tiling:
	def __init__(self, absolute_tiles):
		"""
		absolute_tiles is an array of Polygons with Point vertices,
		bool vertex_properties with True if there will be a node in the twist graph
		at that vertex and False if there will not (ie it is face-internal),
		and int edge_properties corresponding to an id indicating which edges should be paired up.
		Tiles will not be automatically mirrored.
		"""
		edge_pastings = []
		edge_deficit = 0
		for i, tile in enumerate(self.absolute_tiles):
			for j, edge_id in enumerate(tile.edge_properties):
				if edge_id >= len(edge_pastings):
					new_edge_ids = edge_id - len(edge_pastings) + 1
					edge_pastings += [()]*new_edge_ids
					edge_deficit += 2*new_edge_ids
				if len(edge_pastings[edge_id]) == 2:
					pass #Oh no, that's an error
				edge_pastings[edge_id] += ((i, j),)
				edge_deficit -= 1
		if edge_deficit:
			pass #Oh no, that's an error
		self.edge_pastings = edge_pastings
		self.absolute_tiles = absolute_tiles
		self.checkVertices()
	
	def checkVertices(self):
		checked = [set() for _ in self.absolute_tiles]
		for i, tile in enumerate(self.absolute_tiles):
			if len(checked[i]) == len(tile.vertices):
				continue
			for j, edge_id in enumerate(tile.edge_properties):
				if j in checked[i]:
					continue
				checked[i].add(j)
				angle = 0
				first_edge = (i, j)
				cw_edge = first_edge
				ccw_edge = (i, (j + len(tile.vertices) - 1)%len(tile.vertices))
				i1, j1 = i, j
				tile1 = tile
				while True:
					a = tile1.vertices[(j1 + 1)%len(tile1.vertices)]
					b = tile1.vertices[j]
					c = tile1.vertices[ccw_edge[1]]
					dax = a.x - b.x
					day = a.y - b.y
					dcx = c.x - b.x
					dcy = c.y - b.y
					mda = hypot(dax, day)
					dax /= mda
					day /= mda
					mdc = hypot(dcx, dcy)
					dcx /= mdc
					dcy /= mdc
					angle += acos(dax*dcx + day*dcy)
					edge_a, edge_b = self.edge_pastings[edge_id]
					cw_edge = edge_a if edge_b == ccw_edge else edge_b
					if cw_edge == first_edge:
						break
					i1, j1 = cw_edge
					checked[i1].add(j1)
					tile1 = self.absolute_tiles[i1]
					ccw_edge = (i1, (j1 + len(tile1.vertices) - 1)%len(tile1.vertices))
				period = 2*pi/angle
				if period%1 > 0.000001:
					pass #Oh no, that's an error

