import itertools as itt
import pygame
import cmath
from collections import defaultdict, deque
import copy
from enum import IntEnum
from bitarray import bitarray
import abc

from geometry import *

class CreasePatternRegion(ComplexTransformable):
	def __init__(self, vs, is_exterior=False):
		self.vs = vs
		self.edges = [None]*len(vs)
		self.is_exterior = is_exterior
	
	def transform_points(self, f):
		for i, z in enumerate(self.vs):
			self.vs[i] = f(z)
		return self

class CreasePatternEdge(Drawable):
	def __init__(self, a, ai, b, bi, kind):
		self.a = a
		self.ai = ai
		self.b = b
		self.bi = bi
		self.kind = kind
		a.edges[ai] = self
		b.edges[bi] = self
	
	def draw(self, screen):
		pygame.draw.line(screen, (0xFFFF00FF, 0xFF00FFFF, 0xFFFFFFFF)[self.kind], toScreen(self.a.vs[self.ai]), toScreen(self.a.vs[(self.ai+1)%len(self.a.vs)]))

class PartialCreasePattern(ComplexTransformable, Drawable):
	def __init__(self):
		self.regions = []
		self.edges = []
	
	def transform_points(self, f):
		for region in self.regions:
			region.transform_points(f)
		return self
	
	def draw(self, screen):
		for edge in self.edges:
			edge.draw(screen)

class PlaneGraphNode:
	def __init__(self, graph, point, idx):
		self.graph = graph
		self.point = point
		self.idx = idx
		self.neighbors = []
		self.is_sorted = False
		self.num_unvisited_neighbors = 0
	
	def add_neighbor(self, j):
		self.neighbors.append(j)
		self.num_unvisited_neighbors += 1
		return self
	
	def _sort_neighbors_ccw_from(self, f):
		cut_phase = 0 if f == -1 else cmath.phase(self.graph.vertices[f].point - self.point)
		def key(j):
			return (cmath.phase(self.graph.vertices[j].point - self.point) - cut_phase)%(2*cmath.pi)
		self.neighbors.sort(key=key)
		self.is_sorted = True
	
	def has_unvisited_neighbors(self):
		return self.num_unvisited_neighbors != 0
	
	def visit_from(self, f, is_new_face):
		self.num_unvisited_neighbors -= 1
		if self.num_unvisited_neighbors == 0:
			self.graph.done_vertices[self.idx] = True
		if not self.is_sorted:
			self._sort_neighbors_ccw_from(f)
		if f == -1:
			for j in self.neighbors:
				if (self.idx, j) not in self.graph.visited_edges:
					self.graph.visited_edges.add((self.idx, j))
					return j
		j = self.neighbors.index(f) - 1
		while is_new_face and (self.idx, self.neighbors[j]) in self.graph.visited_edges:
			j -= 1
		self.graph.visited_edges.add((self.idx, self.neighbors[j]))
		return self.neighbors[j]

class PartialCreasePatternBuilder:
	def __init__(self):
		self.vertices = []
		self.edges = {}
	
	def with_vertices(self, vs):
		self.vertices = [PlaneGraphNode(self, v, i) for (i, v) in enumerate(vs)]
		return self
	
	def with_edges(self, edges):
		for i, j, kind in edges:
			self.vertices[i].add_neighbor(j)
			self.vertices[j].add_neighbor(i)
			if complex2Cartesian(self.vertices[j].point) < complex2Cartesian(self.vertices[i].point):
				i, j = j, i
			self.edges[(i, j)] = (None, -1, kind)
		return self
	
	def _add_res_face(self, face):
		interior = CreasePatternRegion([self.vertices[i].point for i in face])
		self.res.regions.append(interior)
		face.append(face[0])
		for e in range(len(face) - 1):
			i = face[e]
			j = face[e+1]
			if complex2Cartesian(self.vertices[j].point) < complex2Cartesian(self.vertices[i].point):
				i, j = j, i
			back, k, kind = self.edges[(i, j)]
			if back is not None:
				self.res.edges.append(CreasePatternEdge(interior, e, back, k, kind))
			else:
				self.edges[(i, j)] = (interior, e, kind)
	
	def _check_result(self):
		edge_occs = set()
		z2i = {node.point: idx for (idx, node) in enumerate(self.vertices)}
		#ensure each edge occurs at most once in a given direction
		#this ensures each edge occurs exactly once in both directions after later checks
		for region in self.res.regions:
			for i in range(len(region.vs)):
				edge = (z2i[region.vs[i]], z2i[region.vs[i-1]])
				if edge in edge_occs:
					print("!!! edge occurs multiple times")
					return False
				edge_occs.add(edge)
		#ensure each edge in an output face corresponds to some edge in the input list of edges
		for i, j in edge_occs:
			if not ((i, j) in self.edges or (j, i) in self.edges):
				print("!!! output edge is not an input edge")
				return False
		#ensure there are twice as many edges in the output faces as the input list of edges,
		#and since we now know each occurs at most once in a given direction and corresponds
		#to some input edge, we know each input edge occurs exactly once in both directions
		#in the output
		status = len(edge_occs) == 2*len(self.edges)
		if not status:
			print("!!! not all input edges exist in output")
		return status
	
	def build(self):
		#print("Building PCP")
		#print("points: ", [v.point for v in self.vertices])
		#print("starting from point 0")
		self.res = PartialCreasePattern()
		self.done_vertices = bitarray(len(self.vertices))
		self.done_vertices.setall(False)
		self.visited_edges = set()
		face = []
		i1 = 0
		v1 = self.vertices[i1]
		f = -1
		while True:
			#print(" vertex: ", v1.point)
			i2 = v1.visit_from(f, not face)
			face.append(i1)
			#print(" face: ", face)
			v2 = self.vertices[i2]
			if i2 == face[0]:
				#print(" face completed")
				self._add_res_face(face)
				if v2.has_unvisited_neighbors():
					f = face[1]
					face = []
					v1 = v2
					i1 = i2
				else:
					try:
						i1 = self.done_vertices.index(False)
					except ValueError:
						if not self._check_result():
							print("!!! PartialCreasePatternBuilder._check_result failed")
						return self.res
					f = -1
					face = []
					v1 = self.vertices[i1]
			else:
				f = i1
				v1 = v2
				i1 = i2
		#TODO: identify exterior region (and maybe do convexity check in general)
		if not self._check_result():
			print("!!! PartialCreasePatternBuilder._check_result failed")
		return self.res

