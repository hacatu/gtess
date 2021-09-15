import itertools as itt
import pygame
from collections import defaultdict, deque
import copy
from enum import IntEnum
from bitarray import bitarray
import abc
from scipy.spatial import KDTree

from geometry import *

class CreasePatternRegion(M3Transformable):
	def __init__(self, vs : np.ndarray, is_exterior : bool = False):
		if len(vs.shape) != 2 or vs.shape[0] != 2:
			raise TypeError("vs must be a 2 by x array")
		self.vs = vs
		self.edges : 'list[CreasePatternEdge]' = [None]*vs.shape[1]
		self.is_exterior = is_exterior
	
	def apply_M3(self, m : np.ndarray) -> 'CreasePatternRegion':
		"""Implements M3Transformalbe.apply_M3."""
		if not isinstance(m, np.ndarray) or m.shape != (2, 3):
			raise TypeError("m must be a 2 by 3 array")
		self.vs = m @ np.vstack([self.vs, np.ones((1, self.vs.shape[1]))])
		return self

class CreasePatternEdge(Drawable):
	def __init__(self, a : CreasePatternRegion, ai : int, b : CreasePatternRegion, bi : int, kind : EdgeKind):
		self.a = a
		self.ai = ai
		self.b = b
		self.bi = bi
		self.kind = kind
		a.edges[ai] = self
		b.edges[bi] = self
	
	def draw(self, scene : Scene) -> None:
		if self.kind != EdgeKind.MOUNTAIN:
			pass
		color = (0xFFFF00FF, 0xFF00FFFF, 0xFFFFFFFF)[self.kind]
		#pygame.draw.line(scene.screen, color, scene.toScreen(self.a.vs[:,self.ai]), scene.toScreen(self.b.vs[:,self.bi]))
		pygame.draw.line(scene.screen, color, scene.toScreen(self.a.vs[:,self.ai]), scene.toScreen(self.a.vs[:,(self.ai+1)%self.a.vs.shape[1]]))

class PartialCreasePattern(M3Transformable, Drawable):
	def __init__(self):
		self.regions : list[CreasePatternRegion] = []
		self.edges : list[CreasePatternEdge] = []
	
	def apply_M3(self, m : np.ndarray) -> 'PartialCreasePattern':
		"""Implements M3Transformalbe.apply_M3."""
		for region in self.regions:
			region.apply_M3(m)
		return self
	
	def draw(self, scene : Scene) -> None:
		for edge in self.edges:
			edge.draw(scene)

class PlaneGraphNode:
	def __init__(self, graph : 'PartialCreasePatternBuilder', point : np.ndarray, idx : int):
		if not isinstance(point, np.ndarray) or point.shape != (2,):
			raise TypeError("Invalid argument: point must be an array with 2 elements")
		self.graph = graph
		self.point = point
		self.idx = idx
		self.neighbors : list[int] = []
		self.is_sorted = False
		self.num_unvisited_neighbors : int = 0
	
	def add_neighbor(self, j : int) -> 'PlanarGraphNode':
		self.neighbors.append(j)
		self.num_unvisited_neighbors += 1
		return self
	
	def _sort_neighbors_ccw_from(self, f : int) -> None:
		if f == -1:
			cut_phase = 0.
		else:
			x, y = self.graph.vertices[f].point - self.point
			cut_phase = np.arctan2(x, y)
		def key(j : int) -> float:
			x, y = self.graph.vertices[j].point - self.point
			return (np.arctan2(x, y) - cut_phase)%(2*np.pi)
		self.neighbors.sort(key=key)
		self.is_sorted = True
	
	def has_unvisited_neighbors(self) -> bool:
		return self.num_unvisited_neighbors != 0
	
	def visit_from(self, f : int, is_new_face : bool) -> int:
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

class PartialCreasePatternBuilder(M3Transformable, Drawable, Steppable):
	def __init__(self):
		self.vertices : list[PlaneGraphNode] = []
		self.edges : dict[(int, int), (CreasePatternRegion, CreasePatternRegion, EdgeKind)] = {}
		self.current_transform = np.array([[1, 0, 0], [0, 1, 0]])
		self.res = PartialCreasePattern()
		self.done_building = False
	
	def apply_M3(self, m : np.ndarray) -> 'PartialCreasePatternBuilder':
		self.current_transform = m @ np.vstack([self.current_transform, np.array([[0, 0, 1]])])
	
	def draw(self, scene : Scene) -> None:
		for edge in self.res.edges:
			color = (0xFFFF00FF, 0xFF00FFFF, 0xFFFFFFFF)[edge.kind]
			start = scene.toScreen(edge.a.vs[edge.ai].apply_M3(self.current_transform))
			end = scene.toScreen(edge.a.vs[(edge.ai+1)%len(edge.a.vs)].apply_M3(self.current_transform))
			pygame.draw.line(scene.screen, color, start, end)

	def with_vertices(self, vs : np.ndarray) -> 'PartialCreasePatternBuilder':
		self.vertices = [PlaneGraphNode(self, v, i) for (i, v) in enumerate(vs.T)]
		return self
	
	def with_edges(self, edges : list[CreasePatternEdge]) -> 'PartialCreasePatternBuilder':
		for i, j, kind in edges:
			self.vertices[i].add_neighbor(j)
			self.vertices[j].add_neighbor(i)
			if lexCmp(self.vertices[j].point, self.vertices[i].point) < 0:
				i, j = j, i
			self.edges[(i, j)] = (None, -1, kind)
		return self
	
	def start_build(self) -> None:
		#print("Building PCP")
		#print("points: ", [v.point for v in self.vertices])
		#print("starting from point 0")
		self.done_vertices = bitarray(len(self.vertices))
		self.done_vertices.setall(False)
		self.visited_edges : set[(int, int)] = set()
		self.face : list[int] = []
		self.i1 : int = 0
		self.v1 : PlaneGraphNode = self.vertices[self.i1]
		self.f : int = -1

	def step(self) -> None:
		if self.done_building:
			print("No more steps!  Building is done!")
			return
		print(" vertex: ", self.v1.point)
		self.i2 = self.v1.visit_from(self.f, not self.face)
		self.face.append(self.i1)
		#print(" face: ", face)
		self.v2 = self.vertices[self.i2]
		if self.i2 == self.face[0]:
			print(" face completed")
			self._add_res_face(self.face)
			if self.v2.has_unvisited_neighbors():
				self.f = self.face[1]
				self.face = []
				self.v1 = self.v2
				self.i1 = self.i2
			else:
				try:
					self.i1 = self.done_vertices.index(False)
				except ValueError:
					self.done_building = True
					if not self._check_result():
						print("!!! PartialCreasePatternBuilder._check_result failed")
					return
				self.f = -1
				self.face = []
				self.v1 = self.vertices[self.i1]
		else:
			self.f = self.i1
			self.v1 = self.v2
			self.i1 = self.i2

	def _add_res_face(self, face : list[int]) -> None:
		interior = CreasePatternRegion(np.hstack([self.vertices[i].point.reshape((2,1)) for i in face]))
		self.res.regions.append(interior)
		face.append(face[0])
		for e in range(len(face) - 1):
			i = face[e]
			j = face[e+1]
			if lexCmp(self.vertices[j].point, self.vertices[i].point) < 0:
				i, j = j, i
			back, k, kind = self.edges[(i, j)]
			if back is not None:
				self.res.edges.append(CreasePatternEdge(interior, e, back, k, kind))
			else:
				self.edges[(i, j)] = (interior, e, kind)
	
	def _check_result(self) -> bool:
		edge_occs = set()
		points = np.vstack([node.point.reshape((1,2)) for node in self.vertices])
		print(points.shape)
		z2i = KDTree(points)
		#ensure each edge occurs at most once in a given direction
		#this ensures each edge occurs exactly once in both directions after later checks
		for region in self.res.regions:
			for i in range(len(region.vs)):
				_, edge = z2i.query([region.vs[:,i], region.vs[:,i-1]])
				edge = tuple(edge)
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
	
	def build(self) -> PartialCreasePattern:
		self.start_build()
		while not self.done_building:
			self.step()
		return self.res

