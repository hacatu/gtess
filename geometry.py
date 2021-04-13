import itertools as itt
import pygame
import cmath
from collections import defaultdict, deque
import copy
from enum import IntEnum
from bitarray import bitarray
import abc

GT_EPSILON = 1e-8
WIDTH = 800

def fst(t):
	return t[0]

def snd(t):
	return t[1]

def toScreen(z):
	return (int(z.real*WIDTH + WIDTH/2), int(WIDTH/2 - z.imag*WIDTH))

def fromScreen(p):
	return p[0]/WIDTH - 0.5 - (p[1]/WIDTH - 0.5)*1J

def complexCross(a, b):
	return a.real*b.imag - a.imag*b.real

def intersectLines(a_p, a_o, b_p, b_o):
	det = complexCross(a_o, b_o);
	if abs(det) < GT_EPSILON:
		if abs(complexCross(b_p - a_p, a_o)) < GT_EPSILON:#collinear
			return a_p
		return None
	t = complexCross(b_p - a_p, a_o)/det
	return b_p + b_o*t

def complex2Cartesian(z):
	return (z.real, z.imag)

class Drawable(abc.ABC):
	@abc.abstractmethod
	def draw(self, screen):
		raise NotImplementedError("This default abstract method should not be called")

class ComplexTransformable:
	def __imul__(self, a):
		self.transform_points(lambda z: z*a)
		return self
	
	def __iadd__(self, a):
		self.transform_points(lambda z: z + a)
		return self
	
	def __isub__(self, a):
		self.transform_points(lambda z: z - a)
		return self
	
	def __idiv__(self, a):
		self.__imul__(1/a)
		return self
	
	def __mul__(self, a):
		res = copy.deepcopy(self)
		res.__imul__(a)
		return res
	
	def __rmul__(self, a):
		res = copy.deepcopy(self)
		res.__imul__(a)
		return res
	
	def __add__(self, a):
		res = copy.deepcopy(self)
		res.__iadd__(a)
		return res
	
	def __radd__(self, a):
		res = copy.deepcopy(self)
		res.__iadd__(a)
		return res
	
	def __sub__(self, a):
		res = copy.deepcopy(self)
		res.__isub__(a)
		return res
	
	def __rsub__(self, a):
		res = copy.deepcopy(self)
		res.__isub__(a)
		return res
	
	def __div__(self, a):
		res = copy.deepcopy(self)
		res.__idiv__(a)
		return res
	
	def i_transpose(self):
		self.transform_points(lambda z: (-1J*z).conjugate())
		return self
	
	def transpose(self):
		res = copy.deepcopy(self)
		res.i_transpose()
		return res

class ComposablePointList(ComplexTransformable, Drawable):
	def __init__(self):
		self.vs = []
	
	def transform_points(self, f):
		for i, z in enumerate(self.vs):
			self.vs[i] = f(z)
		return self
	
	def add_points(self, vs):
		self.vs += vs
		return self
	
	def __imatmul__(self, other):
		self.vs = [z*a for a in other.vs for z in self.vs]
		return self
	
	def __matmul__(self, other):
		res = copy.deepcopy(self)
		res.__imatmul__(other)
		return res
	
	def extend(self, *others):
		self.vs += list(itt.chain(*(other.vs for other in others)))
		return self
	
	def draw(self, screen):
		for z in self.vs:
			pygame.draw.circle(screen, 0xFFFFFFFF, toScreen(z), 2)

class EdgeKind(IntEnum):
	VALLEY = 0
	MOUNTAIN = 1
	RAW = 2

class ComposableEdgeList(Drawable):
	def __init__(self, reference_cpl=None):
		self.edges = []
		self.cpl = reference_cpl
	
	def add_edges(self, edges):
		self.edges += edges
		return self
	
	def add_edges_MVR(self, mountains=(), valleys=(), raws=()):
		self.edges += [(i, j, kind) for edges, kind in ((mountains, EdgeKind.MOUNTAIN), (valleys, EdgeKind.VALLEY), (raws, EdgeKind.RAW)) for (i, j) in edges]
		return self
	
	def extend(self, *others):
		self.vs.extend(other.vs for other in others)
		return self
	
	def replicate(self, stride, copies):
		self.edges += [((i + stride*k)%(stride*copies), (j + stride*k)%(stride*copies), kind) for k in range(1, copies) for (i, j, kind) in self.edges]
		return self
	
	def draw(self, screen):
		if self.cpl is None:
			raise ValueError("Can't draw an C. Edge List without a reference C. Point List")
		for i, j, kind in self.edges:
			pygame.draw.line(screen, (0xFFFF00FF, 0xFF00FFFF, 0xFFFFFFFF)[kind], toScreen(self.cpl.vs[i]), toScreen(self.cpl.vs[j]))

