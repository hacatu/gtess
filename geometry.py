import itertools as itt
from collections.abc import Iterable, Sequence
from typing import Union, TypeVar
import copy
from enum import IntEnum
import abc
import numbers
import numpy as np
import pygame

GT_EPSILON = 1e-8
WIDTH = 800

T = TypeVar("T")

def fst(t: Sequence[T]) -> T:
	"""Get t[0]"""
	return t[0]

def snd(t: Sequence[T]) -> T:
	"""Get t[1]"""
	return t[1]

class Scene:
	"""Handles translating between screen coordinates and model/paper coordinates, as well as drawing to the screen"""
	def __init__(self, screen : pygame.Surface, screen_size : (int, int), view_transform : np.ndarray = None):
		if (
			not isinstance(screen, pygame.Surface) or
			not (isinstance(screen_size, Iterable) and len(screen_size := tuple(screen_size)) == 2 and all(isinstance(d, int) and d > 0 for d in screen_size)) or
			not (
				view_transform is None or
				(isinstance(view_transform, np.ndarray) and view_transform.shape == (2, 3) and view_transform.dtype == np.float64))
		):
			raise TypeError("Invalid argument type")
		self.screen = screen
		self.screen_size = screen_size
		if view_transform is None:
			a : float = min(screen_size)/2
			view_transform = np.array([[a, 0, screen_size[0]/2], [0, -a, screen_size[1]/2]])
		self.setTransform(view_transform)
		self.objects : 'list[Drawable]' = []
		self.dirty = False
	
	def setTransform(self, view_transform : np.ndarray = None) -> None :
		"""Update the transform and view transform"""
		if not (
			view_transform is None or
			(isinstance(view_transform, np.ndarray) and view_transform.shape == (2, 3) and view_transform.dtype == np.float64)
		):
			raise TypeError("Invalid argument type")
		self.view_transform = view_transform
		self.inv_transform = np.linalg.inv(np.vstack([view_transform, np.array([[0, 0, 1]])]))[:2]
		self.dirty = True

	def toScreen(self, x : np.ndarray) -> (int, int):
		"""Convert a point from model coordinates to screen coordinates"""
		if not (isinstance(x, np.ndarray) and x.shape == (2,) and x.dtype == np.float64):
			raise TypeError("Input coordinate must be a length 2 ndarray of float64")
		x, y = self.view_transform @ np.vstack([x.reshape(2,1), np.ones((1,1))])
		return (int(x), int(y))
	
	def fromScreen(self, x : (int, int)) -> np.ndarray:
		"""Convert a point from screen coordinates to model coordinates"""
		x = np.array(x).reshape(2,1)
		return self.inv_transform @ np.vstack([x, np.ones((1,1))])
	
	def addObject(self, d : 'Drawable') -> None:
		self.objects.append(d)
		self.dirty = False

	def draw(self) -> None :
		if self.dirty:
			self.dirty = False
			for d in self.objects:
				d.draw(self)

def intersectLines(a_p : np.ndarray, a_o : np.ndarray, b_p : np.ndarray, b_o : np.ndarray) -> np.ndarray:
	"""Find the intersection of the line a_p + a_o*t with the line b_p + b_o*s.
	Note that this intersection may lie outside the segments (a_p, a_p + a_o) and (b_p, b_p + b_o).
	If the lines are collinear, a_p is returned.
	If the lines do not intersect (because they are parallel but not collinear), None is returned."""
	det : float = np.cross(a_o, b_o)
	if abs(det) < GT_EPSILON:
		if abs(np.cross(b_p - a_p, a_o)) < GT_EPSILON:#collinear
			return a_p
		return None
	t : float = np.cross(b_p - a_p, a_o)/det
	return b_p + b_o*t

def reflectOverUnit(a : np.ndarray, v : np.ndarray) -> np.ndarray:
	"""Reflect a 2d vector a over a unit vector v"""
	vp = np.array([v[1], -v[0]])
	return a.dot(v)*v - a.dot(vp)*vp

class M3:
	"""Convenience wrapper for a 2x3 matrix representing an affine transformation on 2D points.
	Complex numbers (including real numbers), pairs of real numbers, and 2x3 iterables of real numbers
	are coerced to transformation matrices.  Complex numbers represent rotation and/or rigid scaling,
	and pairs of real numbers represent heterogeneous scaling (separate x/y values).
	Instances of this class may be called on numpy 2x? arrays of points to transform them."""
	def __init__(self, m : Union[numbers.Number, Iterable[numbers.Real], Iterable[Iterable[numbers.Real]]]):
		if isinstance(m, numbers.Number):
			if isinstance(m, numbers.Real):
				self.m = np.array([[m, 0, 0], [0, m, 0]])
			else:
				self.m = np.array([[m.real, -m.imag, 0], [m.imag, m.real, 0]])
			return
		if not isinstance(m, Iterable):
			raise TypeError("Cannot convert value of type \"" + type(m) + "\" to a transformation matrix")
		m = list(m)
		if len(m) != 2:
			raise TypeError("Iterable has wrong number of terms (must be 2) to convert to a transformation matrix")
		if isinstance(m[0], Iterable):
			for i in range(2):
				m[i] = list(m[i])
				if len(m[i]) < 2:
					raise TypeError("Nested iterable is too short (must have 2 or 3 terms); cannot convert to a transformation matrix")
				if len(m[i]) == 2:
					m[i].append(0)
				elif len(m[i]) > 3:
					m[i].append(0)
				if not all(isinstance(x, numbers.Real) for x in m[i]):
					raise TypeError("Nested iterable contains a value that is not a real number; cannot convert to a transformation matrix")
			self.m = np.array(m)
			return
		if not all(isinstance(x, numbers.Real) for x in m):
			raise TypeError("Iterable contains a value that is not a real number; cannot convert to a transformation matrix")
		self.m = np.array([[m[0], 0, 0], [0, m[1], 0]])
	
	def __call__(self, vs : np.ndarray) -> np.ndarray:
		if len(vs.shape) == 1:
			return self.m @ np.vstack([vs.reshape((2,1)), np.ones(shape=(1,1))])
		return self.m @ np.vstack([vs, np.ones(shape=(1,vs.shape[1]))])

M3.Transpose = M3([[0, 1], [1, 0]])
M3.Conjugate = M3([1, -1])

def ucircPoint(a : float) -> np.ndarray:
	"""Return a point on the unit circle with a given angle, aka (cos(a), sin(a))."""
	return (np.cos(a), np.sin(a))

def v2(v : Union[numbers.Number, Iterable[numbers.Real]]) -> np.ndarray:
	"""Try to coerce a value into a numpy array representing a 2D point.
	In particular, complex numbers (including real numbers) are converted to (real, imag) points,
	and iterables containing 2 real numbers are converted to numpy arrays."""
	if isinstance(v, numbers.Number):
		return np.array([v.real, v.imag])
	if not isinstance(v, Iterable):
		raise TypeError("Cannot convert value of type \"" + type(v) + "\" to a vector")
	v = list(v)
	if len(v) != 2:
		raise TypeError("Iterable has wrong number of terms (must be 2) to convert to a vector")
	if not all(isinstance(x, numbers.Real) for x in v):
		raise TypeError("Iterable contains a value that is not a real number; cannot convert to a vector")
	return np.array(v)

def lexCmp(a : np.ndarray, b : np.ndarray) -> int:
	"""Compare two numpy arrays lexicographically, ignoring their shape and
	any trailing elements of the longer one if applicable.
	Returns -1 if a<b, 0 if a == b, and 1 if a > b"""
	for x, y in zip(np.nditer(a), np.nditer(b)):
		if x != y:
			return -1 if x < y else 1
	return 0

class Drawable(abc.ABC):
	"""Interface for 2D objects which can be drawn to the screen."""
	@abc.abstractmethod
	def draw(self, scene : Scene) -> None:
		"""Abstract method to draw an object to a scene"""
		raise NotImplementedError("This default abstract method should not be called")

class Steppable(abc.ABC):
	"""Quick and dirty interface for creating steppable geometry for debugging purposes"""
	@abc.abstractmethod
	def step(self):
		"""Abstract method to advance an object to its next state"""
		raise NotImplementedError("This default abstract method should not be called")

class M3Transformable(abc.ABC):
	"""Interface for 2D objects which should be able to be transformed by a "2D" matrix.
	Convenience mixin methods are provided to automatically convert multiplication by a scalar etc into a matrix transformation."""
	@abc.abstractmethod
	def apply_M3(self, m : np.ndarray) -> 'M3Transformable':
		"""Abstract method to apply a 2x3 transformation matrix to an object.
		This method must be overriden by M3Transformable classes.
		Note that transformation matrices are 2x3, with the third column containing coefficients for an implicit 1 coordinate.
		This allows any affine transformation to be encoded rather than just linear transformations.  Look up homogeneous coordinates for more detail."""
		raise NotImplementedError("This default abstract method should not be called")

	def __imul__(self, a : Union[numbers.Number, Iterable[numbers.Real]]) -> 'M3Transformable':
		"""If a is a scalar, multiply in place by it (note complex scalars have a rotation component).
		If a is a pair of real numbers, scale the x coordinates by the first and the y coordinates by the second in place."""
		if isinstance(a, numbers.Complex):
			a = np.array([[a.real, -a.imag, 0], [a.imag, a.real, 0]])
		elif isinstance(a, Iterable) and len(_a := list(a)) == 2 and isinstance(_a[0], numbers.Real):
			a = np.array([[_a[0], 0, 0], [0, _a[1], 0]])
		else:
			raise TypeError("Can't infer transformation matrix for multiplying by value of type \"" + type(a) + "\"")
		self.apply_M3(a)
		return self
	
	def __iadd__(self, a : Union[numbers.Number, Iterable[numbers.Real]]) -> 'M3Transformable':
		"""If a is a scalar, translate all points by (a.real, a.imag) in place (simplifies to x translation if a is real).
		If a is a pair of real numbers, translate all points by this pair in place."""
		if isinstance(a, numbers.Number):
			a = np.array([[1, 0, a.real], [0, 1, a.imag]])
		elif isinstance(a, Iterable) and len(_a := list(a)) == 2 and isinstance(_a[0], numbers.Real):
			a = np.array([[1, 0, _a[0], [0, 1, _a[1]]]])
		else:
			raise TypeError("Can't infer transformation matrix for adding a value of type \"" + type(a) + "\"")
		self.apply_M3(a)
		return self
	
	def __isub__(self, a : Union[numbers.Number, Iterable[numbers.Real]]) -> 'M3Transformable':
		"""Like M3Transformable.__iadd__ but reverses the translation."""
		if isinstance(a, numbers.Number):
			a = np.array([[1, 0, -a.real], [0, 1, -a.imag]])
		elif isinstance(a, Iterable) and len(_a := list(a)) == 2 and isinstance(_a[0], numbers.Real):
			a = np.array([[1, 0, -_a[0], [0, 1, -_a[1]]]])
		else:
			raise TypeError("Can't infer transformation matrix for subtracting a value of type \"" + type(a) + "\"")
		self.apply_M3(a)
		return self
	
	def __idiv__(self, a : Union[numbers.Number, Iterable[numbers.Real]]) -> 'M3Transformable':
		"""Like M3Transformable.__imul__ but performs the inverse multiplication.  Note that dividing by 0 will cause an exception."""
		if isinstance(a, numbers.Complex):
			a = 1/a
			a = np.array([[a.real, -a.imag, 0], [a.imag, a.real, 0]])
		elif isinstance(a, Iterable) and len(_a := list(a)) == 2 and isinstance(_a[0], numbers.Real):
			a = np.array([[1/_a[0], 0, 0], [0, 1/_a[1], 0]])
		else:
			raise TypeError("Can't infer transformation matrix for multiplying by value of type \"" + type(a) + "\"")
		self.apply_M3(a)
		return self
	
	def __mul__(self, a : Union[numbers.Number, Iterable[numbers.Real]]) -> 'M3Transformable':
		"""Like M3Transformable.__imul__ but creates a copy."""
		res = copy.deepcopy(self)
		res.__imul__(a)
		return res
	
	def __rmul__(self, a : Union[numbers.Number, Iterable[numbers.Real]]) -> 'M3Transformable':
		"""Like M3Transformable.__imul__ but creates a copy and has reversed arguments."""
		res = copy.deepcopy(self)
		res.__imul__(a)
		return res
	
	def __add__(self, a : Union[numbers.Number, Iterable[numbers.Real]]) -> 'M3Transformable':
		"""Like M3Transformable.__iadd__ but creates a copy."""
		res = copy.deepcopy(self)
		res.__iadd__(a)
		return res
	
	def __radd__(self, a : Union[numbers.Number, Iterable[numbers.Real]]) -> 'M3Transformable':
		"""Like M3Transformable.__iadd__ but creates a copy and has reversed arguments."""
		res = copy.deepcopy(self)
		res.__iadd__(a)
		return res
	
	def __sub__(self, a : Union[numbers.Number, Iterable[numbers.Real]]) -> 'M3Transformable':
		"""Like M3Transformable.__iadd__ but creates a copy and reverses the translation."""
		res = copy.deepcopy(self)
		res.__isub__(a)
		return res
	
	def __rsub__(self, a : Union[numbers.Number, Iterable[numbers.Real]]) -> 'M3Transformable':
		"""Like M3Transformable.__iadd__ but creates a copy, reverses the translation, and has reversed arguments."""
		res = copy.deepcopy(self)
		res.__isub__(a)
		return res
	
	def __div__(self, a : Union[numbers.Number, Iterable[numbers.Real]]) -> 'M3Transformable':
		"""Like M3Transformable.__imul__ but creates a copy and performs the inverse multiplication."""
		res = copy.deepcopy(self)
		res.__idiv__(a)
		return res
	
	def i_transpose(self) -> 'M3Transformable':
		"""Swaps the x and y coordinates of all points in the object in place.
		Equivalent to reflecting over the line y=x or applying the matrix [[0,1,0],[1,0,0]]."""
		return self.apply_M3(np.array[[0, 1, 0], [1, 0, 0]])
	
	def transpose(self) -> 'M3Transformable':
		"""Like M3Transformable.i_transpose but creates a copy."""
		res = copy.deepcopy(self)
		res.i_transpose()
		return res

class ComposablePointList(M3Transformable, Drawable):
	"""A list of 2D points that can be concatenated, transormed by a 2x3 matrix, or replicated by a set of polar or rectangular offsets."""
	def __init__(self):
		"""Create a new empy ComposablePointList."""
		self.vs = np.ndarray((2, 0))
	
	def apply_M3(self, m : np.ndarray) -> 'ComposablePointList':
		if isinstance(m, np.ndarray) and m.shape == (2, 3):
			self.vs = m @ np.vstack([self.vs, np.ones((1, self.vs.shape[1]))])
		return self

	def add_points(self, vs : Union[Iterable[Iterable[numbers.Real]], Iterable[np.ndarray]]) -> 'ComposablePointList':
		"""Add points to the list in place and return the list."""
		if not isinstance(vs, np.ndarray):
			vs1 : list[list[numbers.Real]] = [[None]*len(vs), [None]*len(vs)]
			for i, p in enumerate(vs):
				if isinstance(p, numbers.Complex):
					vs1[0][i] = p.real
					vs1[1][i] = p.imag
				elif isinstance(p, Iterable):
					if isinstance(p, np.ndarray) and len(p.shape) != 1:
						if len(p.shape) != 2:
							raise TypeError("Point list should be a list of pairs of real numbers and/or scalar complex numbers")
						if p.shape[0] < p.shape[1]:
							p = p.T
						if p.shape != (2,1):
							raise TypeError("Point list should be a list of pairs of real numbers and/or scalar complex numbers")
						vs1[0][i] = p[0,0]
						vs1[1][i] = p[1,0]
					else:
						p = list(p)
						if len(p) != 2 or not all(isinstance(x, numbers.Real) for x in p):
							raise TypeError("Point list should be a list of pairs of real numbers and/or scalar complex numbers")
						vs1[0][i] = p[0]
						vs1[1][i] = p[1]
				else:
					raise TypeError("Point list should be a list of pairs of real numbers and/or scalar complex numbers")
			vs = np.array(vs1)
		self.vs = np.hstack([self.vs, vs])
		return self
	
	def i_replicate_polar(self, other : 'ComposablePointList') -> 'ComposablePointList':
		"""Replicate the list and transform each replica by multiplying by a different point in other.
		The points in other are treated as complex scalars, with the x coordinate being the real part and the y coordinate the imaginary part.
		By specifying the nth roots of unity, n evenly spaced rotated copies of the list may be created."""
		self.vs = np.hstack([(self*complex(a, b)).vs for (a, b) in zip(*other.vs)])
		return self
	
	def replicate_polar(self, other : 'ComposablePointList') -> 'ComposablePointList':
		"""Like ComposablePointList.i_replicate_polar but creates a copy."""
		res = copy.deepcopy(self)
		res.i_replicate_polar(other)
		return res
	
	def i_replicate_rect(self, other : 'ComposablePointList') -> 'ComposablePointList':
		"""Replicate the list and transform each replica by adding a different point in other.
		By specifying an evenly spaced n by m grid of points, an n by m grid of evenly spaced copies of the list may be created."""
		self.vs = np.hstack([(self + a).vs for a in other.vs])
		return self
	
	def replicate_rect(self, other : 'ComposablePointList') -> 'ComposablePointList':
		"""Like ComposablePointList.i_replicate_rect but creates a copy."""
		res = copy.deepcopy(self)
		res.i_replicate_rect(other)
		return res
	
	def extend(self, *others: 'tuple[ComposablePointList]') -> 'ComposablePointList':
		"""Add all the points from each of the other given lists in the order given.
		Modifies this list in place and returns it."""
		self.vs = np.hstack(itt.chain(self.vs, *(other.vs for other in others)))
		return self
	
	def draw(self, scene : Scene) -> None:
		for z in self.vs.T:
			pygame.draw.circle(scene.screen, 0xFFFFFFFF, scene.toScreen(z), 2)

class EdgeKind(IntEnum):
	VALLEY = 0
	MOUNTAIN = 1
	RAW = 2

class ComposableEdgeList(Drawable):
	"""A list of crease pattern edges that can be concatenated or replicated.
	Edges are stored as (i, j, kind) triples where i is the index of the start vertex in some ComposablePointList, j is the index of the end vertex, and kind is an EdgeKind.
	The edge is not directed, but i should typically be less than j for consistency.
	Note that crease pattern edges can have even more properties such as length and dihedral angle, but these are beyond the scope of this class."""
	def __init__(self, reference_cpl : ComposablePointList = None):
		"""Create a new empty ComposableEdgeList.
		Some operations, in particular rendering, require a reference ComposablePointList be specified, but many do not."""
		self.edges : list[(int, int, EdgeKind)] = []
		self.cpl = reference_cpl
	
	def add_edges(self, edges : list[(int, int, EdgeKind)]) -> 'ComposableEdgeList':
		"""Add a list of edges given as (i, j, kind) tuples to the list in place and return the list."""
		self.edges += edges
		return self
	
	def add_edges_MVR(self, mountains : Iterable[(int, int)] = (), valleys : Iterable[(int, int)] = (), raws : Iterable[(int, int)] = ()):
		"""Add new mountain, valley, and raw edges to the list in place and return the list.
		This method accepts the edges to add as three iterables of (i, j) pairs for the mountain edges, valley edges, and raw edges to add respectively.
		Any or all of these iterables may be omitted if no edges of that type need to be added."""
		self.edges += [(i, j, kind) for edges, kind in ((mountains, EdgeKind.MOUNTAIN), (valleys, EdgeKind.VALLEY), (raws, EdgeKind.RAW)) for (i, j) in edges]
		return self
	
	def extend(self, *others: 'tuple[ComposableEdgeList]') -> 'ComposableEdgeList':
		"""Add all the edges from each of the other given lists in the order given.
		Modifies this list in place and returns it."""
		self.edges.extend(other.edges for other in others)
		return self
	
	def i_replicate(self, stride: int, copies: int) -> 'ComposableEdgeList':
		"""Replicate the list a number of times equal to copies.
		This method should be used when the crease pattern has inherent polar or rectangular symmetery.
		The underlying list of points should have been created via ComposablePointList.i_replicate_polar or similar functions, and
		stride is the number of points in a single unit of the symmetric pattern.
		stride is simply added to all indices in the given edge list to generate each replica edge list."""
		self.edges += [((i + stride*k)%(stride*copies), (j + stride*k)%(stride*copies), kind) for k in range(1, copies) for (i, j, kind) in self.edges]
		return self
	
	def draw(self, scene):
		if self.cpl is None:
			raise ValueError("Can't draw a C. Edge List without a reference C. Point List")
		for i, j, kind in self.edges:
			color = (0xFFFF00FF, 0xFF00FFFF, 0xFFFFFFFF)[kind]
			pygame.draw.line(scene.screen, color, scene.toScreen(self.cpl.vs[i]), scene.toScreen(self.cpl.vs[j]))

