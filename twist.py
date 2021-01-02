import numpy as np
import itertools as itt
import pygame
import cmath
import sys
import time
from collections import defaultdict
import operator
import functools
import copy
from enum import IntEnum

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

pygame.init()
screen = pygame.display.set_mode((WIDTH, WIDTH))
pygame.display.set_caption("Gtess")

class Drawable:
	pass

class ComplexTransformable:
	def __imul__(self, a):
		self.transform_points(lambda z: z*a)
	
	def __iadd__(self, a):
		self.transform_points(lambda z: z + a)
	
	def __isub__(self, a):
		self.transform_points(lambda z: z - a)
	
	def __idiv__(self, a):
		self.__imul__(1/a)
	
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

class CreasePatternRegion(ComplexTransformable):
	def __init__(self, vs, is_exterior=False):
		self.vs = vs
		self.edges = [None]*len(vs)
		self.is_exterior = is_exterior
	
	def transform_points(self, f):
		for i, z in enumerate(self.vs):
			self.vs[i] = f(z)

class EdgeKind(IntEnum):
	VALLEY = 0
	MOUNTAIN = 1
	RAW = 2

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
	
	def draw(self, screen):
		for edge in self.edges:
			edge.draw(screen)

class FlowerTowerBuilder:
	def __init__(self, n):
		self.n = n
		self.a = 1 + 0J
		self.z = 0J
	
	def scale(self, a):
		self.a *= a
		self.z *= a
		return self
	
	def translate(self, z):
		self.z += z
		return self
	
	def build(self):
		n = self.n
		twist_vs = [cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		inner_vs = [0.5/cmath.cos(cmath.pi/n)*cmath.exp(2*cmath.pi*(i/n + 0.5/n)*1J) for i in range(n)]
		m_inner_vs = [(2*cmath.cos(cmath.pi/n) - 0.5/cmath.cos(cmath.pi/n))*cmath.exp(2*cmath.pi*(i/n + 0.5/n)*1J) for i in range(n)]
		outer_vs = [2*cmath.cos(cmath.pi/n)*cmath.exp(2*cmath.pi*(i/n + 0.5/n)*1J) for i in range(n)]
		twist_cs = [cmath.cos(cmath.pi/n)*cmath.exp(2*cmath.pi*(i/n + 0.5/n)*1J) for i in range(n)]
		outer_cs = [2*cmath.cos(cmath.pi/n)**2*cmath.exp(2*cmath.pi*(i/n)*1J) for i in range(n)]
		
		x_outer_ps = intersectLines(m_inner_vs[0], outer_cs[0] - m_inner_vs[0], twist_vs[0], outer_vs[0] - twist_vs[0])
		x_outer_ps = [x_outer_ps*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		
		lpleat_vs = outer_cs[0] + abs(outer_vs[0] - outer_cs[0])/2**.5*cmath.exp(cmath.pi*0.25J)
		rpleat_vs = lpleat_vs.conjugate()
		lpleat_vs = [lpleat_vs*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		rpleat_vs = [rpleat_vs*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		
		#stella_ps = outer_cs[0] + abs(outer_vs[0] - outer_vs[-1])
		stella_ps = outer_cs[0] + abs(lpleat_vs[0] - rpleat_vs[0])
		stella_ls = stella_ps + abs(inner_vs[0] - inner_vs[-1])/2**.5*cmath.exp(cmath.pi*0.75J)
		ext_lpleat_vs = complex(stella_ls.real, lpleat_vs[0].imag)
		ext_lpleat_cs = complex(stella_ls.real, outer_vs[0].imag)
		stella_rs = stella_ls.conjugate()
		ext_rpleat_vs = ext_lpleat_vs.conjugate()
		ext_rpleat_cs = ext_lpleat_cs.conjugate()
		stella_ps = [stella_ps*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		stella_ls = [stella_ls*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		stella_rs = [stella_rs*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		ext_lpleat_vs = [ext_lpleat_vs*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		ext_rpleat_vs = [ext_rpleat_vs*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		ext_lpleat_cs = [ext_lpleat_cs*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		ext_rpleat_cs = [ext_rpleat_cs*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		
		ext_vs = intersectLines(stella_ls[0], stella_ls[0] - stella_rs[0], stella_rs[1], stella_rs[1] - stella_ls[1])
		ext_vs = [ext_vs*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		
		inner_region = CreasePatternRegion(inner_vs)
		inter_quads = [CreasePatternRegion([inner_vs[i-1], twist_vs[i], twist_cs[i], inner_vs[i]]) for i in range(n)]
		inter_tris = [CreasePatternRegion([inner_vs[i-1], twist_cs[i-1], twist_vs[i]]) for i in range(n)]
		ext1_quads = [CreasePatternRegion([twist_vs[i], x_outer_ps[i], m_inner_vs[i], twist_cs[i]]) for i in range(n)]
		ext1_tris = [CreasePatternRegion([twist_vs[i], twist_cs[i-1], m_inner_vs[i-1]]) for i in range(n)]
		ext2_tris = [CreasePatternRegion([twist_vs[i], m_inner_vs[i-1], outer_cs[i]]) for i in range(n)]
		ext3_tris = [CreasePatternRegion([twist_vs[i], outer_cs[i], x_outer_ps[i]]) for i in range(n)]
		ext2_quads = [CreasePatternRegion([m_inner_vs[i-1], x_outer_ps[i-1], outer_vs[i-1], outer_cs[i]]) for i in range(n)]
		ext3_quads = [CreasePatternRegion([x_outer_ps[i], outer_cs[i], lpleat_vs[i], outer_vs[i]]) for i in range(n)]
		swivel_tris = [CreasePatternRegion([outer_cs[i], outer_vs[i-1], rpleat_vs[i]]) for i in range(n)]
		rabbit_tris = [CreasePatternRegion([outer_cs[i], rpleat_vs[i], lpleat_vs[i]]) for i in range(n)]
		
		stella_hepts = [CreasePatternRegion([lpleat_vs[i], rpleat_vs[i], ext_rpleat_vs[i], stella_rs[i], stella_ps[i], stella_ls[i], ext_lpleat_vs[i]]) for i in range(n)]
		swivel_traps = [CreasePatternRegion([outer_vs[i-1], ext_rpleat_cs[i], ext_rpleat_vs[i], rpleat_vs[i]]) for i in range(n)]
		rabbit_traps = [CreasePatternRegion([outer_vs[i], lpleat_vs[i], ext_lpleat_vs[i], ext_lpleat_cs[i]]) for i in range(n)]
		kites = [CreasePatternRegion([outer_vs[i], ext_lpleat_cs[i], ext_vs[i], ext_rpleat_cs[(i+1)%n]]) for i in range(n)]
		
		exterior_region = []
		for i in range(n):
			exterior_region += [ext_vs[-i], ext_lpleat_cs[-i], ext_lpleat_vs[-i], stella_ls[-i], stella_ps[-i], stella_rs[-i], ext_rpleat_vs[-i], ext_rpleat_cs[-i]]
		exterior_region = CreasePatternRegion(exterior_region, True)
		
		edges = []
		
		for i in range(n):
			edges.append(CreasePatternEdge(inner_region, i, inter_quads[i], 3, EdgeKind.VALLEY))
			edges.append(CreasePatternEdge(inter_quads[i], 0, inter_tris[i], 2, EdgeKind.MOUNTAIN))
			edges.append(CreasePatternEdge(inter_quads[i], 1, ext1_quads[i], 3, EdgeKind.MOUNTAIN))
			edges.append(CreasePatternEdge(inter_tris[i], 1, ext1_tris[i], 0, EdgeKind.MOUNTAIN))
			edges.append(CreasePatternEdge(ext1_tris[i], 2, ext2_tris[i], 0, EdgeKind.VALLEY))
			edges.append(CreasePatternEdge(ext1_quads[i], 0, ext3_tris[i], 2, EdgeKind.VALLEY))
			edges.append(CreasePatternEdge(ext2_tris[i], 2, ext3_tris[i], 0, EdgeKind.MOUNTAIN))
			edges.append(CreasePatternEdge(ext2_tris[i], 1, ext2_quads[i], 3, EdgeKind.MOUNTAIN))
			edges.append(CreasePatternEdge(ext3_tris[i], 1, ext3_quads[i], 0, EdgeKind.VALLEY))
			edges.append(CreasePatternEdge(ext2_quads[i], 2, swivel_tris[i], 0, EdgeKind.VALLEY))
			edges.append(CreasePatternEdge(ext3_quads[i], 1, rabbit_tris[i], 2, EdgeKind.MOUNTAIN))
			edges.append(CreasePatternEdge(swivel_tris[i], 2, rabbit_tris[i], 0, EdgeKind.MOUNTAIN))
			
			edges.append(CreasePatternEdge(inter_tris[i], 0, inter_quads[i-1], 2, EdgeKind.VALLEY))
			edges.append(CreasePatternEdge(ext1_tris[i], 1, ext1_quads[i-1], 2, EdgeKind.MOUNTAIN))
			edges.append(CreasePatternEdge(ext2_quads[i], 0, ext1_quads[i-1], 1, EdgeKind.MOUNTAIN))
			edges.append(CreasePatternEdge(ext2_quads[i], 1, ext3_quads[i-1], 3, EdgeKind.VALLEY))
			
			edges.append(CreasePatternEdge(rabbit_tris[i], 1, stella_hepts[i], 0, EdgeKind.VALLEY))
			edges.append(CreasePatternEdge(swivel_tris[i], 1, swivel_traps[i], 3, EdgeKind.MOUNTAIN))
			edges.append(CreasePatternEdge(ext3_quads[i], 2, rabbit_traps[i], 0, EdgeKind.MOUNTAIN))
			edges.append(CreasePatternEdge(kites[i], 0, rabbit_traps[i], 3, EdgeKind.VALLEY))
			edges.append(CreasePatternEdge(kites[i], 3, swivel_traps[i], 0, EdgeKind.VALLEY))
			edges.append(CreasePatternEdge(rabbit_traps[i], 1, stella_hepts[i], 6, EdgeKind.MOUNTAIN))
			edges.append(CreasePatternEdge(swivel_traps[i], 2, stella_hepts[i], 1, EdgeKind.MOUNTAIN))
			
			if i == 0:
				edges.append(CreasePatternEdge(kites[i], 1, exterior_region, 0, EdgeKind.RAW))
			else:
				edges.append(CreasePatternEdge(kites[i], 1, exterior_region, 8*(n - i) - 0, EdgeKind.RAW))
			edges.append(CreasePatternEdge(rabbit_traps[i], 2, exterior_region, 8*(n - i) - 1, EdgeKind.RAW))
			for k in range(2, 6):
				edges.append(CreasePatternEdge(stella_hepts[i], 7 - k, exterior_region, 8*(n - i) - k, EdgeKind.RAW))
			edges.append(CreasePatternEdge(swivel_traps[i], 1, exterior_region, 8*(n - i) - 6, EdgeKind.RAW))
			edges.append(CreasePatternEdge(kites[i-1], 2, exterior_region, 8*(n - i) - 7, EdgeKind.RAW))
		
		res = PartialCreasePattern()
		res.regions = [inner_region, *inter_quads, *inter_tris, *ext1_quads, *ext1_tris, *ext2_tris, *ext3_tris, *ext2_quads, *ext3_quads, *swivel_tris, *rabbit_tris, *stella_hepts, *swivel_traps, *rabbit_traps, *kites, exterior_region]
		res.edges = edges
		if self.a != 1:
			res.transform_points(lambda z: z*self.a)
			#res *= self.a
		if self.z != 0:
			res.transform_points(lambda z: z + self.z)
			#res += self.z
		return res

twists = [FlowerTowerBuilder(8).scale(0.06).translate(-0.32-0.32J).build()]
for i in range(3):
	twists.append(twists[-1] + (0.22+0.22J))
twists.append(twists[1] + (-0.22+0.22J))
twists.append(twists[2] + (0.22-0.22J))

while True:
	for e in pygame.event.get():
		if e.type == pygame.QUIT or (e.type == pygame.KEYDOWN and e.key == pygame.K_ESCAPE):
			pygame.quit()
			sys.exit(0)
	screen.fill(0xFF000000)
	for twist in twists:
		twist.draw(screen)
	pygame.display.flip()
	#time.sleep(100)

