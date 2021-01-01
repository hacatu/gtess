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


class FlowerTowerUnit:
	def __init__(self, n, r=1, z=0J):
		self.twist_vs = [cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		self.inner_vs = [0.5/cmath.cos(cmath.pi/n)*cmath.exp(2*cmath.pi*(i/n + 0.5/n)*1J) for i in range(n)]
		self.m_inner_vs = [(2*cmath.cos(cmath.pi/n) - 0.5/cmath.cos(cmath.pi/n))*cmath.exp(2*cmath.pi*(i/n + 0.5/n)*1J) for i in range(n)]
		self.outer_vs = [2*cmath.cos(cmath.pi/n)*cmath.exp(2*cmath.pi*(i/n + 0.5/n)*1J) for i in range(n)]
		self.twist_cs = [cmath.cos(cmath.pi/n)*cmath.exp(2*cmath.pi*(i/n + 0.5/n)*1J) for i in range(n)]
		self.outer_cs = [2*cmath.cos(cmath.pi/n)**2*cmath.exp(2*cmath.pi*(i/n)*1J) for i in range(n)]
		x_outer_ps = intersectLines(self.m_inner_vs[0], self.outer_cs[0] - self.m_inner_vs[0], self.twist_vs[0], self.outer_vs[0] - self.twist_vs[0])
		self.x_outer_ps = [x_outer_ps*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		lpleat_vs = self.outer_cs[0] + abs(self.outer_vs[0] - self.outer_cs[0])/2**.5*cmath.exp(cmath.pi/4*1J)
		rpleat_vs = lpleat_vs.conjugate()
		self.lpleat_vs = [lpleat_vs*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		self.rpleat_vs = [rpleat_vs*cmath.exp(2*cmath.pi*i/n*1J) for i in range(n)]
		if r != 1:
			FlowerTowerUnit.__imul__(self, r)
		if z != 0J:
			FlowerTowerUnit.__iadd__(self, z)
	
	def transform_points(self, f):
		for points in (self.twist_vs, self.inner_vs, self.m_inner_vs, self.outer_vs, self.twist_cs, self.outer_cs, self.x_outer_ps, self.lpleat_vs, self.rpleat_vs):
			for i, z in enumerate(points):
				points[i] = f(z)
	
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
	
	def draw(self, screen):
		N = len(self.twist_vs)
		for i in range(N):
			pygame.draw.line(screen, 0xFF00FFFF, toScreen(self.twist_vs[i]), toScreen(self.twist_vs[(i+1)%N]))
			pygame.draw.line(screen, 0xFFFF00FF, toScreen(self.inner_vs[i]), toScreen(self.inner_vs[(i+1)%N]))
			pygame.draw.line(screen, 0xFF00FFFF, toScreen(self.inner_vs[i]), toScreen(self.twist_vs[(i+1)%N]))
			pygame.draw.line(screen, 0xFFFF00FF, toScreen(self.inner_vs[i]), toScreen(self.twist_cs[i]))
			pygame.draw.line(screen, 0xFF00FFFF, toScreen(self.m_inner_vs[i]), toScreen(self.twist_cs[i]))
			pygame.draw.line(screen, 0xFF00FFFF, toScreen(self.twist_vs[i]), toScreen(self.outer_cs[i]))
			pygame.draw.line(screen, 0xFFFF00FF, toScreen(self.m_inner_vs[i]), toScreen(self.twist_vs[(i+1)%N]))
			pygame.draw.line(screen, 0xFFFF00FF, toScreen(self.twist_vs[i]), toScreen(self.outer_vs[i]))
			pygame.draw.line(screen, 0xFFFF00FF, toScreen(self.outer_vs[i]), toScreen(self.outer_cs[(i+1)%N]))
			pygame.draw.line(screen, 0xFF00FFFF, toScreen(self.m_inner_vs[i]), toScreen(self.x_outer_ps[i]))
			pygame.draw.line(screen, 0xFFFF00FF, toScreen(self.outer_cs[i]), toScreen(self.x_outer_ps[i]))
			pygame.draw.line(screen, 0xFF00FFFF, toScreen(self.m_inner_vs[i]), toScreen(self.outer_cs[(i+1)%N]))
			pygame.draw.line(screen, 0xFF00FFFF, toScreen(self.lpleat_vs[i]), toScreen(self.outer_cs[i]))
			pygame.draw.line(screen, 0xFF00FFFF, toScreen(self.lpleat_vs[i]), toScreen(self.outer_vs[i]))
			pygame.draw.line(screen, 0xFF00FFFF, toScreen(self.rpleat_vs[i]), toScreen(self.outer_vs[(i-1)%N]))
			pygame.draw.line(screen, 0xFF00FFFF, toScreen(self.rpleat_vs[i]), toScreen(self.outer_cs[i]))
			pygame.draw.line(screen, 0xFFFF00FF, toScreen(self.rpleat_vs[i]), toScreen(self.lpleat_vs[i]))

twists = [FlowerTowerUnit(8, 0.06, -0.32-0.32J)]
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
	for i in range(3):
		pygame.draw.line(screen, 0xFF00FFFF, toScreen(twists[i].rpleat_vs[1]), toScreen(twists[i+1].lpleat_vs[5]))
		pygame.draw.line(screen, 0xFF00FFFF, toScreen(twists[i].lpleat_vs[1]), toScreen(twists[i+1].rpleat_vs[5]))
		pygame.draw.line(screen, 0xFFFF00FF, toScreen(twists[i].outer_vs[0]), toScreen(twists[i+1].outer_vs[5]))
		pygame.draw.line(screen, 0xFFFF00FF, toScreen(twists[i].outer_vs[1]), toScreen(twists[i+1].outer_vs[4]))
	for (i, j) in ((4, 1), (2, 5)):
		pygame.draw.line(screen, 0xFF00FFFF, toScreen(twists[i].rpleat_vs[7]), toScreen(twists[j].lpleat_vs[3]))
		pygame.draw.line(screen, 0xFF00FFFF, toScreen(twists[i].lpleat_vs[7]), toScreen(twists[j].rpleat_vs[3]))
		pygame.draw.line(screen, 0xFFFF00FF, toScreen(twists[i].outer_vs[7]), toScreen(twists[j].outer_vs[2]))
		pygame.draw.line(screen, 0xFFFF00FF, toScreen(twists[i].outer_vs[6]), toScreen(twists[j].outer_vs[3]))
	pygame.display.flip()
	#time.sleep(100)

