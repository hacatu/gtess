import pygame.freetype
from collections import defaultdict, deque
import sys

from geometry import *
from cp_builders import *

class UserLine(M3Transformable, Drawable):
	def __init__(self, start : np.ndarray):
		start = start.reshape((2,1))
		self.vs = np.hstack([start, start])
		self.active = True
	
	def move_terminus(self, i : int, z : np.ndarray) -> 'UserLine':
		self.vs[:,i] = z.reshape(2,)
		return self
	
	def finalize(self):
		self.active = False

	def apply_M3(self, m : np.ndarray) -> 'UserLine':
		if isinstance(m, np.ndarray) and m.shape == (2, 3):
			self.vs = m @ np.vstack([self.vs, np.ones((1, self.vs.shape[0]))])
		return self
	
	def draw(self, scene : Scene) -> None:
		pygame.draw.line(scene.screen, 0xFF00FF00, scene.toScreen(self.vs[:,0]), scene.toScreen(self.vs[:,1]))

"""
points1 = ComposablePointList().add_points([(0, 3**-.5/2), (.25, 3**-.5/2), (.5, 3**-.5/2)])
points1.i_replicate_polar(ComposablePointList().add_points([(1, 0), ucircPoint(2*np.pi/3), ucircPoint(-2*np.pi/3)]))
edges1 = ComposableEdgeList(reference_cpl=points1).add_edges_MVR(
	mountains=[(0,6)],
	valleys=[(1,6)],
	raws=[(0,1),(1,2),(2,6)]
).i_replicate(3, 3)
pcp_builder = PartialCreasePatternBuilder().with_vertices(points1.vs).with_edges(edges1.edges)
"""

"""
points1 = ComposablePointList().add_points([(.1, 0), (1, 0)])
points1.i_replicate_polar(ComposablePointList().add_points([ucircPoint(2*i*np.pi/3) for i in range(3)]))
edges1 = ComposableEdgeList(points1).add_edges_MVR(
	mountains=[(0, 1)],
	valleys=[(0,2)],
	raws=[(1,3)]
).i_replicate(2, 3)
pcp_builder = PartialCreasePatternBuilder().with_vertices(points1.vs).with_edges(edges1.edges)
"""

points1 = ComposablePointList().add_points([(.5, 0), (0, 1)])
points1.i_replicate_polar(ComposablePointList().add_points([ucircPoint(2*i*np.pi/2) for i in range(2)]))
edges1 = ComposableEdgeList(points1).add_edges_MVR(
	mountains=[(1,2)],
	valleys=[(0,1)]
).i_replicate(2, 2)
pcp_builder = PartialCreasePatternBuilder().with_vertices(points1.vs).with_edges(edges1.edges)

class FTFPoints(M3Transformable, Drawable, Steppable):
	def __init__(self, num_sides):
		cn = np.cos(np.pi/num_sides)
		diag = np.array(ucircPoint(np.pi/num_sides)).T
		points1 = ComposablePointList().add_points([(1, 0), cn*diag, 0.5/cn*diag, 2, 2/cn*diag])
		def stepper():
			points1.add_points([2*points1.vs[:,1] - points1.vs[:,2]])
			yield 
			points1.add_points([intersectLines(points1.vs[:,0], points1.vs[:,4] - points1.vs[:,0], points1.vs[:,5], points1.vs[:,3] - points1.vs[:,5])])
			yield
			points1.add_points([M3(1+1J)(points1.vs[:,3]) - M3(1J)(points1.vs[:,4])])
			yield
			points1.add_points([(points1.vs[:,4] + points1.vs[:,7])/2])
			yield
			points1.add_points([M3.Conjugate(points1.vs[:,8])])
			yield
			tmp = v2(-1+1J)*np.linalg.norm(points1.vs[:,1] - points1.vs[:,0])/2**.5
			yield
			points1.add_points([points1.vs[:,7] + tmp])
			yield
			points1.add_points([M3.Conjugate(points1.vs[:,10]), v2((points1.vs[0,10], points1.vs[1,8])), v2((points1.vs[0,10], points1.vs[1,9])), v2((points1.vs[0,10], points1.vs[1,4]))])
			yield
			points1.add_points([M3.Conjugate(points1.vs[:,14]), points1.vs[0,10]/cn*diag])
		self.points1 = points1
		self._stepper = stepper()

	def apply_M3(self, m):
		self.points1.apply_M3(m)
		return self
	
	def draw(self, scene):
		self.points1.draw(scene)
	
	def step(self):
		next(self._stepper)
	
def makeFlowerTowerBuilder(num_sides : int) -> PartialCreasePatternBuilder:
	cn = np.cos(np.pi/num_sides)
	diag = np.array(ucircPoint(np.pi/num_sides)).T
	points1 = ComposablePointList().add_points([(1, 0), cn*diag, 0.5/cn*diag, 2, 2/cn*diag])
	points1.add_points([2*points1.vs[:,1] - points1.vs[:,2]])
	points1.add_points([intersectLines(points1.vs[:,0], points1.vs[:,4] - points1.vs[:,0], points1.vs[:,5], points1.vs[:,3] - points1.vs[:,5])])
	points1.add_points([M3(1+1J)(points1.vs[:,3]) - M3(1J)(points1.vs[:,4])])
	points1.add_points([(points1.vs[:,4] + points1.vs[:,7])/2])
	points1.add_points([M3.Conjugate(points1.vs[:,8])])
	tmp = v2(-1+1J)*np.linalg.norm(points1.vs[:,1] - points1.vs[:,0])/2**.5
	#print(tmp)
	points1.add_points([points1.vs[:,7] + tmp])
	points1.add_points([M3.Conjugate(points1.vs[:,10]), v2((points1.vs[0,10], points1.vs[1,8])), v2((points1.vs[0,10], points1.vs[1,9])), v2((points1.vs[0,10], points1.vs[1,4]))])
	points1.add_points([M3.Conjugate(points1.vs[:,14]), points1.vs[0,10]/cn*diag])
	ostride = points1.vs.shape[1]*(num_sides - 1)
	points1.i_replicate_polar(ComposablePointList().add_points([ucircPoint(2*np.pi*i/num_sides) for i in range(num_sides)]))
	edges1 = ComposableEdgeList(points1).add_edges_MVR(
		mountains=[(2, 2 + ostride), (1, 2), (0, 6), (0, 5 + ostride), (3, 6), (4, 6), (3, 4 + ostride), (15, 4 + ostride), (8, 9), (4, 14)],
		valleys=[(0, 2 + ostride), (0, 1), (0, 1 + ostride), (1, 5), (0, 3), (3, 5 + ostride), (5, 6), (4, 8), (3, 8), (8, 12), (3, 9), (9, 13), (9, 4 + ostride)],
		raws=[(7, 10), (10, 12), (12, 14), (14, 16), (7, 11), (11, 13), (13, 15), (15, 16 + ostride)]
	).i_replicate(17, num_sides)
	return PartialCreasePatternBuilder().with_vertices(points1.vs).with_edges(edges1.edges)

num_sides = 12
pcp_builder = makeFlowerTowerBuilder(num_sides)
dirty = True
#crease_pattern = pcp_builder.build()
crease_pattern = pcp_builder.build()
#crease_pattern = FTFPoints(12)
crease_pattern *= .2
#crease_pattern = FlowerTowerBuilder(num_sides).scale(.2).build()
user_lines = []
ribbon_text = "Current tool: draw line"

pygame.init()
screen = pygame.display.set_mode((WIDTH, WIDTH))
pygame.display.set_caption("Gtess")
pgft_font = pygame.freetype.SysFont("Monospace", 28)
scene = Scene(screen, (WIDTH, WIDTH))

while True:
	for e in pygame.event.get():
		if e.type == pygame.QUIT or (e.type == pygame.KEYDOWN and e.key == pygame.K_ESCAPE):
			pygame.quit()
			sys.exit(0)
		elif e.type == pygame.VIDEOEXPOSE:
			dirty = True
		elif e.type == pygame.KEYDOWN:
			if e.key == pygame.K_RIGHT:
				num_sides += 1
				pcp_builder = makeFlowerTowerBuilder(num_sides)
				crease_pattern = pcp_builder.build()
				crease_pattern *= .2
				dirty = True
			elif e.key == pygame.K_LEFT and num_sides > 7:
				num_sides -= 1
				pcp_builder = makeFlowerTowerBuilder(num_sides)
				crease_pattern = pcp_builder.build()
				crease_pattern *= .2
				dirty = True
			elif e.key == pygame.K_SPACE:
				if isinstance(crease_pattern, Steppable):
					crease_pattern.step()
					dirty = True
		elif e.type == pygame.MOUSEBUTTONDOWN and e.button == 1:
			user_lines.append(UserLine(scene.fromScreen(e.pos)))
			dirty = True
		elif e.type == pygame.MOUSEMOTION and 1 in e.buttons:
			if user_lines and user_lines[-1].active:
				user_lines[-1].move_terminus(1, scene.fromScreen(e.pos))
				dirty = True
		elif e.type == pygame.MOUSEBUTTONUP and e.button == 1:
			if user_lines:
				user_lines[-1].finalize()
	if dirty:
		screen.fill(0xFF000000)
		crease_pattern.draw(scene)
		pgft_font.render_to(screen, (0, WIDTH - 28), ribbon_text, 0xFFFFFFFF)
		for line in user_lines:
			line.draw(scene)
		pygame.display.flip()
		dirty = False

