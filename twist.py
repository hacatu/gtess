import sys
import pygame.freetype

from geometry import *
from cp_builders import *
from application import *

points1 = ComposablePointList().add_points([(.5, 0), (0, 1)])
points1.i_replicate_polar(ComposablePointList().add_points([ucircPoint(2*i*np.pi/2) for i in range(2)]))
edges1 = ComposableEdgeList(points1).add_edges_MVR(
	mountains=[(1,2)],
	valleys=[(0,1)]
).i_replicate(2, 2)
pcp_builder = PartialCreasePatternBuilder().with_vertices(points1.vs).with_edges(edges1.edges)

def makeFlowerTowerBuilder(num_sides : int) -> PartialCreasePatternBuilder:
	"""Create a PartialCreasePatternBuilder for a flower tower with a given number of sides (should be at least 7 sides)."""
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

pygame.init()
screen = pygame.display.set_mode((WIDTH, WIDTH), pygame.RESIZABLE)
pygame.display.set_caption("Gtess")
pgft_font = pygame.freetype.SysFont("Monospace", 28)
scene = Scene(screen, (WIDTH, WIDTH))

scene.num_sides = 12
pcp_builder = makeFlowerTowerBuilder(scene.num_sides)
crease_pattern = pcp_builder.build()
crease_pattern *= .2

scene.addObject(crease_pattern)

class RootInputContext(InputContext):
	def __init__(self, scene):
		super().__init__()
		self.scene = scene
	
	def _handleWindowEvent(self, e : pygame.event.EventType) -> bool :
		if e.type == pygame.QUIT or (e.type == pygame.KEYDOWN and e.key == pygame.K_ESCAPE):
			pygame.quit()
			sys.exit(0)
		elif e.type == pygame.VIDEOEXPOSE:
			self.scene.dirty = True
		elif e.type == pygame.VIDEORESIZE:
			self.scene.screen_size = e.size
		else:
			return False
		return True

	def handleEvent(self, e):
		if self._handleWindowEvent(e) or super().handleEvent(e):
			return True
		if e.type == pygame.KEYDOWN:
			if e.key == pygame.K_RIGHT:
				self.scene.num_sides += 1
				pcp_builder = makeFlowerTowerBuilder(scene.num_sides)
				crease_pattern = pcp_builder.build()
				crease_pattern *= .2
				self.scene.objects[0] = crease_pattern
				self.scene.dirty = True
			elif e.key == pygame.K_LEFT and self.scene.num_sides > 7:
				self.scene.num_sides -= 1
				pcp_builder = makeFlowerTowerBuilder(self.scene.num_sides)
				crease_pattern = pcp_builder.build()
				crease_pattern *= .2
				self.scene.objects[0] = crease_pattern
				self.scene.dirty = True
			elif e.key == pygame.K_SPACE:
				if isinstance(self.scene.objects[0], Steppable):
					self.scene.objects[0].step()
					self.scene.dirty = True
		if e.type == pygame.MOUSEMOTION and e.buttons[1]:
			dsx, dsy = e.rel
			transform = self.scene.view_transform
			transform = np.vstack([transform, [0, 0, 1]])
			transform = np.array([[1, 0, dsx], [0, 1, dsy], [0, 0, 1]]) @ transform
			self.scene.setTransform(transform[:2])
		elif e.type == pygame.MOUSEWHEEL and e.y:
			# TODO; translate before & after so the zoom is centered at the mouse
			scale = 1.1**e.y
			transform = self.scene.view_transform
			transform = np.vstack([transform, [0, 0, 1]])
			transform = transform @ np.array([[scale, 0, 0], [0, scale, 0], [0, 0, 1]])
			self.scene.setTransform(transform[:2])
		else:
			return False
		return True

class ToolLineInputContext(InputContext):
	def __init__(self, scene):
		super().__init__()
		self.scene = scene
		self.ribbon_text = "Current tool: draw line"
	
	def handleEvent(self, e):
		if e.type == pygame.MOUSEBUTTONDOWN and e.button == MouseButton.LEFT:
			self.scene.addObject(UserLine(self.scene.fromScreen(e.pos)))
			self.scene.dirty = True
		elif e.type == pygame.MOUSEMOTION and e.buttons[0]:
			user_line = next((d for d in reversed(self.scene.objects) if isinstance(d, UserLine)), None)
			if user_line is not None and user_line.active:
				user_line.move_terminus(1, self.scene.fromScreen(e.pos))
				self.scene.dirty = True
		elif e.type == pygame.MOUSEBUTTONUP and e.button == MouseButton.LEFT:
			user_line = next((d for d in reversed(self.scene.objects) if isinstance(d, UserLine)), None)
			if user_line is not None:
				user_line.finalize()
		else:
			return False
		return True

toolInputContextMap = {
	"tool.line": ToolLineInputContext(scene)
}

root_input_ctx = RootInputContext(scene)
root_input_ctx.setActiveChild(toolInputContextMap["tool.line"])

while True:
	for e in pygame.event.get():
		root_input_ctx.handleEvent(e)
	if scene.dirty:
		screen.fill(0xFF000000)
		scene.draw()
		if root_input_ctx.active_child is not None:
			pgft_font.render_to(screen, (0, scene.screen_size[1] - 28), root_input_ctx.active_child.ribbon_text, 0xFFFFFFFF)
		pygame.display.flip()

