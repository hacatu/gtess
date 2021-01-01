import numpy as np
import itertools as itt
import pygame
import inspect
import sys
import cmath
from collections import defaultdict

class AutoDict(dict):
	def __init__(self, make):
		self.make = make
		super().__init__()
	def __missing__(self, key):
		item = self.make(key)
		super().__setitem__(key, item)
		return item

def _loadResource(filename):
	#print("Loading " + filename)
	res = pygame.image.load(filename)
	res.convert()
	#print(res)
	return res

class UIStateTransition(Exception):
	def __init__(self, target, event):
		self.target = target
		self.event = event

class UIStateManager:
	def __init__(self, viewport):
		self.viewport = viewport
		self.state_map = {
			name: clazz(self) for (name, clazz) in
				inspect.getmembers(sys.modules[__name__], inspect.isclass) if
				issubclass(clazz, UIState)
		}
		self.resource_map = AutoDict(_loadResource)
	
	def getStateInstance(self, clazz):
		return self.state_map[clazz.__name__]
	
	def getResource(self, filename):
		return self.resource_map[filename]

class UIDrawable:
	def __init__(self):
		self.dirty = True
		self.children = []
	
	def draw(self, sman):
		self.dirty = False
		for drawable in self.children:
			if drawable.dirty:
				drawable.draw(sman)
	
	def markDirty(self):
		if not self.dirty:
			self.dirty = True
			for drawable in self.children:
				drawable.markDirty()

class ModelRect(UIDrawable):
	def __init__(self, bl, tr):
		super().__init__()
		self.bl = bl
		self.tr = tr
	
	def draw(self, sman):
		left, top = sman.viewport.toScreen(self.bl)
		right, bottom = sman.viewport.toScreen(self.tr)
		pygame.draw.rect(sman.viewport.screen, 0xFFFF00FF, (left, top, right - left, bottom - top))
		super().draw(sman)

class UIIcon(UIDrawable):
	def __init__(self, bl, filename):
		super().__init__()
		self.bl = bl
		self.filename = filename
	
	def draw(self, sman):
		sman.viewport.screen.blit(sman.getResource(self.filename), self.bl)
		super().draw(sman)

class UILayer(UIDrawable):
	def __init__(self, sman):
		super().__init__()
		self.sman = sman
	
	def draw(self):
		super().draw(self.sman)

class UIViewport:
	def __init__(self, screen, offset, scale):
		self.screen = screen
		self.width = screen.get_width()
		self.height = screen.get_height()
		self.offset = offset
		self.scale = scale
	
	def toScreen(self, p):
		ps = p*self.scale + self.offset
		return (int(ps.real), self.height - int(ps.imag))
	
	def fromScreen(self, p):
		ps = complex(p[0], self.height - p[1])
		return (ps - self.offset)/self.scale

class UIState:
	def __init__(self, sman):
		self.sman = sman
		self.layers = []
	
	def handleEvent(self, e):
		if e.type == pygame.QUIT or (e.type == pygame.KEYDOWN and e.key == pygame.K_ESCAPE):
			raise UIStateTransition(UIState_Stop, e)
		return False
	
	def enter(self, event):
		pass
	
	def exit(self):
		pass

class UIState_Normal(UIState):
	def handleEvent(self, e):
		if super().handleEvent(e):
			return True
		if e.type == pygame.MOUSEBUTTONDOWN and e.button == 3:
			raise UIStateTransition(UIState_ShowMenu, e)
		return False

class UIState_ShowMenu(UIState):
	def handleEvent(self, e):
		if super().handleEvent(e):
			return True
		return False
	
	def enter(self, event):
		menu_layer = UILayer(self.sman)
		#menu_layer.children.append(ModelRect(.1+.1j, .11+.11j))
		#menu_layer.children.append(ModelRect(.4+.4j, .41+.41j))
		#menu_layer.children.append(ModelRect(.9+.9j, .91+.91j))
		menu_layer.children.append(UIIcon((10, 10), "uwu.png"))
		self.layers = [menu_layer]
	
class UIState_Stop(UIState):
	def enter(self, event):
		pygame.quit()
		sys.exit()

def main():
	pygame.init()
	screen = pygame.display.set_mode((600, 400))
	pygame.display.set_caption("Gtess")
	
	viewport = UIViewport(screen, 0+0j, 400.)
	sman = UIStateManager(viewport)
	uiState = sman.getStateInstance(UIState_Normal)
	uiState.enter(None)
	
	while True:
		transition = None
		try:
			for e in pygame.event.get():
				uiState.handleEvent(e)
		except UIStateTransition as t:
			transition = t
		if transition is not None:
			uiState.exit()
			uiState = sman.getStateInstance(transition.target)
			uiState.enter(transition.event)
		updated = False
		for layer in uiState.layers:
			if layer.dirty:
				layer.draw()
				updated = True
		if updated:
			pygame.display.flip()

if __name__ == "__main__":
	try:
		main()
	except SystemExit:
		pass #good golly gravey how daft is it that this prints a stack trace

