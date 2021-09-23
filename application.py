import abc
from enum import IntEnum
import pygame
import numpy as np

from geometry import *

class MouseButton(IntEnum):
	LEFT = 1
	MIDDLE = 2
	RIGHT = 3
	WHEEL_UP = 4
	WHEEL_DOWN = 5
	NAV_BACK = 6
	NAV_FORWARD = 7

class InputContext(abc.ABC):
	"""Abstract class to implement nestable input handling.
	Each input handling "layer" (root layer for exit, window reveal/resize,
	etc; layer for active tool; layer for ui toolbars; etc) should have its
	own subclass of this class, and InputContexts should be nested down
	from the root using setActiveChild"""
	def __init__(self):
		self.active_child = None

	def setActiveChild(self, ctx : 'InputContext' = None) -> None :
		"""Set another input context as the active child of this one.
		Any input events not handled by this context will be passed down
		to the child, recursively if need be."""
		self.active_child = ctx
	
	def handleEvent(self, e : pygame.event.EventType) -> bool :
		"""Default event handler.  Subclasses should override this, handle events,
		and then call this default implementation through super().  Return True
		if the event was handled (including if it is ignored and should not be passed down),
		or False if the event was not handled (or if it should be passed down in general)."""
		return self.active_child is not None and self.active_child.handleEvent(e)

class UserLine(M3Transformable, Drawable):
	"""A line drawn by the user in the scene.  Currently useful only to visualize things."""
	def __init__(self, start : np.ndarray):
		start = start.reshape((2,1))
		self.vs = np.hstack([start, start])
		self.active = True
	
	def move_terminus(self, i : int, z : np.ndarray) -> 'UserLine':
		"""Change the coordinates of the start (i==0) or end (i==1) of the line.
		z is given in model coordinates, NOT screen coordinates."""
		self.vs[:,i] = z.reshape(2,)
		return self
	
	def finalize(self):
		"""Mark the line as done (when the user releases the mouse while using this tool."""
		self.active = False

	def apply_M3(self, m : np.ndarray) -> 'UserLine':
		if isinstance(m, np.ndarray) and m.shape == (2, 3):
			self.vs = m @ np.vstack([self.vs, np.ones((1, self.vs.shape[0]))])
		return self
	
	def draw(self, scene : Scene) -> None:
		pygame.draw.line(scene.screen, 0xFF00FF00, scene.toScreen(self.vs[:,0]), scene.toScreen(self.vs[:,1]))

