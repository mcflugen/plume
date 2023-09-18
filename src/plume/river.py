#! /usr/bin/env python


class River:
    def __init__(self, velocity=1.0, width=1.0, depth=1.0, angle=0.0, loc=(0.0, 0.0)):
        self._velocity = velocity
        self._width = width
        self._depth = depth
        self._angle = angle
        self._x0, self._y0 = loc

    @property
    def x0(self):
        return self._x0

    @property
    def y0(self):
        return self._y0

    @property
    def velocity(self):
        return self._velocity

    @property
    def width(self):
        return self._width

    @property
    def depth(self):
        return self._depth

    @property
    def angle(self):
        return self._angle

    @property
    def discharge(self):
        return self.velocity * self.width * self.depth
