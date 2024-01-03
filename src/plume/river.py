#! /usr/bin/env python
from __future__ import annotations


class River:
    def __init__(
        self,
        velocity: float = 1.0,
        width: float = 1.0,
        depth: float = 1.0,
        angle: float = 0.0,
        loc: tuple[float, float] = (0.0, 0.0),
    ):
        self._velocity = velocity
        self._width = width
        self._depth = depth
        self._angle = angle
        self._x0, self._y0 = loc

    @property
    def x0(self) -> float:
        return self._x0

    @property
    def y0(self) -> float:
        return self._y0

    @property
    def velocity(self) -> float:
        return self._velocity

    @property
    def width(self) -> float:
        return self._width

    @property
    def depth(self) -> float:
        return self._depth

    @property
    def angle(self) -> float:
        return self._angle

    @property
    def discharge(self) -> float:
        return self.velocity * self.width * self.depth
