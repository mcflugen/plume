#! /usr/bin/env python
import numpy as np

EPSILON = np.finfo(float).eps
PI_4 = np.pi / 4.0


def unit_vector(v):
    v = np.asarray(v).reshape((2, -1))
    v_abs = np.linalg.norm(v, axis=0)
    return np.divide(v, v_abs, where=v_abs > 0.0, out=np.zeros_like(v)).squeeze()


def nearest_point_on_ray(points, origin=(0.0, 0.0), angle=PI_4, out=None):
    """Find the nearest point on a ray.

    Parameters
    ----------
    points : tuple of float
        Points as (x, y).
    origin : tuple of float, optional
        The start of the ray.
    angle : float, optional
        Angle of the ray.
    """
    points = np.asarray(points, dtype=float).reshape((2, -1))
    origin = np.asarray(origin, dtype=float).reshape((2, -1))

    if out is None:
        out = np.empty_like(points)
    out.shape = (2, -1)
    out.fill(0.0)

    u = points - origin

    u_bar = unit_vector(u)
    v_bar = np.asarray((np.cos(angle), np.sin(angle)))

    v_dot_u = np.asarray(np.dot(v_bar, u_bar)).reshape((-1,))

    close_to_line = v_dot_u > EPSILON

    r = v_dot_u[close_to_line] * np.linalg.norm(u[:, close_to_line], axis=0)

    if np.any(close_to_line):
        out[0, close_to_line] = r * v_bar[0]
        out[1, close_to_line] = r * v_bar[1]

    out += origin

    return out.squeeze()


class PlumeCenterline:

    N = 0.37

    def __init__(
        self,
        river_width,
        river_velocity=1.0,
        ocean_velocity=1.0,
        river_angle=0.0,
        river_loc=(0.0, 0.0),
    ):
        self._river_width = river_width
        self._river_angle = river_angle
        self._river_velocity = river_velocity
        self._river_x0, self._river_y0 = river_loc
        self._ocean_velocity = ocean_velocity
        if self.is_straight:
            self._c = None
        else:
            self._c = np.fabs(1.53 * 0.909 * river_velocity / ocean_velocity)

    @property
    def is_straight(self):
        return np.fabs(self.ocean_velocity) < 1e-12

    @property
    def river_angle(self):
        """Angle river mouth makes with the coast."""
        return self._river_angle

    @property
    def river_width(self):
        return self._river_width

    @property
    def river_velocity(self):
        return self._river_velocity

    @property
    def x0(self):
        return self._river_x0

    @property
    def y0(self):
        return self._river_y0

    @property
    def ocean_velocity(self):
        return self._ocean_velocity

    def x(self, y):
        if self.is_straight:
            return y / np.tan(self.river_angle)
        else:
            return (
                self.river_width
                * self._c
                * np.power(
                    np.sign(self.ocean_velocity) * (y - self.y0) / self.river_width,
                    self.N,
                )
            ) + self.x0

    def y(self, x):
        if self.is_straight:
            return np.tan(self.river_angle) * x
        else:
            return (
                np.sign(self.ocean_velocity)
                * (
                    self.river_width
                    * np.power((x - self.x0) / self.river_width / self._c, 1.0 / self.N)
                )
                + self.y0
            )

    def is_function_of_x(self):
        if self.is_straight:
            return np.fabs(np.cos(self.river_angle)) > 1e-12

        if (
            np.fabs(np.cos(self.river_angle)) <= 1e-12
            or np.tan(self.river_angle) * self.ocean_velocity <= 0.0
        ):
            return True
        else:
            return False

    def is_function_of_y(self):
        if self.is_straight:
            return np.fabs(np.sin(self.river_angle)) > 1e-12

        if (
            np.fabs(np.cos(self.river_angle)) <= 1e-12
            or np.tan(self.river_angle) * self.ocean_velocity >= 0.0
        ):
            return True
        else:
            return False

    def path_length(self, bounds):
        from .ext.centerline import path_lengths

        bounds = np.asarray(bounds, dtype=float).reshape((-1, 2))
        lengths = np.empty(len(bounds), dtype=float)

        if self.is_function_of_x():
            angle = self.river_angle
            shift = self.x0
        else:
            angle = self.river_angle - np.pi * 0.5
            shift = self.y0

        if self.ocean_velocity < 0.0:
            angle *= -1.0

        if self.is_straight:
            lengths[:] = np.abs(np.diff(bounds, axis=1) / np.cos(angle)).squeeze()
        else:
            path_lengths(
                bounds - shift, self.river_width, self._c, self.N, angle, lengths
            )

        return lengths

    def r(self, x, x0, y0):
        return np.sqrt((x - x0) ** 2.0 + (self.y(x) - y0) ** 2.0)

    @staticmethod
    def rotate_points(x, y, angle):
        return (
            x * np.cos(angle) - y * np.sin(angle),
            x * np.sin(angle) + y * np.cos(angle),
        )

    def nearest_point(self, points):
        from .ext.centerline import nearest_points

        points = np.asarray(points).reshape((-1, 2))
        nearest = np.empty_like(points)

        if self.is_straight:
            nearest_point_on_ray(
                points.T,
                angle=self.river_angle,
                origin=(self.x0, self.y0),
                out=nearest.T,
            )
            return nearest

        angle = self.river_angle * np.sign(self.ocean_velocity)

        nearest_points(
            (points - (self.x0, self.y0)) * (1.0, np.sign(self.ocean_velocity)),
            self.river_width,
            self._c,
            self.N,
            angle,
            nearest,
        )

        return nearest * (1.0, np.sign(self.ocean_velocity)) + (self.x0, self.y0)

    def distance_to(self, points):
        points_on_centerline = self.nearest_point(points)
        distance = np.sqrt(np.power(points_on_centerline - points, 2).sum(axis=1))
        return distance
