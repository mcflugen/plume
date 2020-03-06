#! /usr/bin/env python

import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

from plume.centerline import nearest_point_on_ray, unit_vector


def test_unit_vector_with_scalars():
    """Test passing scalars."""
    u = unit_vector((1.0, 1.0))
    assert_array_equal(u, np.full(2, 1.0 / np.sqrt(2.0)))


def test_unit_vector_with_arrays():
    """Test passing scalars."""
    u = unit_vector([[1.0, -1], [1.0, 1.0]])
    assert_array_equal(u, np.array([[1, -1], [1.0, 1.0]]) / np.sqrt(2.0))


def test_unit_vector_on_x_axis():
    """Test vector on +- x-axis."""
    u = unit_vector([[10.0, -2.0], [0.0, 0.0]])
    assert_array_equal(u, ((1.0, -1.0), (0.0, 0.0)))


def test_unit_vector_on_y_axis():
    """Test vector on +- y-axis."""
    u = unit_vector([[0.0, 0.0], [0.5, -1e6]])
    assert_array_equal(u, ((0.0, 0.0), (1.0, -1.0)))


def test_nearest_with_scalars():
    """Test passing scalars."""
    pts = nearest_point_on_ray((1.0, 0.0), (0.0, 0.0), np.deg2rad(45.0))
    assert_array_almost_equal(pts, [0.5, 0.5])


def test_nearest_with_arrays():
    angle = np.deg2rad(45.0)

    pts = nearest_point_on_ray(((1.0,), (0.0,)), (0.0, 0.0), angle)
    assert_array_almost_equal(pts, [0.5, 0.5])

    pts = nearest_point_on_ray(((1.0, 0.0), (0.0, 2.0)), (0.0, 0.0), angle)
    assert_array_almost_equal(pts, [[0.5, 1.0], [0.5, 1.0]])


def test_nearest_moving_origin():
    """Test using the origin keyword."""
    pts = nearest_point_on_ray(((1.0, 1.0, 1.0), (1.0, 2.0, 10.0)), origin=(1.0, 2.0))
    assert_array_almost_equal(pts, [[1.0, 1.0, 5.0], [2.0, 2.0, 6.0]])


def test_nearest_changing_angle():
    """Test using the angle keyword."""
    pts = nearest_point_on_ray(((1.0, 1.0, 0.0), (1.0, -1.0, -10)), angle=-np.pi / 4.0)
    assert_array_almost_equal(pts, [[0.0, 1.0, 5.0], [0.0, -1.0, -5.0]])


def test_nearest_with_defaults():
    """Test with default optional arguments."""
    pts = nearest_point_on_ray(((0.0, 3.0), (0.0, 1.0)))
    assert_array_almost_equal(pts, [[0.0, 2.0], [0.0, 2.0]])


def test_nearest_to_x_axis():
    """Test if the ray is the +- x-axis."""
    for angle in (0.0, 2 * np.pi):
        pts = nearest_point_on_ray(((1.0, -1.0), (1.0, -1.0)), angle=angle)
        assert_array_almost_equal(pts, [[1.0, 0.0], [0.0, 0.0]])

    for angle in (-np.pi, np.pi):
        pts = nearest_point_on_ray(((1.0, -1.0), (1.0, -1.0)), angle=angle)
        assert_array_almost_equal(pts, [[0.0, -1.0], [0.0, 0.0]])


def test_nearest_to_y_axis():
    """Test if the ray is the +- y-axis."""
    for angle in (np.pi / 2.0, 5.0 * np.pi / 2.0):
        pts = nearest_point_on_ray(((1.0, -1.0), (1.0, -1.0)), angle=angle)
        assert_array_almost_equal(pts, [[0.0, 0.0], [1.0, 0.0]])

    for angle in (-np.pi / 2.0, 3.0 * np.pi / 2.0):
        pts = nearest_point_on_ray(((1.0, -1.0), (1.0, -1.0)), angle=angle)
        assert_array_almost_equal(pts, [[0.0, 0.0], [0.0, -1.0]])


def test_nearest_on_line():
    """Test if the point is already on the line."""
    pts = nearest_point_on_ray(((1.0, 0.0), (1.0, 0.0)))
    assert_array_almost_equal(pts, [[1.0, 0.0], [1.0, 0.0]])


def test_nearest_behind_origin():
    """Test if the nearest point is the origin."""
    pts = nearest_point_on_ray((-1.0, 0.0))
    assert_array_almost_equal(pts, [0.0, 0.0])
