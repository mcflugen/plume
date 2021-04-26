#! /usr/bin/env python

from timeit import default_timer as timer

import numpy as np
from landlab import RasterModelGrid
from numpy.testing import assert_array_almost_equal, assert_array_equal

from plume import Plume


def test_caching():
    """Test caching."""
    grid = RasterModelGrid((500, 500), xy_spacing=(100.0, 100.0))

    params = {
        "river_velocity": 2.5,
        "river_width": 50.0,
        "river_depth": 5.0,
        "river_loc": (0.0, 25000.0),
        "ocean_velocity": 0.015,
        "river_angle": np.deg2rad(0.0),
    }
    plume = Plume(grid, **params)
    time_0 = timer()
    plume.run_one_step()
    time_1 = timer()
    plume.run_one_step()
    time_2 = timer()

    first_time = time_1 - time_0
    second_time = time_2 - time_1

    assert plume._run_one_step.cache_info().hits == 1
    assert second_time < 0.5 * first_time


def test_change_velocity():
    params = {
        "river_velocity": 2.5,
        "river_width": 50.0,
        "river_depth": 5.0,
        "river_loc": (0.0, 25000.0),
        "ocean_velocity": 0.015,
        "river_angle": np.deg2rad(0.0),
    }

    conc = {}

    params["river_velocity"] = 2.5
    grid = RasterModelGrid((500, 500), xy_spacing=(100.0, 100.0))

    plume = Plume(grid, **params)
    plume.run_one_step()
    conc[2.5] = grid.at_node["sediment~suspended__mass_concentration"].copy()

    params["river_velocity"] = 0.5
    grid = RasterModelGrid((500, 500), xy_spacing=(100.0, 100.0))
    plume = Plume(grid, **params)
    plume.run_one_step()
    conc[0.5] = grid.at_node["sediment~suspended__mass_concentration"].copy()

    assert np.any(conc[0.5] != conc[2.5])

    plume.river_velocity = 2.5
    plume.run_one_step()
    c = grid.at_node["sediment~suspended__mass_concentration"].copy()

    assert np.any(c != conc[0.5])
    assert_array_equal(c, conc[2.5])
