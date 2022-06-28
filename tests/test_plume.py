import pytest
from landlab import RasterModelGrid

from plume import Plume


def test_plume_deposit_thickness():
    grid = RasterModelGrid((3, 500), xy_spacing=(100.0, 10000.0))

    plume = Plume(grid)

    assert grid.at_node["sediment_deposit__thickness"] == pytest.approx(0.0)

    plume.run_one_step()
    actual = grid.at_node["sediment_deposit__thickness"].copy()

    assert not any(actual < 0.0)
    assert any(actual > 0.0)

    grid.at_node["sediment_deposit__thickness"][:] = 1.0

    plume = Plume(grid)
    assert grid.at_node["sediment_deposit__thickness"] == pytest.approx(1.0)

    plume.run_one_step()
    assert grid.at_node["sediment_deposit__thickness"] == pytest.approx(actual + 1.0)
