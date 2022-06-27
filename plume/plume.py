#! /usr/bin/env python
import numpy as np
from landlab import Component

from .centerline import PlumeCenterline
from .river import River

SQRT_PI = np.sqrt(np.pi)
SQRT_TWO = np.sqrt(2.0)
SECONDS_PER_DAY = 60 * 60 * 24.0


class Plume(Component):

    _name = "Plume"
    _unit_agnostic = True
    _info = {
        "sediment__removal_rate": {
            "dtype": "float",
            "intent": "in",
            "optional": True,
            "units": "1 / d",
            "mapping": "grid",
            "doc": "removal rate of sediment carried by the plume",
        },
        "sediment__bulk_density": {
            "dtype": "float",
            "intent": "in",
            "optional": True,
            "units": "kg / m^3",
            "mapping": "grid",
            "doc": "bulk density of sediment deposited by the plume",
        },
        "tracer~conservative__mass_concentration": {
            "dtype": "float",
            "intent": "out",
            "optional": False,
            "units": "kg / m^3",
            "mapping": "node",
            "doc": "concentration of a conservative tracer",
        },
        "sediment~suspended__mass_concentration": {
            "dtype": "float",
            "intent": "out",
            "optional": False,
            "units": "kg / m^3",
            "mapping": "node",
            "doc": "concentration of suspended sediment in the plume",
        },
        "sediment_deposit__thickness": {
            "dtype": "float",
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "amount of sediment deposited by the plume",
        },
    }

    PLUG_WIDTH = 5.17605
    CONST_ALBERTSON = 0.109

    PLUG_FLOW = 1
    ESTABLISHING_FLOW = 2
    ESTABLISHED_FLOW = 3

    def __init__(
        self,
        grid,
        river_velocity=1.0,
        river_width=1.0,
        river_depth=1.0,
        river_angle=0.0,
        ocean_velocity=0.0,
        river_concentration=1.0,
        river_loc=(0.0, 0.0),
        sediment_removal_rate=1.0,
        sediment_bulk_density=1600.0,
    ):
        """Simulate a hypopycnal sediment plume.

        Parameters
        ----------
        grid : RasterModelGrid
            The solution grid.
        river_velocity: float, optional
            Velocity of the river (m/s).
        river_width: float, optional
            Width of the river (m).
        river_depth: float, optional
            Depth of the river (m).
        river_angle: float, optional
            Direction that river flows into the ocean (rads).
        ocean_velocity: float, optional
            Along-shore current velocity (m/s).
        river_concentration: float, optional
            Concentration of sediment in the river.
        river_loc: tuple of float, optional
            Location of the river mouth.
        """
        self._grid = grid
        self._river = River(
            velocity=river_velocity,
            width=river_width,
            depth=river_depth,
            angle=river_angle,
            loc=river_loc,
        )

        self._river_concentration = river_concentration
        self._centerline = PlumeCenterline(
            river_width,
            river_velocity=river_velocity,
            ocean_velocity=ocean_velocity,
            river_angle=river_angle,
            river_loc=river_loc,
        )

        self._sediment_removal_rate = sediment_removal_rate
        self._sediment_bulk_density = sediment_bulk_density

        self.grid.at_grid["sediment__removal_rate"] = sediment_removal_rate
        self.grid.at_grid["sediment__bulk_density"] = sediment_bulk_density

        for name in (
            "tracer~conservative__mass_concentration",
            "sediment~suspended__mass_concentration",
        ):
            if name not in self.grid.at_node:
                self.grid.add_empty(name, at="node")
        if "sediment_deposit__thickness" not in self.grid.at_node:
            self.grid.add_zeros("sediment_deposit__thickness", at="node")

        super().__init__(grid)

        self._albertson_velocity = self.grid.zeros(at="node")

    @property
    def centerline(self):
        return self._centerline

    @property
    def river(self):
        return self._river

    @property
    def ocean_sed_concentration(self):
        return 0.0

    @property
    def shore_normal(self):
        return self.river.angle

    @property
    def sediment_removal_rate(self):
        return self._sediment_removal_rate

    @property
    def plug_width(self):
        return self.river.width * self.PLUG_WIDTH

    @property
    def established_flow(self):
        try:
            self._established_flow
        except AttributeError:
            self._established_flow = self.where_established_flow()
        return self._established_flow

    @property
    def establishing_flow(self):
        try:
            self._establishing_flow
        except AttributeError:
            self._establishing_flow = self.where_establishing_flow()
        return self._establishing_flow

    @property
    def plug_flow(self):
        try:
            self._plug_flow
        except AttributeError:
            self._plug_flow = self.where_plug_flow()
        return self._plug_flow

    @property
    def distance_to_river(self):
        try:
            self._distance_to_river
        except AttributeError:
            self._distance_to_river = np.sqrt(
                np.power(self.grid.x_of_node - self.river.x0, 2)
                + np.power(self.grid.y_of_node - self.river.y0, 2)
            )
        return self._distance_to_river

    @property
    def concentration(self):
        try:
            self._concentration
        except AttributeError:
            self._concentration = self.calc_concentration()
        return self._concentration

    def calc_concentration(self):
        """Calculate the concentration of a conservative tracer."""
        conc = self.grid.at_node["tracer~conservative__mass_concentration"]
        conc.fill(0.0)

        u_albertson = self._albertson_velocity

        y = self.distance_to_centerline
        x = self.distance_along_centerline

        conc[self.plug_flow] = 1.0
        u_albertson[self.plug_flow] = 1.0

        a = (
            y[self.establishing_flow]
            + 0.5 * SQRT_PI * self.CONST_ALBERTSON * x[self.establishing_flow]
            - 0.5 * self.river.width
        )
        b = np.clip(
            SQRT_TWO * self.CONST_ALBERTSON * x[self.establishing_flow], 0.01, None
        )
        conc[self.establishing_flow] = np.exp(-np.sqrt(a / b))
        u_albertson[self.establishing_flow] = np.exp(-((a / b) ** 2))

        v1 = self.river.width / (
            SQRT_PI * self.CONST_ALBERTSON * x[self.established_flow]
        )
        v2 = y[self.established_flow] / (
            SQRT_TWO * self.CONST_ALBERTSON * x[self.established_flow]
        )
        conc[self.established_flow] = np.sqrt(v1) * np.exp(-np.sqrt(v2))
        u_albertson[self.established_flow] = np.sqrt(v1) * np.exp(-(v2**2))

        return conc

    def calc_sediment_concentration(self, removal_rate):
        # removal_rate /= SECONDS_PER_DAY

        conc_sed = self.grid.at_node["sediment~suspended__mass_concentration"]
        conc_sed.fill(0.0)

        conc_sed = conc_sed.reshape(self.grid.shape)
        concentration = self.concentration.reshape(self.grid.shape)
        u_albertson = self._albertson_velocity.reshape(self.grid.shape)

        uu = 0.2 * self.river.velocity * (1.0 + u_albertson[0] + 3.0 * u_albertson)

        inds = np.where(u_albertson > 0.05)

        conc_sed.fill(self.ocean_sed_concentration)
        conc_sed[inds] = (
            concentration[inds]
            * np.exp(
                -removal_rate
                / SECONDS_PER_DAY
                * self.distance_to_river.reshape(self.grid.shape)[inds]
                / uu[inds]
            )
            + self.ocean_sed_concentration
        )

        return conc_sed

    def calc_deposit_thickness(self, removal_rate, out=None):
        if out is None:
            out = self.grid.zeros(at="node")

        deposit = out.reshape(self.grid.shape)

        # deposit = self.grid.at_node["sediment_deposit__thickness"]
        deposit.fill(0.0)

        bulk_density = self.grid.at_grid["sediment__bulk_density"]
        removal_rate = self.grid.at_grid["sediment__removal_rate"]

        # ocean_cw = 0.0
        # sed_rho = 1600.0
        dl = 0.5 * (self.grid.dx + self.grid.dy)

        conc_sed = self.calc_sediment_concentration(removal_rate)
        u_albertson = (
            self._albertson_velocity.reshape(self.grid.shape) * self.river.velocity
        )
        # deposit = deposit.reshape(self.grid.shape)

        # inds = np.where((conc_sed >= ocean_cw) & (u_albertson > 0.05))
        inds = np.where(
            (u_albertson > 0.05) & (conc_sed > self.ocean_sed_concentration)
        )

        # removal_rate /= SECONDS_PER_DAY
        deposit[inds] = (
            conc_sed[inds]
            * (np.exp(removal_rate / SECONDS_PER_DAY * dl / u_albertson[inds]) - 1.0)
            * (self.river.depth * SECONDS_PER_DAY * u_albertson[inds])
            / (bulk_density * dl)
        )

        return out

    @property
    def xy_at_nearest_centerline(self):
        try:
            self._xy_at_nearest_centerline
        except AttributeError:
            self._xy_at_nearest_centerline = self.calc_nearest_centerline_point()
        return self._xy_at_nearest_centerline

    @property
    def distance_to_centerline(self):
        try:
            self._distance_to_centerline
        except AttributeError:
            self._distance_to_centerline = self.calc_distance_to_centerline()
        return self._distance_to_centerline

    @property
    def distance_along_centerline(self):
        try:
            self._distance_along_centerline
        except AttributeError:
            self._distance_along_centerline = self.calc_distance_along_centerline()
        return self._distance_along_centerline

    @property
    def zones(self):
        try:
            self._zones
        except AttributeError:
            self._zones = self.calc_zones()
        return self._zones

    def calc_nearest_centerline_point(self):
        return self._centerline.nearest_point(
            tuple(zip(self.grid.x_of_node, self.grid.y_of_node))
        )

    def calc_distance_to_centerline(self):
        xy_at_node = tuple(zip(self.grid.x_of_node, self.grid.y_of_node))
        return np.sqrt(
            np.power(self.xy_at_nearest_centerline - xy_at_node, 2).sum(axis=1)
        )

    def calc_distance_along_centerline(self):
        bounds = np.empty((self.grid.number_of_nodes, 2))
        if self.centerline.is_function_of_x():
            bounds[:, 0] = self.river.x0
            bounds[:, 1] = self.xy_at_nearest_centerline[:, 0]
        else:
            bounds[:, 0] = self.river.y0
            bounds[:, 1] = self.xy_at_nearest_centerline[:, 1]
        return self._centerline.path_length(bounds)

    def calc_zones(self):
        zones = np.full(self.grid.number_of_nodes, self.ESTABLISHED_FLOW, dtype=int)
        zones[self.where_plug_flow()] = self.PLUG_FLOW
        zones[self.where_establishing_flow()] = self.ESTABLISHING_FLOW

        return zones

    def where_plug_flow(self):
        lengths = self.distance_along_centerline
        return np.where(
            (lengths < self.plug_width)
            & (
                self.distance_to_centerline
                < (self.river.width * 0.5) * (1 - lengths / (self.plug_width))
            )
        )

    def where_established_flow(self):
        lengths = self.distance_along_centerline
        return np.where(lengths > self.plug_width)

    def where_establishing_flow(self):
        lengths = self.distance_along_centerline
        return np.where(
            (lengths < self.plug_width)
            & (
                self.distance_to_centerline
                >= (self.river.width * 0.5) * (1.0 - lengths / self.plug_width)
            )
        )

    def where_ocean(self):
        v_river = self.unit_vector(np.cos(self.shore_angle), np.sin(self.shore_angle))
        v_point = self.unit_vector(
            self.river.x0 - self.grid.x_of_node, self.river.y0 - self.grid.y_of_node
        )

        return np.where(np.dot(v_river.squeeze(), v_point) <= 0.0)

    def where_land(self):
        v_river = self.unit_vector(np.cos(self.shore_angle), np.sin(self.shore_angle))
        v_point = self.unit_vector(
            self.river.x0 - self.grid.x_of_node, self.river.y0 - self.grid.y_of_node
        )

        return np.where(np.dot(v_river.squeeze(), v_point) > 0.0)

    def is_land(self):
        mask = np.full(self.grid.number_of_nodes, False, dtype=bool)
        mask[self.where_land()] = True
        return mask

    @staticmethod
    def unit_vector(x, y):
        v = np.asarray((x, y)).reshape((2, -1))
        v_abs = np.linalg.norm(v, axis=0)
        return np.divide(v, v_abs, where=v_abs > 0.0, out=np.zeros_like(v))

    def run_one_step(self):
        removal_rate = self.grid.at_grid["sediment__removal_rate"]
        # bulk_density = self.grid.at_grid["sediment__bulk_density"]
        try:
            needs_updating = np.fabs(self._removal_rate - removal_rate) > 1e-12
        except AttributeError:
            needs_updating = True

        if needs_updating:
            deposit = self.calc_deposit_thickness(removal_rate)
            self.grid.at_node["sediment_deposit__thickness"] += deposit

            self._removal_rate = removal_rate
