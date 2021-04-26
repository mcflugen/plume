#! /usr/bin/env python
from functools import lru_cache

import numpy as np
import scipy
from landlab import Component, FieldError

from .centerline import PlumeCenterline
from .river import River

SQRT_PI = np.sqrt(np.pi)
SQRT_TWO = np.sqrt(2.0)
SECONDS_PER_DAY = 60 * 60 * 24.0


from collections import UserDict, namedtuple

FieldInfo = namedtuple(
    "FieldInfo", ["dtype", "intent", "optional", "units", "mapping", "doc"]
)


class ChoiceError(Exception):
    def __init__(self, value, choices):
        self.value = value
        self.choices = list(choices)

    def __str__(self):
        return "{value} is not one of {choices}".format(
            value=repr(self.value),
            choices=", ".join([repr(choice) for choice in self.choices]),
        )


class FieldInfo(UserDict):
    # def __init__(self, name, dtype="float", intent="inout", optional=False, units="", mapping="node", doc=""):
    def __init__(self, *args, **kwds):
        info = dict(*args, **kwds)

        self.data = {
            "name": info.get("name", None),
            "dtype": self._validate_dtype(info.get("dtype", float)),
            "intent": self._validate_intent(info.get("intent", "in")),
            "optional": self._validate_optional(info.get("optional", False)),
            "units": info.get("units", ""),
            "mapping": self._validate_mapping(info.get("mapping", "node")),
            "doc": info.get("doc", ""),
        }

        if info.keys() - self.data.keys():
            raise ValueError(
                "unknown keywords ({0})".format(
                    ", ".join([repr(key) for key in info.keys() - self.data.keys()])
                )
            )

    @property
    def name(self):
        return self.data["name"]

    @property
    def intent(self):
        return self.data["intent"]

    def _validate_intent(self, intent):
        valid_intents = ["in", "out", "inout"]
        if intent not in valid_intents:
            raise ChoiceError(intent, valid_intents)
        return intent

    @property
    def optional(self):
        return self.data["optional"]

    def _validate_optional(self, optional):
        if optional not in (True, False):
            raise ChoiceError(optional, (True, False))
        return optional

    @property
    def dtype(self):
        return self.data["dtype"]

    def _validate_dtype(self, dtype):
        try:
            return str(np.dtype(dtype))
        except ValueError:
            raise ValueError("invalid value for 'dtype'")

    @property
    def units(self):
        return self.data["units"]

    @property
    def mapping(self):
        return self.data["mapping"]

    def _validate_mapping(self, mapping):
        valid_mappings = ["node", "link", "patch", "corner", "face", "cell", "grid"]
        if mapping not in valid_mappings:
            raise ChoiceError(mapping, valid_mappings)
        return mapping

    @property
    def doc(self):
        return self.data["doc"]

    def __setitem__(self, key, value):
        raise NotImplementedError("FieldInfo is read-only")

    def __delitem__(self, key, value):
        raise NotImplementedError("FieldInfo is read-only")


class Plume(Component):

    _name = "Plume"

    _input_var_names = (
        "sediment__removal_rate",
        "sediment__bulk_density",
        "channel_exit_water_flow__speed",
        "channel_exit__width",
    )

    _output_var_names = (
        "tracer~conservative__mass_concentration",
        "sediment~suspended__mass_concentration",
        "sediment_deposit__thickness",
    )

    _var_units = {
        "tracer~conservative__mass_concentration": "kg / m^3",
        "sediment~suspended__mass_concentration": "kg / m^3",
        "sediment_deposit__thickness": "m",
        "sediment__removal_rate": "1 / d",
        "sediment__bulk_density": "kg / m^3",
        "channel_exit_water_flow__speed": "m / s",
        "channel_exit__width": "m",
    }

    _var_mapping = {
        "tracer~conservative__mass_concentration": "node",
        "sediment~suspended__mass_concentration": "node",
        "sediment_deposit__thickness": "node",
        "sediment__removal_rate": "grid",
        "sediment__bulk_density": "grid",
        "channel_exit_water_flow__speed": "grid",
        "channel_exit__width": "grid",
    }

    _var_doc = {
        "tracer~conservative__mass_concentration": "concentration of a conservative tracer",
        "sediment~suspended__mass_concentration": "concentration of suspended sediment in the plume",
        "sediment_deposit__thickness": "amount of sediment deposited by the plume",
        "sediment__removal_rate": "removal rate of sediment carried by the plume",
        "sediment__bulk_density": "bulk density of sediment deposited by the plume",
    }

    _info = {
        "sediment__removal_rate": FieldInfo(
            dtype=float,
            intent="in",
            optional=False,
            units="1 / d",
            mapping="grid",
            doc="removal rate of sediment carried by the plume",
        ),
        "sediment__bulk_density": FieldInfo(
            dtype=float,
            intent="in",
            optional=False,
            units="kg / m^3",
            mapping="grid",
            doc="bulk density of sediment deposited by the plume",
        ),
        "channel_exit_water_flow__speed": FieldInfo(
            dtype=float,
            intent="in",
            optional=False,
            units="m / s",
            mapping="grid",
            doc="flow velocity at the river mouth",
        ),
        "channel_exit__width": FieldInfo(
            dtype=float,
            intent="in",
            optional=False,
            units="m",
            mapping="grid",
            doc="channel width at the river mouth",
        ),
        "tracer~conservative__mass_concentration": FieldInfo(
            dtype=float,
            intent="out",
            units="kg / m^3",
            mapping="node",
            doc="concentration of a conservative tracer",
        ),
        "sediment~suspended__mass_concentration": FieldInfo(
            dtype=float,
            intent="out",
            units="kg / m^3",
            mapping="node",
            doc="concentration of suspended sediment in the plume",
        ),
        "sediment_deposit__thickness": FieldInfo(
            dtype=float,
            intent="out",
            units="m",
            mapping="node",
            doc="amount of sediment deposited by the plume",
        ),
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
        # self._river_width = river_width
        # self._river_velocity = river_velocity
        self._ocean_velocity = ocean_velocity
        self._river_angle = river_angle
        self._river_loc = river_loc

        # self._sediment_removal_rate = sediment_removal_rate
        # self._sediment_bulk_density = sediment_bulk_density

        self.grid.at_grid["sediment__removal_rate"] = sediment_removal_rate
        self.grid.at_grid["sediment__bulk_density"] = sediment_bulk_density
        self.grid.at_grid["channel_exit_water_flow__speed"] = river_velocity
        self.grid.at_grid["channel_exit__width"] = river_width

        super(Plume, self).__init__(grid)

        for name in (
            "sediment_deposit__thickness",
            "sediment~suspended__mass_concentration",
            "tracer~conservative__mass_concentration",
        ):
            self.grid.add_zeros(
                name,
                at=self._info[name].mapping,
                units=self._info[name].units,
                clobber=False,
            )

        self._albertson_velocity = self.grid.zeros(at="node")

    @property
    def sediment_removal_rate(self):
        return float(self.grid.at_grid["sediment__removal_rate"])

    @property
    def sediment_bulk_density(self):
        return float(self.grid.at_grid["sediment__bulk_density"])

    @property
    def river_velocity(self):
        return float(self.grid.at_grid["channel_exit_water_flow__speed"])

    @river_velocity.setter
    def river_velocity(self, new_val):
        if new_val != self.river_velocity:
            self.clear_cache()
            self.grid.at_grid["channel_exit_water_flow__speed"] = new_val

    @property
    def river_width(self):
        return float(self.grid.at_grid["channel_exit__width"])

    @river_width.setter
    def river_width(self, new_val):
        if new_val != self.river_width:
            self.clear_cache()
            self.grid.at_grid["channel_exit__width"] = new_val

    def clear_cache(self):
        attrs = (
            "_centerline",
            "_concentration",
            "_established_flow",
            "_establishing_flow",
            "_plug_flow",
            "_distance_to_river",
            "_xy_at_nearest_centerline",
            "_distance_to_centerline",
            "_distance_along_centerline",
            "_zones",
        )
        for attr in attrs:
            try:
                del self.__dict__[attr]
            except KeyError:
                pass

    @property
    def centerline(self):
        try:
            self._centerline
        except AttributeError:
            self._centerline = PlumeCenterline(
                self.river_width,
                river_velocity=self.river_velocity,
                ocean_velocity=self._ocean_velocity,
                river_angle=self._river_angle,
                river_loc=self._river_loc,
            )
        finally:
            return self._centerline

    def _calc_centerline(
        self,
        river_width,
        river_velocity,
        river_angle,
        xy_of_river,
        ocean_velocity,
    ):
        return PlumeCenterline(
            river_width,
            river_velocity=river_velocity,
            ocean_velocity=ocean_velocity,
            river_angle=river_angle,
            river_loc=xy_of_river,
        )

    @property
    def ocean_sed_concentration(self):
        return 0.0

    @property
    def shore_normal(self):
        return self._river.angle

    @property
    def plug_width(self):
        return self.river_width * self.PLUG_WIDTH

    @property
    def established_flow(self):
        try:
            self._established_flow
        except AttributeError:
            self._established_flow = self.where_established_flow()
        finally:
            return self._established_flow

    @property
    def establishing_flow(self):
        try:
            self._establishing_flow
        except AttributeError:
            self._establishing_flow = self.where_establishing_flow()
        finally:
            return self._establishing_flow

    @property
    def plug_flow(self):
        try:
            self._plug_flow
        except AttributeError:
            self._plug_flow = self.where_plug_flow()
        finally:
            return self._plug_flow

    @property
    def distance_to_river(self):
        try:
            self._distance_to_river
        except AttributeError:
            self._distance_to_river = np.sqrt(
                np.power(self.grid.x_of_node - self._river.x0, 2)
                + np.power(self.grid.y_of_node - self._river.y0, 2)
            )
        finally:
            return self._distance_to_river

    @property
    def concentration(self):
        try:
            self._concentration
        except AttributeError:
            self._concentration = self.calc_concentration()
        finally:
            return self._concentration

    def calc_concentration(self):
        """Calculate the concentration of a conservative tracer."""
        # depends on:
        # *  river_width
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
            - 0.5 * self.river_width
        )
        b = np.clip(
            SQRT_TWO * self.CONST_ALBERTSON * x[self.establishing_flow], 0.01, None
        )
        conc[self.establishing_flow] = np.exp(-np.sqrt(a / b))
        u_albertson[self.establishing_flow] = np.exp(-((a / b) ** 2))

        v1 = self.river_width / (
            SQRT_PI * self.CONST_ALBERTSON * x[self.established_flow]
        )
        v2 = y[self.established_flow] / (
            SQRT_TWO * self.CONST_ALBERTSON * x[self.established_flow]
        )
        conc[self.established_flow] = np.sqrt(v1) * np.exp(-np.sqrt(v2))
        u_albertson[self.established_flow] = np.sqrt(v1) * np.exp(-(v2 ** 2))

        return conc

    def calc_sediment_concentration(self, removal_rate):
        # depends on:
        # *  removal_rate
        # *  river_velocity
        # *  ocean_sed_concentration
        # removal_rate /= SECONDS_PER_DAY

        conc_sed = self.grid.at_node["sediment~suspended__mass_concentration"]
        conc_sed.fill(0.0)

        conc_sed = conc_sed.reshape(self.grid.shape)
        concentration = self.calc_concentration().reshape(self.grid.shape)
        # concentration = self.concentration.reshape(self.grid.shape)
        u_albertson = self._albertson_velocity.reshape(self.grid.shape)

        uu = 0.2 * self.river_velocity * (1.0 + u_albertson[0] + 3.0 * u_albertson)

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

    def calc_deposit_thickness(self, removal_rate):
        # depends on:
        # *  bulk_density
        # *  removal_rate
        deposit = self.grid.at_node["sediment_deposit__thickness"]
        deposit.fill(0.0)

        bulk_density = self.grid.at_grid["sediment__bulk_density"]
        removal_rate = self.grid.at_grid["sediment__removal_rate"]

        ocean_cw = 0.0
        sed_rho = 1600.0
        dl = 0.5 * (self.grid.dx + self.grid.dy)

        conc_sed = self.calc_sediment_concentration(removal_rate)
        u_albertson = (
            self._albertson_velocity.reshape(self.grid.shape) * self.river_velocity
        )
        deposit = deposit.reshape(self.grid.shape)

        # inds = np.where((conc_sed >= ocean_cw) & (u_albertson > 0.05))
        inds = np.where(
            (u_albertson > 0.05) & (conc_sed > self.ocean_sed_concentration)
        )

        # removal_rate /= SECONDS_PER_DAY
        deposit[inds] = (
            conc_sed[inds]
            * (np.exp(removal_rate / SECONDS_PER_DAY * dl / u_albertson[inds]) - 1.0)
            * (self._river.depth * SECONDS_PER_DAY * u_albertson[inds])
            / (bulk_density * dl)
        )

        return deposit

    @property
    def xy_at_nearest_centerline(self):
        try:
            self._xy_at_nearest_centerline
        except AttributeError:
            self._xy_at_nearest_centerline = self.calc_nearest_centerline_point()
        finally:
            return self._xy_at_nearest_centerline

    @property
    def distance_to_centerline(self):
        try:
            self._distance_to_centerline
        except AttributeError:
            self._distance_to_centerline = self.calc_distance_to_centerline()
        finally:
            return self._distance_to_centerline

    @property
    def distance_along_centerline(self):
        try:
            self._distance_along_centerline
        except AttributeError:
            self._distance_along_centerline = self.calc_distance_along_centerline()
        finally:
            return self._distance_along_centerline

    @property
    def zones(self):
        try:
            self._zones
        except AttributeError:
            self._zones = self.calc_zones()
        finally:
            return self._zones

    def calc_nearest_centerline_point(self):
        return self.centerline.nearest_point(self.grid.xy_of_node)

    def calc_distance_to_centerline(self):
        xy_at_node = self.grid.xy_of_node
        return np.sqrt(
            np.power(self.xy_at_nearest_centerline - xy_at_node, 2).sum(axis=1)
        )

    def calc_distance_along_centerline(self):
        bounds = np.empty((self.grid.number_of_nodes, 2))
        if self.centerline.is_function_of_x():
            bounds[:, 0] = self._river.x0
            bounds[:, 1] = self.xy_at_nearest_centerline[:, 0]
        else:
            bounds[:, 0] = self._river.y0
            bounds[:, 1] = self.xy_at_nearest_centerline[:, 1]
        return self.centerline.path_length(bounds)

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
                < (self.river_width * 0.5) * (1 - lengths / (self.plug_width))
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
                >= (self.river_width * 0.5) * (1.0 - lengths / self.plug_width)
            )
        )

    def where_ocean(self):
        v_river = self.unit_vector(np.cos(self.shore_angle), np.sin(self.shore_angle))
        v_point = self.unit_vector(
            self._river.x0 - self.grid.x_of_node, self._river.y0 - self.grid.y_of_node
        )

        return np.where(np.dot(v_river.squeeze(), v_point) <= 0.0)

    def where_land(self):
        v_river = self.unit_vector(np.cos(self.shore_angle), np.sin(self.shore_angle))
        v_point = self.unit_vector(
            self._river.x0 - self.grid.x_of_node, self._river.y0 - self.grid.y_of_node
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
        self._run_one_step(
            self.sediment_removal_rate,
            self.sediment_bulk_density,
            self.river_velocity,
            self.river_width,
        )

    @lru_cache(1)
    def _run_one_step(self, *args):
        print("recalculating...")
        self.clear_cache()
        self.calc_deposit_thickness(self.sediment_removal_rate)

    def update(self):
        return self.run_one_step()
