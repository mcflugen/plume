#! /usr/bin/env python
import os
import pathlib
import sys
from collections import defaultdict
from functools import partial
from io import StringIO
from typing import Optional, TextIO

import click
import numpy as np
import yaml
from landlab import RasterModelGrid, load_params
from landlab.io.netcdf import write_raster_netcdf
from landlab.plot import imshow_grid

from .plume import Plume

out = partial(click.secho, bold=True, err=True)
err = partial(click.secho, fg="red", err=True)

DEFAULT_PARAMS = {
    "grid": {
        "shape": [50, 100],
        "xy_spacing": [100.0, 100.0],
        "xy_of_lower_left": [0.0, 0.0],
    },
    "river": {
        "filepath": "river.csv",
        "width": 200.0,
        "depth": 1.0,
        "velocity": 1.0,
        "location": [0.0, 0.0],
        "angle": 0.0,
    },
    "sediment": {"removal_rate": 60.0, "bulk_density": 1600.0,},
    "ocean": {"along_shore_velocity": 0.1, "sediment_concentration": 0.0,},
    "output": {"filepath": "plume.nc"},
}


LONG_NAME = {
    "c": "tracer~conservative__mass_concentration",
    "cs": "sediment~suspended__mass_concentration",
    "dz": "sediment_deposit__thickness",
}


def load_config(file: Optional[TextIO] = None):
    """Load plume config file.

    Parameters
    ----------
    fname : file-like, optional
        Opened config file or ``None``. If ``None``, return default
        values.

    Returns
    -------
    dict
        Config parameters.
    """
    conf = {
        "grid": {
            "shape": [50, 100],
            "xy_spacing": [100.0, 100.0],
            "xy_of_lower_left": [0.0, 0.0],
        },
        "river": {
            "filepath": "river.csv",
            "width": 200.0,
            "depth": 1.0,
            "velocity": 1.0,
            "location": [0.0, 0.0],
            "angle": 0.0,
        },
        "sediment": {"removal_rate": 60.0, "bulk_density": 1600.0,},
        "ocean": {"along_shore_velocity": 0.1, "sediment_concentration": 0.0,},
        "output": {"filepath": "plume.nc"},
    }
    if file is not None:
        conf.update(yaml.safe_load(file))
    return conf


def _contents_of_input_file(infile: str) -> str:
    params = load_config()

    def as_csv(data, header=None):
        with StringIO() as fp:
            np.savetxt(fp, data, header=header, delimiter=",", fmt="%.1f")
            contents = fp.getvalue()
        return contents

    contents = {
        "plume": yaml.dump(params, default_flow_style=False),
        "river": as_csv(
            [[0.0, 200.0, 1.0, 1.0]],
            header="Time [d], Width [m], Depth [m], Velocity [m/s]",
        ),
    }

    return contents[infile]


def load_params_from_strings(values):
    params = defaultdict(dict)

    for param in values:
        group_dot_name, value = param.split("=")
        value = yaml.load(value)
        try:
            group, name = group_dot_name.split(".")
        except ValueError:
            name = group_dot_name
            params[name] = value
        else:
            params[group][name] = value

    return params


@click.group()
@click.version_option()
def plume() -> None:
    """Simulate a hypopycnal plume.

    \b
    Examples:

      Create a folder with example input files,

        $ plume setup example

      Run a simulation using the examples input files,

        $ plume run example
    """
    pass  # pragma: no cover


@plume.command()
@click.option("-v", "--verbose", is_flag=True, help="Emit status messages to stderr.")
@click.option("--dry-run", is_flag=True, help="Do not actually run the model")
@click.argument(
    "run_dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
)
def run(run_dir: str, dry_run: bool, verbose: bool) -> None:
    os.chdir(run_dir)

    params = load_params("plume.yaml")
    if verbose:
        out(yaml.dump(params, default_flow_style=False))

    if dry_run:
        out("Nothing to do. 😴")
    else:
        grid = RasterModelGrid(
            params["grid"]["shape"],
            xy_spacing=params["grid"]["xy_spacing"],
            xy_of_lower_left=params["grid"]["xy_of_lower_left"],
        )

        river = np.loadtxt(
            params["river"]["filepath"], delimiter=",", comments="#"
        ).reshape((-1, 4))

        for day, width, depth, velocity in river:
            params["river"]["angle"] = np.deg2rad(params["river"]["angle"])

            plume = Plume(
                grid,
                river_width=width,
                river_depth=depth,
                river_velocity=velocity,
                river_angle=params["river"]["angle"],
                river_loc=params["river"]["location"],
                ocean_velocity=params["ocean"]["along_shore_velocity"],
            )

            plume.grid.at_grid["sediment__removal_rate"] = params["sediment"][
                "removal_rate"
            ]
            plume.grid.at_grid["sediment__bulk_density"] = params["sediment"][
                "bulk_density"
            ]
            deposit = plume.calc_deposit_thickness(params["sediment"]["removal_rate"])

        write_raster_netcdf(params["output"]["filepath"], plume.grid)

        out("💥 Finished! 💥")
        out("Output written to {0}".format(params["output"]["filepath"]))

    sys.exit(0)


@plume.command()
@click.argument("infile", type=click.Choice(["plume", "river"]))
def show(infile: str) -> None:
    """Show example input files."""
    print(_contents_of_input_file(infile))


@plume.command()
@click.argument(
    "dest", type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
)
def setup(dest: str) -> None:
    """Setup a folder of input files for a simulation."""
    folder = pathlib.Path(dest)

    files = [pathlib.Path(fname) for fname in ["plume.yaml", "river.csv"]]

    existing_files = [folder / name for name in files if (folder / name).exists()]
    if existing_files:
        for name in existing_files:
            err(
                f"{name}: File exists. Either remove and then rerun or choose a different destination folder",
            )
    else:
        for fname in files:
            with open(folder / fname, "w") as fp:
                print(_contents_of_input_file(fname.stem), file=fp)
        print(str(folder))

    sys.exit(len(existing_files))
