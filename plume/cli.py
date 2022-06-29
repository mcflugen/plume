#! /usr/bin/env python
import os
import pathlib
import textwrap
import warnings
from collections import defaultdict
from functools import partial
from io import StringIO
from typing import Optional, TextIO

import rich_click as click
import numpy as np
import yaml
from landlab import RasterModelGrid
from landlab.io.netcdf import write_raster_netcdf
from packaging.version import parse as parse_version
from scipy import interpolate

from ._version import __version__
from .plume import Plume

import tomlkit as tomllib

click.rich_click.ERRORS_SUGGESTION = (
    "Try running the '--help' flag for more information."
)
click.rich_click.ERRORS_EPILOGUE = (
    "To find out more, visit https://github.com/csdms/bmi-wavewatch3"
)
click.rich_click.STYLE_ERRORS_SUGGESTION = "yellow italic"
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = False
click.rich_click.SHOW_METAVARS_COLUMN = True
click.rich_click.USE_MARKDOWN = True
out = partial(click.secho, bold=True, err=True)
err = partial(click.secho, fg="red", err=True)


# class RiverParameters:
#     name = "river"
#     width = Parameter(
#         "width", default=50.0, valid=Range(0.0, max=None), help="river width", units="m",
#     )
#     filepath = Parameter(
#         "filepath", default="river.csv", valid=Path(exists=True, file_okay=True, dir_okay=False), help="river file",
#     )
#
#
# class SedimentParameters:
#     name = "sediment"
#     removal_rate = Parameter(
#     )
#     bulk_density = Parameter(
#     )
#
#
# class OceanParameters:
#     name = "ocean"
#     filepath = Parameter(
#         "filepath", default="ocean.csv", valid=CsvFile(exists=True, nrows="+", ncols=3), help="river file",
#     )
#
#
# class GridParameters:
#     name = "grid"
#     shape = Parameter("shape", [500, 500], valid=Array(length=2, dtype=int), help="number of grid rows and columns")
#     shape = Parameter("shape", [500, 500], valid=Length(2), help="number of grid rows and columns")
#     xy_spacing = Parameter("xy_spacing", [100.0, 100.0], valid=Length(2), help="spacing of grid columns and rows")
#     xy_of_lower_left = Parameter("xy_of_lower_left", [0.0, 0.0], valid=Length(2), help="coordinates of lower-left node of grid")


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
        # "_version": __version__,
        "grid": {
            "shape": [500, 500],
            "xy_spacing": [100.0, 100.0],
            "xy_of_lower_left": [0.0, 0.0],
        },
        "river": {
            # "filepath": "river.csv",
            "width": 50.0,
            "depth": 5.0,
            "velocity": 1.5,
            "location": [0.0, 25000.0],
            "angle": 0.0,
        },
        "sediment": {"removal_rate": 60.0, "bulk_density": 1600.0},
        "ocean": {
            # "filepath": "ocean.csv",
            "along_shore_velocity": 0.1,
            "sediment_concentration": 0.0,
        },
        "output": {"filepath": "plume.nc"},
    }
    if file is not None:
        conf.update(yaml.safe_load(file))

    conf.setdefault("_version", __version__)

    this_version = parse_version(__version__)
    that_version = parse_version(conf["_version"])
    if this_version.major != that_version.major:
        warnings.warn(
            "possible version mismatch. file is version v{}, but you are using plume v{}".format(
                that_version, this_version
            )
        )

    return conf


def _contents_of_input_file(infile: str) -> str:
    params = load_config()

    def as_csv(data, header=None):
        with StringIO() as fp:
            np.savetxt(fp, data, header=header, delimiter=",", fmt="%.1f")
            contents = fp.getvalue()
        return contents

    contents = {
        "plume.toml": tomllib.dumps(dict(plume=params)),
        # "river.csv": as_csv(
        #     [[0.0, 50.0, 5.0, 1.5]],
        #     header=os.linesep.join(
        #         [
        #             f"version: {__version__}",
        #             "Time [d], Width [m], Depth [m], Velocity [m/s]",
        #         ]
        #     ),
        # ),
        # "ocean.csv": as_csv(
        #     [[0.0, 0.1, 0.0]],
        #     header=os.linesep.join(
        #         [
        #             f"version: {__version__}",
        #             "Time [d], Along-shore velocity [m/s], Sediment Concentration [-]",
        #         ]
        #     ),
        # ),
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


@click.group(chain=True)
@click.version_option(package_name="sed-plume")
@click.option(
    "--cd",
    default=".",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="chage to directory, then execute",
)
@click.option(
    "-s",
    "--silent",
    is_flag=True,
    help="Suppress status status messages, including the progress bar.",
)
@click.option(
    "-v", "--verbose", is_flag=True, help="Also emit status messages to stderr."
)
def plume(cd: str, silent: bool, verbose: bool) -> None:
    """Simulate a hypopycnal plume.

    ## Examples

    Create a folder with example input files,

    ```bash
    $ plume setup example
    ```

    Run a simulation using the examples input files,

    ```bash
    $ plume run example
    ```
    """
    if not silent and verbose:
        out(f"run_dir = {cd!r}")
    os.chdir(cd)


class RiverTimeSeries:
    def __init__(self, filepath, kind="linear"):
        data = np.loadtxt(filepath, delimiter=",", comments="#").reshape((-1, 4))
        if len(data) == 1:
            data = np.vstack([data, data])
            data[1, 0] = data[0, 0] + 1

        self._river = interpolate.interp1d(
            data[:, 0],
            data[:, 1:],
            kind=kind,
            axis=0,
            copy=True,
            assume_sorted=True,
            bounds_error=False,
            fill_value=(data[0, 1:], data[-1, 1:]),
        )
        self._time = 0.0
        self._width, self._depth, self._velocity = self._river(self._time)

    @property
    def velocity(self):
        return self._velocity

    @property
    def width(self):
        return self._width

    @property
    def depth(self):
        return self._depth

    def update(self):
        self._time += 1.0
        self._width, self._depth, self._velocity = self._river(self._time)


class OceanTimeSeries:
    def __init__(self, filepath, kind="linear"):
        data = np.loadtxt(filepath, delimiter=",", comments="#").reshape((-1, 3))
        if len(data) == 1:
            data = np.vstack([data, data])
            data[1, 0] = data[0, 0] + 1

        self._ocean = interpolate.interp1d(
            data[:, 0],
            data[:, 1:],
            kind=kind,
            copy=True,
            assume_sorted=True,
            bounds_error=False,
            fill_value=(data[0, 1:], data[-1, 1:]),
        )
        self._time = 0.0
        self._velocity, _ = self._ocean(self._time)

    @property
    def velocity(self):
        return self._velocity

    def update(self):
        self._time += 1.0
        self._velocity, _ = self._ocean(self._time)


@plume.command()
@click.option("--dry-run", is_flag=True, help="Do not actually run the model")
@click.pass_context
def run(ctx, dry_run: bool) -> None:
    """Run the plume simulation on a set of input files.

    ## Examples

    ```bash
    $ mkdir example && cd example
    $ plume setup
    $ plume run
    """
    verbose = ctx.parent.params["verbose"]
    silent = ctx.parent.params["silent"]

    # params = load_params("plume.toml")
    with open("plume.toml") as fp:
        params = tomllib.load(fp)

    if verbose and not silent:
        out("plume.toml = |")
        out(textwrap.indent(tomllib.dumps(dict(plume=params)).strip(), prefix="  "))

    params = params["plume"]

    if dry_run:
        out("Nothing to do. 😴")
    else:
        output_file = pathlib.Path(params["output"]["filepath"])
        if output_file.exists():
            err(f"{output_file}: output files exists")
            raise click.Abort()

        grid = RasterModelGrid(
            params["grid"]["shape"],
            xy_spacing=params["grid"]["xy_spacing"],
            xy_of_lower_left=params["grid"]["xy_of_lower_left"],
        )

        grid.add_zeros("sediment_deposit_thickenss", at="node")  # .reshape(grid.shape)
        grid.at_grid["sediment__removal_rate"] = params["sediment"]["removal_rate"]
        grid.at_grid["sediment__bulk_density"] = params["sediment"]["bulk_density"]

        # river = RiverTimeSeries(params["river"]["filepath"])
        # ocean = OceanTimeSeries(params["ocean"]["filepath"])

        params["river"]["angle"] = np.deg2rad(params["river"]["angle"])

        # n_days = 2
        # for day in range(n_days):

        plume = Plume(
            grid,
            river_width=params["river"]["width"],
            river_depth=params["river"]["depth"],
            river_velocity=params["river"]["velocity"],
            river_angle=params["river"]["angle"],
            river_loc=params["river"]["location"],
            ocean_velocity=params["ocean"]["along_shore_velocity"],
        )
        plume.run_one_step()

        write_raster_netcdf(
            output_file,
            plume.grid,
            time=0.0,
            append=True,
            attrs={"_version": __version__},
        )

        out("💥 Finished! 💥")
        out(f"Output written to {output_file}")

        print(output_file)


@plume.command()
@click.argument("infile", type=click.Choice(sorted(["plume.toml"])))
def generate(infile: str) -> None:
    """Print example input files.

    ## Examples

    To see an example *ocean.csv* input file,
    ```bash
    $ plume generate ocean.csv
    ```
    """
    print(_contents_of_input_file(infile))


@plume.command()
@click.pass_context
def setup(ctx) -> None:
    """Setup a folder of input files for a simulation.

    *plume setup* creates a set of example input files for the *plume* program.
    This set of input files can then be used with the *plume run* command.

    > **NOTE**: *plume setup* will not overwrite any of your files. If you want to overwrite
    existing files, you will have to remove them yourself and then re-run
    *plume setup*.

    ## Examples

    ```bash
    $ mkdir _example
    $ plume --cd=_example setup
    ```
    """
    verbose = ctx.parent.params["verbose"]
    silent = ctx.parent.params["silent"]

    folder = pathlib.Path(".")

    files = [pathlib.Path(fname) for fname in ["plume.toml"]]

    existing_files = [folder / name for name in files if (folder / name).exists()]
    if existing_files:
        for name in existing_files:
            err(
                f"{name}: File exists. Either remove and then rerun or choose a different destination folder",
            )
    else:
        for fname in files:
            with open(folder / fname, "w") as fp:
                print(_contents_of_input_file(str(fname)), file=fp)
        if verbose and not silent:
            out(f"files = {[str(f) for f in files]!r}")
        print(str(folder.absolute()))

    if existing_files:
        raise click.Abort()
