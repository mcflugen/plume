#! /usr/bin/env python
import argparse
import os
from collections import defaultdict
from io import StringIO

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
    "grid": {"shape": [50, 100], "spacing": [100.0, 100.0], "origin": [0.0, 0.0],},
    "river": {
        "width": 200.0,
        "depth": 1.0,
        "velocity": 1.0,
        "location": [0.0, 0.0],
        "angle": 0.0,
    },
    "sediment": {"removal_rate": 60.0, "bulk_density": 1600.0,},
    "ocean": {"along_shore_velocity": 0.1, "sediment_concentration": 0.0,},
}


LONG_NAME = {
    "c": "tracer~conservative__mass_concentration",
    "cs": "sediment~suspended__mass_concentration",
    "dz": "sediment_deposit__thickness",
}


class RiverPlume:
    def __init__(self, plume="plume.yaml", river="river.csv"):
        self._river = np.loadtxt(river)
        self._plume = Plume()

    def update(self):
        pass


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", nargs="?", help="plume config file")
    parser.add_argument("--output", help="output file")
    parser.add_argument("--verbose", action="store_true", help="be verbose")
    parser.add_argument(
        "--plot", choices=("c", "cs", "dz"), default=None, help="value to plot"
    )
    parser.add_argument(
        "--set", action="append", default=[], help="set plume parameters"
    )

    args = parser.parse_args()

    params = DEFAULT_PARAMS
    if args.file:
        params_from_file = load_params(args.file)
        for group in params.keys():
            params[group].update(params_from_file.get(group, {}))

    params_from_cl = load_params_from_strings(args.set)
    for group in params.keys():
        params[group].update(params_from_cl.get(group, {}))

    if args.verbose:
        print(yaml.dump(params, default_flow_style=False))

    params["river"]["angle"] = np.deg2rad(params["river"]["angle"])

    grid = RasterModelGrid(
        params["grid"]["shape"],
        spacing=params["grid"]["spacing"],
        origin=params["grid"]["origin"],
    )
    plume = Plume(
        grid,
        river_width=params["river"]["width"],
        river_depth=params["river"]["depth"],
        river_velocity=params["river"]["velocity"],
        river_angle=params["river"]["angle"],
        river_loc=params["river"]["location"],
        ocean_velocity=params["ocean"]["along_shore_velocity"],
    )

    plume.grid.at_grid["sediment__removal_rate"] = params["sediment"]["removal_rate"]
    plume.grid.at_grid["sediment__bulk_density"] = params["sediment"]["bulk_density"]
    deposit = plume.calc_deposit_thickness(params["sediment"]["removal_rate"])
    if args.plot:
        imshow_grid(
            plume.grid, plume.grid.at_node[LONG_NAME[args.plot]], at="node", show=True
        )

    if args.output:
        write_raster_netcdf(args.output, plume.grid)


@click.group(chain=True)
@click.version_option()
@click.option(
    "--cd",
    default=".",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="chage to directory, then execute",
)
def plume(cd) -> None:
    """Simulate a sediment-laden hypopycnal plume

    \b
    Examples:

      Create a folder with example input files,

        $ mkdir example && plume setup

      Run a simulation using the examples input files,

        $ cd example && plume run

      Commands can also be chained together. The folling will setup
      a simulation, run it, and then plot elevations.

        $ mkdir example && plume setup run plot concentration
    """
    os.chdir(cd)


@plume.command()
@click.option("-v", "--verbose", is_flag=True, help="Emit status messages to stderr.")
@click.option("--dry-run", is_flag=True, help="Do not actually run the model")
def run(dry_run: bool, verbose: bool) -> None:
    """Run a simulation."""
    run_dir = pathlib.Path.cwd()
    config_file = run_dir / "plume.yaml"

    message = []
    if not config_file.is_file():
        message.append("missing plume configuration file: {0}".format(config_file))
    if (run_dir / "output").exists():
        message.append(
            "plume output directory already exists: {0}".format(run_dir / "output")
        )
    if message:
        err(os.linesep.join(message))
        raise click.Abort(os.linesep.join(message))

    params = load_params(config_file)
    plume = Plume(**params)

    grid = RasterModelGrid(
        params["grid"]["shape"],
        xy_spacing=params["grid"]["xy_spacing"],
        xy_of_origin=params["grid"]["xy_of_origin"],
    )

    if dry_run:
        out("Nothing to do. 😴")
    else:
        spacing = 10
        n_days = params["days"]  # + years * 365
        n_steps = n_days

        if n_steps == 0:
            out("Nothing to do (days == 0). 😴")
        else:
            with click.progressbar(
                range(n_steps),
                label=" ".join(["🚀", str(run_dir)]),
                item_show_func=lambda step: "day {0} of {1}".format(
                    int(0 if step is None else step), n_days
                ),
            ) as bar:
                for step in bar:
                    plume.update()
                    # output.update(1)

            out("💥 Finished! 💥")
            out("Output written to {0}".format(output.prefix))


@plume.command()
@click.argument("infile", type=click.Choice(["plume.yaml", "river.csv"]))
def generate(infile: str) -> None:
    """Show example input files."""
    print(_format_input_file(infile))


@plume.command()
def setup() -> None:
    """Setup a folder of input files for a simulation."""
    files = [pathlib.Path("plume.yaml"), pathlib.Path("river.csv")]

    existing_files = [name for name in files if name.exists()]
    if existing_files:
        for name in existing_files:
            err(
                f"{name}: File exists. Either remove and then rerun or setup in a different folder",
            )
    else:
        for fname in files:
            with open(fname, "w") as fp:
                print(_format_input_file(fname), file=fp)

    if existing_files:
        raise click.Abort()


def _format_input_file(infile: str) -> str:
    def as_csv(data, header=None):
        with StringIO() as fp:
            np.savetxt(fp, data, header=header, delimiter=",", fmt="%.1f")
            contents = fp.getvalue()
        return contents

    contents = {
        # "plume.yaml": yaml.dump(params, default_flow_style=False),
        "river.csv": as_csv(
            [[0, 1.0, 200.0, 1.0, 0.0, 0.0, 0.0, 0.15]],
            header="Time [d], Velocity [m/s], Width [m], Depth [m], Mouth X [m], Mouth Y [m], Angle [deg], Ocean Velocity [m/s]",
        ),
    }
    return contents[infile]
