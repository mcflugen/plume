#! /usr/bin/env python
from __future__ import print_function

import argparse
from collections import defaultdict

import numpy as np
import yaml

from landlab import RasterModelGrid, load_params
from landlab.plot import imshow_grid
from landlab.io.netcdf import write_raster_netcdf

from .plume import Plume


DEFAULT_PARAMS = {
    'grid': {
        'shape': [50, 100],
        'spacing': [100., 100.],
        'origin': [0., 0.],
    },
    'river': {
        'width': 200.,
        'depth': 1.,
        'velocity': 1.,
        'location': [0., 0.],
        'angle': 0.,
    },
    'sediment': {
        'removal_rate': 60.,
        'bulk_density': 1600.,
    },
    'ocean': {
        'along_shore_velocity': .1,
        'sediment_concentration': 0.,
    }
}


LONG_NAME = {
    'c': 'tracer~conservative__mass_concentration',
    'cs': 'sediment~suspended__mass_concentration',
    'dz': 'sediment_deposit__thickness',
}


def load_params_from_strings(values):
    params = defaultdict(dict)

    for param in values:
        group_dot_name, value = param.split('=')
        value = yaml.load(value)
        try:
            group, name = group_dot_name.split('.')
        except ValueError:
            name = group_dot_name
            params[name] = value
        else:
            params[group][name] = value

    return params


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file', nargs='?', help='plume config file')
    parser.add_argument('--output', help='output file')
    parser.add_argument('--verbose', action='store_true',
                        help='be verbose')
    parser.add_argument('--plot', choices=('c', 'cs', 'dz'), default=None,
                        help='value to plot')
    parser.add_argument('--set', action='append', default=[],
                        help='set plume parameters')

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

    params['river']['angle'] = np.deg2rad(params['river']['angle'])

    grid = RasterModelGrid(params['grid']['shape'],
                           spacing=params['grid']['spacing'],
                           origin=params['grid']['origin'])
    plume = Plume(grid,
                  river_width=params['river']['width'],
                  river_depth=params['river']['depth'],
                  river_velocity=params['river']['velocity'],
                  river_angle=params['river']['angle'],
                  river_loc=params['river']['location'],
                  ocean_velocity=params['ocean']['along_shore_velocity'],
                  )

    plume.grid.at_grid['sediment__removal_rate'] = params['sediment']['removal_rate']
    plume.grid.at_grid['sediment__bulk_density'] = params['sediment']['bulk_density']
    deposit = plume.calc_deposit_thickness(params['sediment']['removal_rate'])
    if args.plot:
        imshow_grid(plume.grid, plume.grid.at_node[LONG_NAME[args.plot]],
                    at='node', show=True)

    if args.output:
        write_raster_netcdf(args.output, plume.grid)
