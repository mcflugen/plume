.. image:: https://travis-ci.org/mcflugen/plume.svg?branch=master
   :target: https://travis-ci.org/mcflugen/plume

.. image:: https://ci.appveyor.com/api/projects/status/yle29j1hl6a8yu8p?svg=true
   :target: https://ci.appveyor.com/project/mcflugen/plume

.. image:: https://coveralls.io/repos/github/mcflugen/plume/badge.svg?branch=mcflugen%2Fadd-unit-tests
   :target: https://coveralls.io/github/mcflugen/plume?branch=master

==================================================
plume: A hypopycnal plume model built with landlab
==================================================


Requirements
------------

*plume* requires Python 3.

Apart from Python, *plume* has a number of other requirements, all of which
can be obtained through either *pip* or *conda*, that will be automatically
installed when you install *plume*.

To see a full listing of the requirements, have a look at the project's
*requirements.txt* file.

If you are a developer of *plume* you will also want to install
additional dependencies for running *plume*'s tests to make sure
that things are working as they should. These dependencies are listed
in *requirements-testing.txt*.

Installation
------------

To install *plume*, first create a new environment in
which *plume* will be installed. This, although not necessary, will
isolate the installation so that there won't be conflicts with your
base *Python* installation. This can be done with *conda* as:

.. code:: bash

    $ conda create -n plume python=3
    $ conda activate plume

Stable Release
--------------

*plume*, and its dependencies, can be installed either with *pip*
or *conda*. Using *pip*:

.. code:: bash

    $ pip install plume

Using *conda*:

.. code:: bash

    $ conda install plume -c conda-forge

From Source
```````````

Before building *plume* from source, you will need an installation of
the GNU Scientific Library (gsl). There are several ways to install
this but the easiest is through *conda*,

.. code:: bash

   $ mamba install gsl

After downloading the *plume* source code, run the following from
*plume*'s top-level folder (the one that contains *setup.py*) to
install *plume* into the current environment:

.. code:: bash

    $ pip install -e .

Input Files
-----------

Configuration File
``````````````````

The main *plume* input file is a yaml-formatted text file that lists
constants used by *plume*. Running the following will print a sample
*plume* configuration file:

.. code:: bash

    $ plume generate plume.toml

This will print something like the following,

.. code:: toml

   [plume]
   _version = "0.2.0.dev0"

   [plume.grid]
   shape = [500, 500]
   xy_spacing = [100.0, 100.0]
   xy_of_lower_left = [0.0, 0.0]

   [plume.river]
   filepath = "river.csv"
   width = 50.0
   depth = 5.0
   velocity = 1.5
   location = [0.0, 25000.0]
   angle = 0.0

   [plume.sediment]
   removal_rate = 60.0
   bulk_density = 1600.0

   [plume.ocean]
   filepath = "ocean.csv"
   along_shore_velocity = 0.1
   sediment_concentration = 0.0

   [plume.output]
   filepath = "plume.nc"

Ocean File
``````````

The *plume* ocean file defines parameters of the ocean for each day of
the simulation. This is a csv-formatted text file to *day*, *along-shore velocity*,
and *sediment concentration*.

.. code:: bash

    $ plume generate ocean.csv

.. code:: csv

   # version: 0.2.0.dev0
   # Time [d], Along-shore velocity [m/s], Sediment Concentration [-]
   0.0,0.1,0.0

River File
``````````

The *plume* river file is a csv-formatted text file that gives river parameters
for each day of the simulation. Columns are *time*, *river width*, *river depth*,
and *river velocity*.

.. code:: bash

  $ plume generate river.csv

.. code:: csv

  # version: 0.2.0.dev0
  # Time [d], Width [m], Depth [m], Velocity [m/s]
  0.0,50.0,5.0,1.5

The *plume* river file defines

Output File
-----------

The only output file of *plume* is a *netCDF* file that contains
sediment concentrations for each day of the simulation.

Examples
--------

To run a simulation using the sample input files described above, you first
need to create a set of sample files:

.. code:: bash

    $ mkdir example
    $ plume --cd=example setup

You can now run the simulation:

.. code:: bash

    $ plume --cd=example run
    # Layer Thickness [m], Porosity [-]
    100.0,0.5
    96.18666488709239,0.4801774231522433
    92.78860257194452,0.4611407154102571
