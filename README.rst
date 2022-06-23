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

    $ compact show config
    c: 5.0e-08
    porosity_max: 0.5
    porosity_min: 0.0
    rho_grain: 2650.0
    rho_void: 1000.0

Porosity File
`````````````

The *plume* porosity file defines initial porosity of each of the
sediment layers to be compacted as a two-column CSV file. The first
column is layer thickness (in meters) and the second the porosity of
the sediment in that layer. A sample porosity file can be obtained with:

.. code:: bash

    $ compact show porosity
    # Layer Thickness [m], Porosity [-]
    100.0,0.5
    100.0,0.5
    100.0,0.5

Output File
-----------

The output file of *plume* is a porosity file of the same form as
the input porosity file - a CSV file of layer thickness and porosity.

Examples
--------

To run a simulation using the sample input files described above, you first
need to create a set of sample files:

.. code:: bash

    $ plume setup example
    example/plume.yaml

You can now run the simulation:

.. code:: bash

    $ plume run example/config.yaml
    # Layer Thickness [m], Porosity [-]
    100.0,0.5
    96.18666488709239,0.4801774231522433
    92.78860257194452,0.4611407154102571
