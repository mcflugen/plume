=============
Release Notes
=============

.. towncrier release notes start

0.2.1 (2022-06-29)
------------------

Bug Fixes
`````````

- Fixed an error when trying to get the version from the command line with
  ``plume --version```. (`#9 <https://github.com/mcflugen/plume/issues/9>`_)


0.2.0 (2022-06-28)
------------------

New Features
````````````

- Updated the *Plume* component to be compatible with *landlab* version 2. This
  will allow *plume* to operate with the newest version of *landlab* and be
  incorporated into other *landlab* frameworks such as the *Sequence* model. (`#5 <https://github.com/mcflugen/plume/issues/5>`_)
- Change the *plume* command-line program to use *toml*-formatted input files. (`#5 <https://github.com/mcflugen/plume/issues/5>`_)
- The *plume* model now builds and runs on Windows! (`#5 <https://github.com/mcflugen/plume/issues/5>`_)


Bug Fixes
`````````

- Fixed a bug where the existing deposit thickness was being zeroed-out before
  the plume was run. The plume now adds sediment to an existing deposit and
  leaves it up to the user to clear an existing deposit. (`#7 <https://github.com/mcflugen/plume/issues/7>`_)


Documentation Enhancements
``````````````````````````

- Updated the documentation to include installation instructions, usage and
  examples. (`#5 <https://github.com/mcflugen/plume/issues/5>`_)


Other Changes and Additions
```````````````````````````

- Setup towncrier to manage the changelog. (`#5 <https://github.com/mcflugen/plume/issues/5>`_)
- Switch from *versioneer* to *zest.releaser* to manage release versions. (`#5 <https://github.com/mcflugen/plume/issues/5>`_)
- Set up continuous integration workflows using GitHub Actions. CI tests include
  testing for code style using *black*, checking for lint with *flake8*,
  testing notebooks and running the test suite. (`#5 <https://github.com/mcflugen/plume/issues/5>`_)
- Updated the package metadata and moved static metaddata from *setup.py*
  to *pyproject.toml*. (`#5 <https://github.com/mcflugen/plume/issues/5>`_)
- Setup up pre-commit hooks for ensuring style, lack of lint, python 3.8+
  syntax, and clean notebooks. (`#5 <https://github.com/mcflugen/plume/issues/5>`_)
- Added a citation file that users of *plume* can use to cite the software. (`#6 <https://github.com/mcflugen/plume/issues/6>`_)
- Removed ocean and river time series from the plume command line interface.
  The plume command now runs just a single plume. (`#7 <https://github.com/mcflugen/plume/issues/7>`_)
- Added GitHub Actions workflows for pushing prereleases to TestPyPI and
  releases to PyPI. (`#8 <https://github.com/mcflugen/plume/issues/8>`_)


0.1.0 (2018-06-13)
------------------

Features
````````

- Create *plume* package that simulates 1D (i.e. centerline), and 2D hypopycnal
  plumes.

Documentation Enhancements
``````````````````````````

- Added a new *Jupyter* notebook that demonstrates how to use simulate 2D and 1D
  plumes using the *plume* package. (`#3 <https://github.com/mcflugen/plume/issues/3>`_)


Other Changes and Additions
```````````````````````````

- Set up continuous integration testing with Travis-CI. CI builds are run
  for Linux and Mac. Windows support is not yet available. (`#1 <https://github.com/mcflugen/plume/issues/1>`_)
- Added a suite of tests, which use *pytest*, for testing the *plume*
  package. (`#2 <https://github.com/mcflugen/plume/issues/2>`_)
- Added *versioneer* to manage package release versions. (`#4 <https://github.com/mcflugen/plume/issues/4>`_)
