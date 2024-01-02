from __future__ import annotations

import os
import pathlib
import sys

import cython_gsl
import numpy
from setuptools import Extension
from setuptools import setup

try:
    from Cython.Build import cythonize
except ModuleNotFoundError:
    WITH_CYTHON = False
else:
    WITH_CYTHON = True

os.environ.setdefault("LIB_GSL", str(pathlib.Path(sys.prefix) / "Library"))

ext = ".pyx" if WITH_CYTHON else ".c"

extension = Extension(
    "plume.ext.centerline",
    ["src/plume/ext/centerline" + ext],
    extra_compile_args=["-O3"],
    libraries=cython_gsl.get_libraries(),
    library_dirs=[cython_gsl.get_library_dir()],
    include_dirs=[
        cython_gsl.get_cython_include_dir(),
        cython_gsl.get_include(),
    ],
)

if WITH_CYTHON:
    extensions = cythonize(extension)

setup(
    include_dirs=[cython_gsl.get_include(), numpy.get_include()],
    ext_modules=extensions,
)
