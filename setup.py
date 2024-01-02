from __future__ import annotations

import os
import pathlib
import sys

import cython_gsl
import numpy
from Cython.Build import cythonize
from setuptools import Extension
from setuptools import setup

os.environ.setdefault("LIB_GSL", str(pathlib.Path(sys.prefix) / "Library"))


setup(
    include_dirs=[cython_gsl.get_include(), numpy.get_include()],
    ext_modules=cythonize(
        [
            Extension(
                "plume.ext.centerline",
                ["src/plume/ext/centerline.pyx"],
                extra_compile_args=["-O3"],
                libraries=cython_gsl.get_libraries(),
                library_dirs=[cython_gsl.get_library_dir()],
                include_dirs=[
                    cython_gsl.get_cython_include_dir(),
                    cython_gsl.get_include(),
                ],
            )
        ]
    ),
)
