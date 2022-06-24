import os
import pathlib
import sys
from distutils.core import setup

import cython_gsl
import pkg_resources
from Cython.Distutils import Extension, build_ext

os.environ.setdefault("LIB_GSL", str(pathlib.Path(sys.prefix) / "Library"))


setup(
    include_dirs=[
        cython_gsl.get_include(),
        pkg_resources.resource_filename("numpy", "core/include"),
    ],
    cmdclass={"build_ext": build_ext},
    ext_modules=[
        Extension(
            "plume.ext.centerline",
            ["plume/ext/centerline.pyx"],
            extra_compile_args=["-O3"],
            libraries=cython_gsl.get_libraries(),
            library_dirs=[cython_gsl.get_library_dir()],
            include_dirs=[
                cython_gsl.get_cython_include_dir(),
                cython_gsl.get_include(),
            ],
        )
    ],
)
