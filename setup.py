from setuptools import setup, Extension

import cython_gsl


setup(
    ext_modules=[
        Extension(
            "plume.ext.centerline",
            ["plume/ext/centerline.pyx"],
            extra_compile_args=["-O3"],
            libraries=cython_gsl.get_libraries(),
            library_dirs=[cython_gsl.get_library_dir()],
            include_dirs=[cython_gsl.get_cython_include_dir()],
        )
    ],
)
