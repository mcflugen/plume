from setuptools import setup, find_packages, Extension

import numpy as np
import cython_gsl
import versioneer


setup(
    name="plume",
    version=versioneer.get_version(),
    description="A hypopycnal sediment-carrying plume entering the ocean",
    long_description=open("README.md", encoding="utf-8").read(),
    author="Eric Hutton",
    author_email="huttone@colorado.edu",
    url="http://csdms.colorado.edu",
    install_requires=open("requirements.txt", "r").read().splitlines(),
    setup_requires=[
        "setuptools",
    ],
    packages=find_packages(),
    include_dirs=[np.get_include(), cython_gsl.get_include()],
    entry_points={
        "console_scripts": [
            "plume=plume.cli:plume",
        ],
    },
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
    cmdclass=versioneer.get_cmdclass(),
)
