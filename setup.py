from setuptools import setup, find_packages
from distutils.extension import Extension

import numpy as np
import cython_gsl
import versioneer


def read_requirements():
    import os

    path = os.path.dirname(os.path.abspath(__file__))
    requirements_file = os.path.join(path, 'requirements.txt')
    try:
        with open(requirements_file, 'r') as req_fp:
            requires = req_fp.read().split()
    except IOError:
        return []
    else:
        return [require.split() for require in requires]


setup(name='plume',
      version=versioneer.get_version(),
      description='A hypopycnal sediment-carrying plume entering the ocean',
      author='Eric Hutton',
      author_email='huttone@colorado.edu',
      url='http://csdms.colorado.edu',
      install_requires=read_requirements(),
      setup_requires=['setuptools', ],
      packages=find_packages(),
      include_dirs = [np.get_include(), cython_gsl.get_include()],
      entry_points={
          'console_scripts': [
              'plume=plume.cli:plume',
          ],
      },
      ext_modules = [
          Extension('plume.ext.centerline',
                    ['plume/ext/centerline.pyx'],
                    extra_compile_args=['-O3'],
                    libraries=cython_gsl.get_libraries(),
                    library_dirs=[cython_gsl.get_library_dir()],
                    include_dirs=[cython_gsl.get_cython_include_dir()])],
      cmdclass=versioneer.get_cmdclass(),
)
