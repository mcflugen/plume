from distutils.core import setup
from Cython.Distutils import Extension
from Cython.Distutils import build_ext
import cython_gsl
import pkg_resources


setup(
    include_dirs=[
        cython_gsl.get_include_dir(),
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
            include_dirs=[cython_gsl.get_cython_include_dir()],
        )
    ],
)
