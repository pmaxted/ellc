#!/usr/bin/env python
import os
import shlex
import subprocess
import tempfile

import setuptools
from setuptools.command.build_py import build_py
from setuptools.dist import Distribution

# Same files and order as the makefile (Fortran modules must be compiled
# before the files that use them).
F_SRC = ['constants.f90', 'utils.f90', 'coords.f90', 'gauss_legendre.f90',
         'solve_real_poly.f90', 'ellipse.f90', 'stellar.f90', 'spots.f90',
         'ellc.f90']


class BuildPyWithLibellc(build_py):
    """Compile f_src/ into ellc/libellc.so for the machine doing the build.

    The compiler is $FC (default gfortran). $FFLAGS replaces the default
    optimisation flags and $LDFLAGS is appended to the link line.
    """

    def run(self):
        super().run()
        here = os.path.dirname(os.path.abspath(__file__))
        target = os.path.join(self.build_lib, 'ellc', 'libellc.so')
        os.makedirs(os.path.dirname(target), exist_ok=True)
        fc = os.environ.get('FC', 'gfortran')
        fflags = shlex.split(os.environ.get('FFLAGS', '-O2'))
        ldflags = shlex.split(os.environ.get('LDFLAGS', ''))
        sources = [os.path.join(here, 'f_src', f) for f in F_SRC]
        with tempfile.TemporaryDirectory() as moddir:
            cmd = ([fc] + fflags + ['-fPIC', '-shared', '-J', moddir]
                   + sources + ['-o', target] + ldflags)
            print(' '.join(cmd))
            subprocess.check_call(cmd)


class BinaryDistribution(Distribution):
    """libellc.so is platform specific, so the wheel must not be tagged
    py3-none-any."""

    def has_ext_modules(self):
        return True


if __name__ == '__main__':
    setuptools.setup(
        name='ellc',
        version='1.8.15',
        author='Pierre Maxted',
        author_email='p.maxted@keele.ac.uk',
        license='GNU GPLv3',
        scripts=['bin/ellc_emcee'],
        url='https://github.com/pmaxted/ellc',
        description='Light curve model for eclipsing binary stars and transiting exoplanets',
        classifiers = [
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering',
          'Programming Language :: Python',
          'Programming Language :: Fortran'],
        install_requires=["numpy >= 1.10.0","astropy >= 1.1.1", "scipy",
                          "emcee", "corner", "matplotlib"],
        packages=['ellc'],
        package_data={
            'ellc': ['data/*', 'doc/*', 'examples/*']
        },
        include_package_data=True,
        distclass=BinaryDistribution,
        cmdclass={
            'build_py': BuildPyWithLibellc,
        }
    )
