#!/usr/bin/env python
import setuptools
from os.path import join

def configuration(parent_package='', top_path=None):
  from numpy.distutils.misc_util import Configuration

  config = Configuration('ellc', parent_package, top_path)

  # Order of source files in order of dependencies.
  sources = [ 'ellc_f.pyf',
     join('f_src','constants.f90'),
     join('f_src','utils.f90'),
     join('f_src','coords.f90'),
     join('f_src','gauss_legendre.f90'),
     join('f_src','solve_real_poly.f90'),
     join('f_src','ellipse.f90'),
     join('f_src','stellar.f90'),
     join('f_src','spots.f90'),   
     join('f_src','ellc.f90')
  ]

  config.add_extension('ellc_f', sources=sources)

  config.add_data_dir('data')

  config.add_data_dir('doc')

  config.add_data_dir('examples')

  return config

if __name__ == '__main__':
    from numpy.distutils.core import setup

    try:
      from version import __version__
    except:
      __version__ = ''

    setup(
        version=__version__,
        author='Pierre Maxted',
        author_email='p.maxted@keele.ac.uk',
        license='GNU GPLv3',
        scripts=['bin/ellc_emcee'],
        url='http://sourceforge.net/projects/goodricke/',
        description='Light curve model for eclipsing binary stars and transiting exoplanets',
        classifiers = [
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering',
          'Programming Language :: Python',
          'Programming Language :: Fortran'],
          install_requires=["numpy >= 1.10.0","astropy >= 1.1.1"], 
        **configuration(top_path='').todict())

