#!/usr/bin/env python
import setuptools
import subprocess
# from os.path import join


def run_make():
  subprocess.check_call(['make'])
  subprocess.check_call(['make', 'install'])


def configuration(parent_package='', top_path=None):
  run_make()

  package_data = {
    'ellc': ['data/*', 'doc/*', 'examples/*']
  }

  return [], package_data


if __name__ == '__main__':
    try:
      from version import __version__
    except:
      __version__ = ''

    ext_modules, package_data = configuration()

    setuptools.setup(
        name='ellc',
        version=__version__,
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
        ext_modules=ext_modules,
        packages=['ellc'],
        package_data=package_data,
        include_package_data=True
    )
