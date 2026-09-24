#!/usr/bin/env python
import setuptools
import subprocess
from setuptools.command.build_ext import build_ext

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

class CustomBuildExt(build_ext):
    def run(self):
        subprocess.check_call(['make'])
        subprocess.check_call(['make', 'install'])
        super().run()

if __name__ == '__main__':
    setuptools.setup(
        name='ellc',
        version='1.8.11',
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
        cmdclass={
            'build_ext': CustomBuildExt,
        }
    )
