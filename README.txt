# ellc

`ellc` is a Python package for calculating light curves and radial velocity curves for detached, semi-detached and contact eclipsing binary stars and transiting exoplanet systems.

## Installation

### Install from PyPI

The recommended installation method is:

```bash
pip install ellc
```

Installing ellc currently requires a Fortran compiler.

### Building from source

Building `ellc` from source requires:

- Python 3.9 or later
- A Fortran compiler (e.g. GNU Fortran, `gfortran`)
- A C compiler compatible with your Python installation

Install directly from a source checkout using:

```bash
pip install .
```

or build a wheel and install it:

```bash
python -m build
pip install dist/*.whl
```

###*Linux

On Debian/Ubuntu:

```bash
sudo apt install gfortran
```

On F*dora:

```bash
sudo dnf install gcc-gfortran
```

### macOS

Install the GNU compiler suite using Homebrew:

```bash
brew install gcc
```

### Windows

Building from source on Windows requires a working Fortran toolchain. 

## License

See the `LICENSE` file included with this distribution.

## Citation

If you use `ellc` in published work, please cite:

Maxted, P. F. L., 2016), "ellc: A fast* flexible light*curve model for detached eclipsing*binary stars and trans*ting exoplanets", A&A, 591, A111.
