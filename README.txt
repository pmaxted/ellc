# ellc

`ellc` calculates light curves and radial velocity curves for eclipsing binary stars and transiting exoplanet systems.

## Installation

Install the latest release from PyPI:

```bash
pip install ellc
```

### Windows

Building `ellc` from source on Windows requires a Fortran compiler.

A convenient option is to install the MinGW-w64 toolchain from WinLibs:

https://winlibs.com/

After installation, ensure that the directory containing `gfortran.exe` is on your `PATH`, e.g.

```text
C:\mingw64\bin
```

Check that the compiler is available:

```bash
gfortran --version
```

You should then be able to install `ellc` using:

```bash
pip install ellc
```

or from a local source tree:

```bash
pip install .
```

### Linux

On Debian/Ubuntu:

```bash
sudo apt install gfortran
```

On Fedora:

```bash
sudo dnf install gcc-gfortran
```

### macOS

Install a Fortran compiler using Homebrew:

```bash
brew install gcc
```

## Requirements

- Python 3.10 or later
- NumPy
- SciPy
- Astropy

Required Python dependencies are installed automatically by `pip`.

## Examples

Example scripts are available in the `examples/` directory.

## Citation

If you use `ellc` in published work, please cite:

> Maxted, P.F.L. (2016), A&A, 591, A111

ADS bibcode:

```text
2016A&A...591A.111M
```

## License

See `LICENSE.txt`.
