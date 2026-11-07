project(
  'ellc',
  'fortran',
)

py = import('python').find_installation()

py.install_sources(
  [
    'ellc/__init__.py',
    'ellc/lc.py',
    'ellc/rv.py',
    'ellc/ldy.py',
    'ellc/fluxes.py',
  ],
  subdir : 'ellc'
)

install_subdir(
  'ellc/data',
  install_dir : py.get_install_dir() / 'ellc'
)

libellc = shared_library(
  'ellc',
  files(
    'f_src/constants.f90',
    'f_src/utils.f90',
    'f_src/coords.f90',
    'f_src/gauss_legendre.f90',
    'f_src/solve_real_poly.f90',
    'f_src/ellipse.f90',
    'f_src/stellar.f90',
    'f_src/spots.f90',
    'f_src/ellc.f90',
  ),
  install: true,
  install_dir: py.get_install_dir() / 'ellc'
)
