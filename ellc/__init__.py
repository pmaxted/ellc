"""
  ellc
  ----

  Generate light curves or radial velocity curves for a binary star system 
  with the ellc light curve model [1].

  Routines in this module:

  lc(t_obs, radius_1, radius_2, sbratio, incl, ... ) 

  rv(t_obs, radius_1, radius_2, sbratio, incl, t_zero, period, a, q, ...)

  fluxes(t_obs, radius_1, radius_2, sbratio, incl, ... ) 

  tmin(t_obs, radius_1, radius_2, sbratio, incl, ... ) 

  Sub-packages in this module

  ldy
   - ldy.list_bands()
   - ldy.LimbGravityDarkeningCoeffs(band)


  Please cite Maxted (2016) if you published results based on this software.

  References
  ----------
  .. [1] Maxted, P.F.L. 2016. A fast, flexible light curve model for detached
      eclipsing binary stars and transiting exoplanets. A&A VOLUME, ARTICLE.

"""

from .version import __version__
from .lc import lc
from .rv import rv
# from .tmin import tmin
from .fluxes import fluxes
from . import ldy