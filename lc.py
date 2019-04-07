# This file is part of the ellc binary star model
# Copyright (C) 2017 Pierre Maxted
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np

from ellc import ellc_f

def lc(t_obs, radius_1, radius_2, sbratio, incl, 
       light_3 = 0, 
       t_zero = 0, period = 1,
       a = None,
       q = 1,
       f_c = None, f_s = None,
       ldc_1 = None, ldc_2 = None,
       gdc_1 = None, gdc_2 = None,
       didt = None,
       domdt = None,
       rotfac_1 = 1, rotfac_2 = 1,
       hf_1 = 1.5, hf_2 = 1.5,
       bfac_1 = None, bfac_2 = None,
       heat_1 = None, heat_2 = None, 
       lambda_1 = None, lambda_2 = None,
       vsini_1 = None, vsini_2 = None,
       t_exp=None, n_int=None, 
       grid_1='default', grid_2='default',
       ld_1=None, ld_2=None,
       shape_1='sphere', shape_2='sphere',
       spots_1=None, spots_2=None, 
       exact_grav=False, verbose=1):
  """
  Calculate the light curve of a binary star

  This function calculates the light curve of a binary star using the ellc
  binary star model [1].

  Parameters
  ----------
  t_obs : array_like
      Times or phases of observation. The units and time system used must be
      consistent with t_zero and period.

  radius_1 : float
      Radius of star 1 in units of the semi-major axis of the binary.
      The radius is defined to be the same as a sphere with the same volume as
      the ellipsoid used to approximate the shape of the star.
      Set radius_1=1 to fix radius at limiting radius in the Roche potential.

  radius_2 : float
      Radius of star 2 in units of the semi-major axis of the binary.
      The radius is defined to be the same as a sphere with the same volume as
      the ellipsoid used to approximate the shape of the star.
      Set radius_2=1 to fix radius at limiting radius in the Roche potential.

  sbratio : float
      Surface brightness ratio, S_2/S_1 

  incl : float
      Inclination in degrees.

  light_3 : float, optional
      Third light contribution relative to total flux from both stars at 
      time t_zero excluding eclipse effects.

  t_zero : float, optional
      Time (or phase) of mid-eclipse for star 1 by star 2.
      The units and time system must be consistent with the values in t_obs.
      Default is 0.

  period : float, optional
      Orbital period of the binary or (for phased data) 1.
      For binary systes with apsidal motion, this is the anomalistic period.
      If light-time effect or Doppler boosting is to be calculated correctly
      then the units of period must be days, otherwise arbitrary units can be
      used provided that these are consistent with t_zero and t_obs.
      Default is 1.

  a : {None, float}, optional
      Semi-major axis in solar radii for calculation of light travel time
      and Doppler boosting.

  q : float, optional
      Mass ratio m_2/m_1.
      Default is 1.

  f_c : {None, float},  optional
      For eccentric orbits with eccentricity e and longitude of periastron w,
      f_c = sqrt(e).cos(w)
      If not None, then f_s must also be specified.
      Default is None.

  f_s : {None, float},  optional
      For eccentric orbits with eccentricity e and longitude of periastron w,
      f_s = sqrt(e).sin(w)
      If not None, then f_c must also be specified.
      Default is None.

  ldc_1 : {None, float, or array_like}, optional
      Limb darkening coefficients for star 1.
      Number of elements must match the selected limb-darkening law.
      For the option ld_1='mugrid', ldc_1 is a grid of specific
      intensity values on a uniform grid from mu=0 (first element) to mu=1
      (last element). The specific intensies are assumed to be normalised,
      i.e., ldc_1[-1] = 1, but this is not checked.
      Default is None

  ldc_2 : {None, float, or array_like}, optional
      Limb darkening coefficients for star 2.
      Number of elements must match the selected limb-darkening law.
      For the option ld_2='mugrid', ldc_2 is a grid of specific
      intensity values on a uniform grid from mu=0 (first element) to mu=1
      (last element). The specific intensies are assumed to be normalised,
      i.e., ldc_2[-1] = 1, but this is not checked.
      Default is None

  gdc_1 : {None, float},  optional
      Gravity darkening exponent for star 1.

  gdc_2 : {None, float},  optional
      Gravity darkening exponent for star 2.

  didt : {None, float},  optional
      Rate of change of inclination [degrees/anomalistic period]
      Default is None.

  domdt : {None, float},  optional
      Apsidal motion rate  [degrees/anomalistic period].
      Default is None.

  rotfac_1 : float, optional
      Asynchronous rotation factor for star 1, F_1
      Default is 1.

  rotfac_2 : float, optional
      Asynchronous rotation factor for star 2, F_2
      Default is 1.

  hf_1 : float, optional
      h_f, fluid second Love number for radial displacement, for star 1.
      Only used if shape_1 = 'love'
      Default is 1.5.

  hf_2 : float, optional
      h_f, fluid second Love number for radial displacement, for star 2.
      Only used if shape_2 = 'love'
      Default is 1.5.


      N.B. Doppler boosting is not calculated is parameter a is None
      Default is None.
      
  bfac_2 : {None, float}, optional
      Doppler boosting factor, star 2
      N.B. Doppler boosting is not calculated is parameter a is None
      Default is None.
      
  heat_1 : {None, scalar, 3-tuple of floats}, optional
      If scalar, coefficient of simplified reflection model - see Notes.
      If 3-tuple, parameters of heating+reflection model, [H_0, H_1, u_H]
      H_0 is the coefficient, H_1 is the exponent and u_H is the linear
      limb-darkening coefficient.
      Default is None.
  
  heat_2 : {None, scalar, 3-tuple of floats}, optional
      If scalar, coefficient of simplified reflection model - see Notes.
      If 3-tuple, parameters of heating+reflection model, [H_0, H_1, u_H]
      H_0 is the coefficient, H_1 is the exponent and u_H is the linear
      limb-darkening coefficient.
      Default is None.
  
  lambda_1 : {None, float},  optional
       Sky-projected angle between orbital and rotation axes, star 1 [degrees]
       N.B. lambda_1 is only used if shape_1='sphere'
       Default is None.
  
  lambda_2 : {None, float},  optional
       Sky-projected angle between orbital and rotation axes, star 2 [degrees]
       N.B. lambda_2 is only used if shape_2='sphere'
       Default is None.

  vsini_1  : {None, float}, optional
      V_rot.sini for calculation of R-M effect for star 1 [km/s]
      See notes below.
      Default is None.

  vsini_2  : {None, float}, optional
      V_rot.sini for calculation of R-M effect for star 2 [km/s]
      See notes below.
      Default is None.

  t_exp : {None, float, array_like}, optional
      Exposure time in the same units as t_obs and, if array_like, with the
      same number of elements.
      Default is None.

  n_int : {None, int, array_like}
      Number of integration points used to account for finite exposure time.
      Set n_int or elements of n_int to 1 to make no correction for finite
      integration time.
      If array_like then set elements to 0 to use linear interpolation of the
      other values in the light curve to estimate the light curve at the
      corresponding values of t_obs. 
      Default is None.

  grid_1 : {"very_sparse", "sparse", "default", "fine", "very_fine"}, optional
      Grid size used to calculate the flux from star 1.
      Default is "default"

  grid_2 : {"very_sparse", "sparse", "default", "fine", "very_fine"}, optional
      Grid size used to calculate the flux from star 2.
      Default is "default"

  ld_1 : {None, "lin", "quad", "sing", "claret", "log", "sqrt", "exp",
          "power-2", "mugrid"} 
   Limb darkening law for star 1
   The "power-2" limb darkening law is taken from  Morello et al., 
   2017AJ....154..111M:
    I_lambda(mu)/I_lambda(1) = 1 - ldc_1[0] * (1 - mu**ldc_1[1])
   For the option "mugrid", ldc_1 must be an array of specific intensity
   values as a function of mu=cos(theta), where theta is the angle between the
   line of sight and the normal to a given point of the stellar surface. The
   value of mu corresponding to element i of ldc_1 is i/(len(ldc_1)-1),  where
   i=0,1, ... len(ldc_1)-1.  
   Default is None

  ld_2 : {None, "lin", "quad", "sing", "claret", "log", "sqrt", "exp",
          "power-2", "mugrid"} 
   Limb darkening law for star 2.
   See ld_1 for description of the power-2 and mugrid options.
   Default is None

  shape_1 : {"roche", "roche_v", "sphere", "poly1p5", "poly3p0", "love"}
      Model used to calculate the shape of star 1 - see Notes. 
      Default is "sphere".

  shape_2 : {"roche", "roche_v", "sphere", "poly1p5", "poly3p0", "love"}
      Model used to calculate the shape of star 2 - see Notes. 
      Default is "sphere".

  spots_1 : (4, n_spots_1) array_like
   Parameters of the spots on star 1. For each spot the parameters, in order,
   are longitude, latitude, size and brightness factor. All three angles are
   in degrees.

  spots_2 : (4, n_spots_2) array_like
   Parameters of the spots on star 2. For each spot the parameters, in order,
   are longitude, latitude, size and brightness factor. All three angles are
   in degrees.

  exact_grav : {True|False}
    Use point-by-point calculation of local surface gravity for calculation of
    gravity darkening is True, otherwise use (much faster) approximation based
    on functional form fit to local gravity at 4 points on the star.

  Returns
  -------
  flux : ndarray
      Flux in arbitrary units.

  Notes
  -----

   The asynchronous rotation factors rotfac_1 and rotfac_2 are used to
  calculate the shapes of the stars (unless spherical stars are specified).
  These rotation factors are relative to the actual synchronous rotation rate
  (rotation period = orbital period) not pseudo-synchronous rotation rate.

   The effect of the spot on the light curve is calculated using the
  algorithm by Eker [2] for circular spots on a spherical star with
  quadratic limb darkening. If the limb-darkening law used for the main
  calculation is not linear or quadratic then the coefficients of the
  limb-darkening law used for the calculation of the effects of the spots
  are set so that the intensity distribution matches at mu = 0, 0.5 and 1.

   N.B. The effect of each spot on the light curve is additive so overlapping
  spots can result in non-physical negative fluxes for some regions of the
  star.
 
   For the calculation of the star shape, the rotation and orbital angular
  momentum vectors are assumed parallel. The shape of each star is
  approximated by a triaxial ellipsoid with semi-major axes (A,B,C) and
  offset towards companion, D, from the centre-of-mass of the star towards
  the companion star. For the option "roche" the definition of the Roche
  potential from Wilson [3] is used and the values of A, B, C, D are set to
  that intersection points of the triaxial ellipsoid with x-, y- and z-axes
  lie on an equipotential surface. For the options "poly1p5" and "poly3p0"
  the star is assumed to behave as a polytrope with index n=1.5 or n=3,
  respectively. The tidal and rotational distortion of the polytrope are 
  assumed to be independent. The tidal distortion for polytropes is from
  Chandrasekhar [4] and the rotational distortion is calculated by
  interpolation in Table 1 of James [5]. For the option "love" the shape of
  the star is calculated using equation (10) and (11) from Correia [6] using
  the fluid second Love number for radial displacement, h_f (Correia equation
  (8)). The offset, D, is calculated using the approximation D=q.(R/d)^4,
  i.e., equation (38) from Chandrasekhar [4] with Delta_3=1 and nu = radius/d, 
  where d is the separation of the stars and q is the mass ratio.

    In eccentric orbits, the volume of the star is assumed to be constant. In
  general, the volume is calculated from the volume of the approximating
  ellipsoid. In the case of synchronous rotation, the volume of the star can
  be calculated using equation (2.18) from Kopal "Dynamics of Close Binary
  Systems" (Springer, 1978) by selecting the star shape model "roche_v".
   
   The simplified reflection model is approximately equivalent to Lambert
  law scattering with the coefficients heat_1 and heat_2  being equal to
  A_g/2, where A_g is the geometric albedo.

  Example
  -------
  >>> import ellc
  >>> import numpy as np
  >>> import matplotlib.pyplot as plt
  >>> t = np.arange(-0.25,0.75, 0.001)
  >>> spots_1 = [[30,180],[45,-45],[25,35],[0.2,0.8]]
  >>> flux = ellc.lc(t,radius_1=0.1,radius_2=0.05,sbratio=0.2,
  ...   incl=89.95,q=0.5,ld_1='quad',ldc_1=[0.65,0.2],ld_2='lin',ldc_2=0.45,
  ...   shape_1='poly3p0',shape_2='poly1p5',spots_1=spots_1)
  >>> plt.plot(t,flux)
  >>> plt.show()

  References
  ----------
  .. [1] Maxted, P.F.L. 2016. A fast, flexible light curve model for detached
      eclipsing binary stars and transiting exoplanets. A&A 591, A111, 2016.
  .. [2] Eker, 1994, ApJ, 420, 373.
  .. [3] Wilson, 1979, ApJ, 234, 1054.
  .. [4] Chandrasekhar, 1933, MNRAS, 93, 449.
  .. [5] James, 1964, ApJ, 140, 552.
  .. [6] Correia, 2014, A&A 570, L5.

  """


  # Copy control parameters into an np.array

  gridname_to_gridsize = {
    "very_sparse" :  4,
    "sparse"       :  8,
    "default"      : 16,
    "fine"         : 24,
    "very_fine"   : 32,
  }
  n1 = gridname_to_gridsize.get(grid_1,None)
  if n1 is None:
    raise Exception("Invalid grid size name")
  n2 = gridname_to_gridsize.get(grid_2,None)
  if n2 is None:
    raise Exception("Invalid grid size name")

  ldstr_to_ldcode = {
    "none"   :  0,
    "lin"    :  1,
    "quad"   :  2,
    "sing"   :  3,
    "claret" :  4,
    "log"    :  5,
    "sqrt"   :  6,
    "exp"    :  7,
    "power-2":  8,
    "mugrid" : -1
  }
  if ld_1 is None:
    ldstr_1 = 'none'
  else:
    ldstr_1 = ld_1

  if ld_2 is None:
    ldstr_2 = 'none'
  else:
    ldstr_2 = ld_2

  l1 = ldstr_to_ldcode.get(ldstr_1,None)
  if l1 is None:
    raise Exception("Invalid limb darkening law name")

  l2 = ldstr_to_ldcode.get(ldstr_2,None)
  if l2 is None:
    raise Exception("Invalid limb darkening law name")

  shapename_to_shapecode = {
    "roche_v" : -2,
    "roche"   : -1,
    "sphere"  :  0,
    "poly1p5" :  1,
    "poly3p0" :  2,
    "love"    :  3
  }
  s1 = shapename_to_shapecode.get(shape_1,None)
  if s1 is None:
    raise Exception("Invalid star shape name")
  s2 = shapename_to_shapecode.get(shape_2,None)
  if s2 is None:
    raise Exception("Invalid star shape name")

  if spots_1 is None:
    spar_1 = np.zeros([1,1])
    n_spots_1 = 0
  else:
    spar_1 = np.array(spots_1)
    if (spar_1.ndim != 2) or (spar_1.shape[0] != 4 ):
      raise Exception("spots_1 is not  (4, n_spots_1) array_like")
    n_spots_1 = spar_1.shape[1]

  if spots_2 is None:
    spar_2 = np.zeros([1,1])
    n_spots_2 = 0
  else:
    spar_2 = np.array(spots_2)
    if (spar_2.ndim != 2) or (spar_2.shape[0] != 4 ):
      raise Exception("spots_2 is not  (4, n_spots_2) array_like")
    n_spots_2 = spar_2.shape[1]

  ipar = np.array([n1,n2,n_spots_1,n_spots_2,l1,l2,s1,s2,1,0+exact_grav],
      dtype=int)

  if ld_1 == 'mugrid':
      try:
          mugrid_1 = np.array(ldc_1)
          n_mugrid_1 = len(mugrid_1)
      except:
          raise Exception("failed to obtain limb darkening data from ldc_1")
      if (n_mugrid_1 < 2):
          raise Exception("ldc_1 must have at least 2 elements")
  else:
      mugrid_1 = np.array([0.])
      n_mugrid_1 = 0

  if ld_2 == 'mugrid':
      try:
          mugrid_2 = np.array(ldc_2)
          n_mugrid_2 = len(mugrid_2)
      except:
          raise Exception("failed to obtain limb darkening data from ldc_2")
      if (n_mugrid_2 < 2):
          raise Exception("ldc_2 must have at least 2 elements")
  else:
      mugrid_2 = np.array([0.])
      n_mugrid_2 = 0

  # Copy binary parameters into an np.array

  if (radius_1 <= 0) or (radius_1 > 1):
    raise ValueError("radius_1 argument out of range")
  if (radius_1 == 1) and (shape_1 != "roche"):
    raise ValueError("radius_1=1 only allowed for Roche potential")

  if (radius_2 <= 0) or (radius_2 > 1):
    raise ValueError("radius_2 argument out of range")
  if (radius_2 == 1) and (shape_2 != "roche"):
    raise ValueError("radius_2=1 only allowed for Roche potential")

  par = np.zeros(39)
  par[0] = t_zero
  par[1] = period
  par[2] = sbratio
  par[3] = radius_1
  par[4] = radius_2
  par[5] = incl
  par[6] = light_3

  if a is not None : par[7] = a

  if (f_c is None) and (f_s is None): 
    pass
  elif (f_c is not None) and (f_s is not None):
    par[8] = f_c
    par[9] = f_s
  else:
    raise Exception("Must specify both f_c and f_s or neither.")

  if q <= 0 :
    raise ValueError("Mass ratio q must be positive.")
  par[10] = q
  
  ld_to_n  = {
    "none"   :  0,
    "lin"    :  1,
    "quad"   :  2,
    "sing"   :  3,
    "claret" :  4,
    "log"    :  2,
    "sqrt"   :  2,
    "exp"    :  2,
    "power-2":  2,
    "mugrid" :  0 
  }
  ld_n_1 = ld_to_n.get(ldstr_1,None)
  try: 
      if (ld_n_1 > 0): 
          par[11:11+ld_n_1] = ldc_1
  except:
    raise Exception("ldc_1 and ld_1 are inconsistent")
  ld_n_2 = ld_to_n.get(ldstr_2,None)
  try: 
      if (ld_n_2 > 0): 
          par[15:15+ld_n_2] = ldc_2
  except:
    raise Exception("ldc_2 and ld_2 are inconsistent")

  if gdc_1 is not None : par[19] = gdc_1

  if gdc_2 is not None : par[20] = gdc_2

  if didt is not None : par[21] = didt

  if domdt is not None : par[22] = domdt

  par[23] = rotfac_1

  par[24] = rotfac_2

  if bfac_1 is not None : par[25] = bfac_1

  if bfac_2 is not None : par[26] = bfac_2

  if heat_1 is not None : 
    t = np.array(heat_1)
    if t.size == 1:
      par[27] = t
    elif t.size == 3:
      par[27:30] = t
    else:
      raise Exception('Invalid size for array heat_1')

  if heat_2 is not None : 
    t = np.array(heat_2)
    if t.size == 1:
      par[30] = t
    elif t.size == 3:
      par[30:33] = t
    else:
      raise Exception('Invalid size for array heat_2')

  if lambda_1 is not None : par[33] = lambda_1

  if lambda_2 is not None : par[34] = lambda_2

  if vsini_1 is not None : par[35] = vsini_1

  if vsini_2 is not None : par[36] = vsini_2

  if shape_1 == "love":
    if rotfac_1 != 1:
      raise Exception(
                'Non-synchronous rotation not implemented for shape_1="love"')
    if (hf_1 <= 2/3.) or (hf_1 > 5):
        raise Exception('Invalid value for hf_1')
    par[37] = hf_1

  if shape_2 == "love":
    if rotfac_2 != 1:
      raise Exception(
                'Non-synchronous rotation not implemented for shape_2="love"')
    if (hf_2 <= 2/3.) or (hf_2 > 5):
        raise Exception('Invalid value for hf_2')
    par[38] = hf_2
     
  t_obs_array = np.array(t_obs)
  n_obs = len(t_obs_array)
  if t_exp is None:
    t_exp_array = np.zeros(n_obs)
  else:
    t_exp_array = np.ones(n_obs)*t_exp

  if n_int is not None :
    if np.amax(n_int) < 1 : raise Exception("No n_int values > 1.") 
    if np.amin(n_int) < 0 : raise Exception("Invalid negative n_int value(s).") 
    n_int_array = np.array(np.ones(n_obs)*n_int, dtype=int)
  else:
    n_int_array = np.ones(n_obs, dtype=int)

  # Create list of times for calculation, weights for integration and
  # indices to relate these both back to the original t_obs array
  i_obs = np.arange(0,n_obs)
  t_calc = t_obs_array[n_int_array == 1]
  w_calc = np.ones_like(t_calc)
  i_calc = i_obs[n_int_array == 1]
  n_int_max = np.amax(n_int_array)
  for i_int in np.unique(n_int_array[n_int_array > 1]) :
    t_obs_i = t_obs_array[n_int_array == i_int]
    t_exp_i = t_exp_array[n_int_array == i_int]
    i_obs_i = i_obs[n_int_array == i_int]
    for j_int in range(0,i_int):
      t_calc = np.append(t_calc, t_obs_i+(j_int/(i_int-1.)-0.5)*t_exp_i)
      i_calc = np.append(i_calc,i_obs_i)
      if (j_int == 0) or (j_int == i_int-1):
        w_calc = np.append(w_calc, 0.5*np.ones_like(t_obs_i)/(i_int-1.))
      else:
        w_calc = np.append(w_calc, np.ones_like(t_obs_i)/(i_int-1.))

  lc_rv_flags = ellc_f.ellc.lc(t_calc,par,ipar,spar_1,spar_2,
                n_mugrid_1, mugrid_1,n_mugrid_2, mugrid_2, verbose)
  if ((np.sum(np.isnan(lc_rv_flags)) > 0 )) & (verbose > 0):
    lc_dummy = ellc_f.ellc.lc(t_calc,par,ipar,spar_1,spar_2,
               n_mugrid_1, mugrid_1,n_mugrid_2, mugrid_2,9)

  flux = np.zeros(n_obs)
  for j in range(0,len(t_calc)):
    flux[i_calc[j]] += lc_rv_flags[j,0]*w_calc[j]

  t_obs_0 = t_obs_array[n_int_array == 0 ] # Points to be interpolated
  n_obs_0 = len(t_obs_0)
  if n_obs_0 > 0 :
    i_sort = np.argsort(t_calc)
    t_int = t_calc[i_sort]
    f_int = lc_rv_flags[i_sort,0]
    flux[n_int_array == 0 ] = np.interp(t_obs_0,t_int,f_int)

  return flux
 

