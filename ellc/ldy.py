""" 
 Return limb-darkening and gravity-darkening coefficients

  The coefficients for a 4-parameter limb-darkening law and the gravity
 darkening coefficient, y, are interpolated from the tables provided by Claret
 and Bloemen (2011A&A...529A..75C) for ATLAS stellar atmosphere models and an
 assumed micro-turbulent velocity xi=2 km/s calculated using a least-squares
 fit. Also available are limb-darkening coefficients for the SDSS u', g', r',
 i' and z' bands from Claret (2004A&A...428.1001C), but no gravity darkening
 coefficient is available for these bands.


  To interpolate the limb-darkening and gravity-darkening coefficients at a
  given effective temperature (T_eff in K), surface gravity (log(g) in cgs
  units) and metallicity ([M/H] logarithmic metal abundance relative to
  solar), create an instance of the class LimbGravityDarkeningCoeffs for the
  desired photometric band and then call this instance with the parameters
  (T_eff, log(g), [M/H]). The first 4 values in the array returned are the
  limb-darkening coefficients a_1 .. a_4 and the final value is the
  gravity-darkening coefficient, y.  If no gravity-darkening coefficient is
  available then the fifth returned value will be None. 

  If the input values of T_eff, log(g) or [M/H] are outside the tabulated
  range then all coefficients are returned as NaN.

   To see the names of the available photometric bands use the function
   ellc.ldy.list_bands(). For the SDSS bands, 'u_' = u', 'g_' = g', etc.

 Example
 -------
  
 >>> from ellc.ldy import LimbGravityDarkeningCoeffs, list_bands
 >>> print(list_bands())
 ['B', 'C', 'H', 'I', 'J', 'K', 'Kp', 'R', 'S1', 'S2', 'S3', 'S4', 'U', 'V',
 'b', 'u', 'v', 'y', 'u_', 'g_', 'r_', 'i_', 'z_']
 >>> ldy_Kp = LimbGravityDarkeningCoeffs('Kp')
 >>> Teff = 10450
 >>> logg = 3.9
 >>> M_H = -0.31
 >>> a1,a2,a3,a4,y = ldy_Kp(Teff, logg, M_H)

"""
  
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from astropy.io import fits
from os.path import join,abspath,dirname
from scipy.interpolate import RegularGridInterpolator

dir_path = dirname(abspath(__file__))
ld_file_1 = join(dir_path,'data','g_Claret2011.fits')
ld_file_2 = join(dir_path,'data','g_Claret2004.fits')

def list_bands():
  hdulist = fits.open(ld_file_1)
  list_1 = (hdulist['band'].data)['BAND']
  hdulist.close()

  hdulist = fits.open(ld_file_2)
  list_2 = (hdulist['band'].data)['BAND']
  hdulist.close()

  return  [x.strip(' ') for x in  np.append(list_1, list_2)]

class LimbGravityDarkeningCoeffs:

  def __init__(self, band):
    self._band = band
    if band in ('u_','g_','r_','i_','z_'):
      hdulist = fits.open(ld_file_2)
      teff = hdulist['Teff'].data
    else:
      hdulist = fits.open(ld_file_1)
      logteff = hdulist['logTeff'].data
    logg = hdulist['logg'].data
    M_H = hdulist['M_H'].data
    bands = hdulist['band'].data
    i_band = (bands['BAND'].rfind(band) > -1).nonzero()[0][0]
    grid = (hdulist[0].data)[i_band,:,:,:,:]
    hdulist.close()
    if band in ('u_','g_','r_','i_','z_') :
      pts = (M_H,logg,teff)
      self._a1 = RegularGridInterpolator(pts,grid[0,:,:,:],fill_value=None)
      self._a2 = RegularGridInterpolator(pts,grid[1,:,:,:],fill_value=None)
      self._a3 = RegularGridInterpolator(pts,grid[2,:,:,:],fill_value=None)
      self._a4 = RegularGridInterpolator(pts,grid[3,:,:,:],fill_value=None)
    else:
      pts = (M_H,logg,logteff)
      self._y  = RegularGridInterpolator(pts,grid[0,:,:,:],fill_value=None)
      self._a1 = RegularGridInterpolator(pts,grid[1,:,:,:],fill_value=None)
      self._a2 = RegularGridInterpolator(pts,grid[2,:,:,:],fill_value=None)
      self._a3 = RegularGridInterpolator(pts,grid[3,:,:,:],fill_value=None)
      self._a4 = RegularGridInterpolator(pts,grid[4,:,:,:],fill_value=None)

  def __call__(self, teff, logg, M_H):
    if teff < 0:
      y,a1,a2,a3,a4  = [np.nan,np.nan,np.nan,np.nan,np.nan]
    else:
      if self._band in ('u_','g_','r_','i_','z_'):
        try: 
          y  = None
          a1 = self._a1([M_H,logg,teff])
          a2 = self._a2([M_H,logg,teff])
          a3 = self._a3([M_H,logg,teff])
          a4 = self._a4([M_H,logg,teff])
        except:
          y,a1,a2,a3,a4  = [np.nan,np.nan,np.nan,np.nan,np.nan]
      else:
        try:
          y  = self._y ([M_H,logg,np.log10(teff)])
          a1 = self._a1([M_H,logg,np.log10(teff)])
          a2 = self._a2([M_H,logg,np.log10(teff)])
          a3 = self._a3([M_H,logg,np.log10(teff)])
          a4 = self._a4([M_H,logg,np.log10(teff)])
        except:
          y,a1,a2,a3,a4  = [np.nan,np.nan,np.nan,np.nan,np.nan]

    return (np.array([a1, a2, a3, a4, y])).flatten()
