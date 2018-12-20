#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function,
                            unicode_literals)
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--eps", help="genereate .eps file", action="store_true")
args = parser.parse_args()

if args.eps:
  import matplotlib
  matplotlib.use('Agg')
import numpy as np
import ellc
import emcee
import scipy.optimize as op
import corner
import matplotlib.pyplot as plt


#------------------------------------------------------------------------------

# Handy class for printing floats...
class prettyfloat(float):
  def __repr__(self):
    return "%0.6f" % self

#------------------------------------------------------------------------------

#  Likelihood function
def lnlike(par,
           lc_wasp,lc_nites,lc_oed,lc_byu,lc_kpno,radvel,
           period, t_zero, q, ldc_white, ldc_iband, ldc_jband, 
           n_int_nites, n_int_kpno):

  r_1, r_2, incl, f_s, f_c, sbratio_jband, K_1 = par

# For the WASP data, only need to calculate the light curve in the region of
# the transit. Phase here puts transit at phase 0.5 (for convenince).
  ph_wasp = (1.5 + (((lc_wasp['HJD']-t_zero)/period) % 1)) % 1
  m = np.zeros_like(lc_wasp['HJD'])
  try:
    m_tr = -2.5*np.log10(ellc.lc((lc_wasp['HJD'])[abs(ph_wasp-0.5) < 0.02], 
      radius_1=r_1, radius_2=r_2, incl=incl, sbratio=0,f_s=f_s, f_c=f_c,
      ld_1='quad', ldc_1 = ldc_white, period=period, t_zero=t_zero ,
      grid_1='sparse',grid_2='sparse'))
    m[abs(ph_wasp-0.5) < 0.02] = m_tr
    res = lc_wasp['dmag']-m
    wt = 1./lc_wasp['e_dmag']**2
    zp = np.sum(res*wt)/np.sum(wt)
    chisq_wasp = np.sum((res-zp)**2*wt)
  except:
    chisq_wasp = 1e20
  
  try:
    m = -2.5*np.log10(ellc.lc(lc_nites['HJD'], radius_1=r_1, radius_2=r_2,
      incl=incl, sbratio=0,f_s=f_s, f_c=f_c, ld_1='quad', ldc_1 = ldc_white,
      period=period, t_zero=t_zero, n_int=n_int_nites ,
      grid_1='sparse',grid_2='sparse') )
    res = lc_nites['dmag']-m
    wt = 1./lc_nites['e_dmag']**2
    zp = np.sum(res*wt)/np.sum(wt)
    chisq_nites = np.sum((res-zp)**2*wt)
  except:
    chisq_nites = 1e20
  
  try:
    m = -2.5*np.log10(ellc.lc(lc_oed['HJD'], radius_1=r_1, radius_2=r_2,
      incl=incl, sbratio=0,f_s=f_s, f_c=f_c, ld_1='quad', ldc_1 = ldc_iband,
      period=period, t_zero=t_zero , grid_1='sparse',grid_2='sparse'))
    res = lc_oed['dmag']-m
    wt = 1./lc_oed['e_dmag']**2
    zp = np.sum(res*wt)/np.sum(wt)
    chisq_oed = np.sum((res-zp)**2*wt)
  except:
    chisq_oed = 1e20
  
  try:
    m = -2.5*np.log10(ellc.lc(lc_byu['HJD'], radius_1=r_1, radius_2=r_2,
      incl=incl, sbratio=0,f_s=f_s, f_c=f_c, ld_1='quad', ldc_1 = ldc_iband,
      period=period, t_zero=t_zero, grid_1='sparse',grid_2='sparse' ))
    res = lc_byu['dmag']-m
    wt = 1./lc_byu['e_dmag']**2
    zp = np.sum(res*wt)/np.sum(wt)
    chisq_byu = np.sum((res-zp)**2*wt)
  except:
    chisq_byu = 1e20
    
  # Calculate semi-major axis
  ecc = f_s**2 + f_c**2
  a_1 = 0.019771142 * K_1 * period * np.sqrt(1 - ecc**2)/np.sin(incl*np.pi/180)
  a = (1+1/q)*a_1

  # kpno secondary eclipse - include semi-major axis for light time correction
  try:
    m = -2.5*np.log10(ellc.lc(lc_kpno['HJD'], radius_1=r_1, radius_2=r_2,
      incl=incl, sbratio=sbratio_jband,f_s=f_s, f_c=f_c, ld_2='quad', a=a, 
      ldc_2 = ldc_jband, period=period, t_zero=t_zero, n_int=n_int_kpno ,
      grid_1='sparse',grid_2='sparse'))
    res = lc_kpno['dmag']-m
    wt = 1./lc_kpno['e_dmag']**2
    zp = np.sum(res*wt)/np.sum(wt)
    chisq_kpno = np.sum((res-zp)**2*wt)
  except:
    chisq_kpno = 1e20

  # Radial velocity
  try:
    rv1,rv2 = ellc.rv(radvel['HJD'], radius_1=r_1, radius_2=r_2, incl=incl, 
                 q=q,sbratio=0,f_s=f_s, f_c=f_c, a=a, period=period, 
                 t_zero=t_zero, flux_weighted=False,
                 grid_1='very_sparse',grid_2='very_sparse')
    res = radvel['RV']-rv1
    wt = 1./radvel['e_RV']**2
    zp = np.sum(res*wt)/np.sum(wt)
    chisq_rv = np.sum((res-zp)**2*wt)
  except:
    chisq_rv = 1e20

  chisq = (chisq_wasp + chisq_nites + chisq_byu + chisq_oed + chisq_kpno + 
           chisq_rv)
  print(list(map(prettyfloat,[r_1, r_2, incl, f_s, f_c, sbratio_jband,
    K_1,chisq])))
  return -0.5*chisq
  
#------------------------------------------------------------------------------
# Generate plots

# Load data from ftp://cdsarc.u-strasbg.fr/pub/cats/J/A%2BA/572/A50/
colnames = (str('HJD'),str('dmag'),str('e_dmag'),str('Band'))
lc_wasp  = np.loadtxt('table5.dat',
                      dtype={'names': colnames,
                             'formats': ('f8','f4','f4','S13')})
lc_nites = np.loadtxt('table6.dat',
                      dtype={'names': colnames,
                             'formats': ('f8','f4','f4','S13')})
lc_oed   = np.loadtxt('table7.dat',
                      dtype={'names': colnames,
                             'formats': ('f8','f4','f4','S13')})
lc_byu   = np.loadtxt('table8.dat',
                      dtype={'names': colnames,
                             'formats': ('f8','f4','f4','S13')})
lc_kpno  = np.loadtxt('table9.dat',
                       dtype={'names': colnames,
                              'formats': ('f8','f4','f4','S13')})
colnames = (str('HJD'),str('RV'),str('e_RV'),str('Inst'))
radvel   = np.loadtxt('table10.dat',
                       dtype={'names': colnames,
                              'formats': ('f8','f4','f4','S8')})

period = 14.2769001
t_zero = 2456023.26988


# Interpolate 9/10 points for nites data.
n_int_nites = np.zeros_like(lc_nites['HJD'])
n_int_nites[::10] = 1
# Make sure first/last points after the large break and the last point in the 
# data are calculated and not interpolated
igap = np.argmax((lc_nites['HJD'])[lc_nites['HJD'] < 2455825])
n_int_nites[igap] = 1
n_int_nites[igap + 1] = 1
n_int_nites[-1] = 1
# Interpolate 9/10 points for kpno data.
n_int_kpno = np.zeros_like(lc_kpno['HJD'])
n_int_kpno[::10] = 1
n_int_kpno[-1] = 1

# Quadratic limb darkening coefficients for primary star from
# Claret & Bloemen, 2011A&A...529A..75C, J/A+A/529/A75/table-af
# Teff = 6000, logg=4.0, [Fe/H] = -0.5, Met=F
# Using Kepler passband for  "white-light" WASP and NITES passbands
ldc_white = [0.3532, 0.2606]
ldc_iband = [0.2552, 0.2598]
# Quadratic limb darkening coefficients for secondary star from
# Claret & Bloemen, 2011A&A...529A..75C, J/A+A/529/A75/table-af
# Teff = 4000, logg=5.0, [Fe/H] = -0.5, Met=F
ldc_jband = [0.1076, 0.2857]

# Initial values
e = 0.3098
om = 278.85*np.pi/180 
f_s = np.sqrt(e)*np.sin(om)
f_c = np.sqrt(e)*np.cos(om)
r_1 = 0.0534
r_2 = 0.0081
incl = 89.084
q = 0.1968 
sbratio_jband = 0.346  # From a previous run
K_1 = 15.84
p_0 = [r_1, r_2, incl, f_s, f_c, sbratio_jband, K_1]
x_e = e + 0.0005*np.random.randn(1000)
x_om = om + 1.29*np.random.randn(1000)*np.pi/180
e_f_s = np.std(np.sqrt(x_e)*np.cos(x_om))
e_f_c = np.std(np.sqrt(x_e)*np.sin(x_om))
# Errors on parameters  - used to set size of initial walkers.
e_p = [0.0021, 0.0004, 0.037, e_f_s, e_f_c, 0.013, 0.01]

# Add estimate of jitter to standard errors for radvel data
jitter = 0.02
radvel['e_RV'] = np.sqrt(radvel['e_RV']**2 + jitter**2)

# Adjust WASP photometry errors
lc_wasp['e_dmag'] = lc_wasp['e_dmag']* 0.84

# Affine-invariant MCMC
ndim, nwalkers, nruns = len(p_0), 50, 1000
pos = [p_0 + e_p*np.random.randn(ndim) for i in range(nwalkers)]
print('Starting MCMC parameter search')
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnlike, 
    args=(lc_wasp,lc_nites,lc_oed,lc_byu,lc_kpno,radvel,
          period, t_zero, q, ldc_white, ldc_iband, ldc_jband, 
          n_int_nites, n_int_kpno),threads=4)
sampler.run_mcmc(pos, nruns)
print('Completed MCMC parameter search')

np.save('chain',sampler.chain)

samples = sampler.chain[:, nruns/2:, :].reshape((-1, ndim))
r_1 = np.median(samples[:,0])
r_2 = np.median(samples[:,1])
incl= np.median(samples[:,2])
f_s = np.median(samples[:,3])
f_c = np.median(samples[:,4])
S_J = np.median(samples[:,5])
K_1 = np.median(samples[:,6])

print('r_1 = ',r_1 ,' +/- ',np.std(samples[:,0]))
print('r_2 = ',r_2 ,' +/- ',np.std(samples[:,1]))
print('incl= ',incl,' +/- ',np.std(samples[:,2]))
print('f_s = ',f_s ,' +/- ',np.std(samples[:,3]))
print('f_c = ',f_c ,' +/- ',np.std(samples[:,4]))
print('S_J = ',S_J ,' +/- ',np.std(samples[:,5]))
print('K_1 = ',K_1 ,' +/- ',np.std(samples[:,6]))
acor = sampler.acor
print('r_1  autocorrelation = ',acor[0])
print('r_2  autocorrelation = ',acor[1])
print('incl autocorrelation = ',acor[2])
print('f_s  autocorrelation = ',acor[3])
print('f_c  autocorrelation = ',acor[4])
print('S_J  autocorrelation = ',acor[5])
print('K_1  autocorrelation = ',acor[6])

af = sampler.acceptance_fraction
print('Median acceptance fraction =',np.median(af))
print('Range of  acceptance fraction =',np.min(af),' to ', np.max(af))
# Plot fit to lightcurves
fig=plt.figure(2,figsize=(5,7))

ph_mod = np.arange(-0.024,0.025,0.001)
offstep = 0.08
offset = 0 

ph_wasp = (1 + (((lc_wasp['HJD']-t_zero)/period) % 1)) % 1
plt.scatter(ph_wasp   ,lc_wasp['dmag'],s=1)
plt.scatter(ph_wasp-1,lc_wasp['dmag'],s=1)
plt.xlim([-0.024,0.024])
plt.ylim([0.35,-0.03])
m = -2.5*np.log10(ellc.lc(ph_mod*period+t_zero, radius_1=r_1, radius_2=r_2,
      incl=incl, sbratio=0,f_s=f_s, f_c=f_c, ld_1='quad', ldc_1 = ldc_white,
      period=period, t_zero=t_zero, grid_1='sparse',grid_2='sparse' ))
plt.plot(ph_mod,m)
plt.text(-0.023,-0.019,'WASP')

offset += offstep
ph_nites = (1 + (((lc_nites['HJD']-t_zero)/period) % 1)) % 1
plt.scatter(ph_nites   ,lc_nites['dmag']+offset,s=1)
plt.scatter(ph_nites-1,lc_nites['dmag']+offset,s=1)
m = -2.5*np.log10(ellc.lc(ph_mod*period+t_zero, radius_1=r_1, radius_2=r_2,
      incl=incl, sbratio=0,f_s=f_s, f_c=f_c, ld_1='quad', ldc_1 = ldc_white,
      period=period, t_zero=t_zero , grid_1='sparse',grid_2='sparse'))
plt.plot(ph_mod,m+offset)
plt.text(-0.023,offset-0.019,'NITES')

offset += offstep
ph_oed = (1 + (((lc_oed['HJD']-t_zero)/period) % 1)) % 1
plt.scatter(ph_oed   ,lc_oed['dmag']+offset,s=1)
plt.scatter(ph_oed-1,lc_oed['dmag']+offset,s=1)
m = -2.5*np.log10(ellc.lc(ph_mod*period+t_zero, radius_1=r_1, radius_2=r_2,
      incl=incl, sbratio=0,f_s=f_s, f_c=f_c, ld_1='quad', ldc_1 = ldc_iband,
      period=period, t_zero=t_zero , grid_1='sparse',grid_2='sparse'))
plt.plot(ph_mod,m+offset)
plt.text(-0.023,offset-0.019,'OED')

offset += offstep
ph_byu = (1 + (((lc_byu['HJD']-t_zero)/period) % 1)) % 1
plt.scatter(ph_byu   ,lc_byu['dmag']+offset,s=1)
plt.scatter(ph_byu-1,lc_byu['dmag']+offset,s=1)
m = -2.5*np.log10(ellc.lc(ph_mod*period+t_zero, radius_1=r_1, radius_2=r_2,
      incl=incl, sbratio=0,f_s=f_s, f_c=f_c, ld_1='quad', ldc_1 = ldc_iband,
      period=period, t_zero=t_zero , grid_1='sparse',grid_2='sparse'))
plt.plot(ph_mod,m+offset)
plt.text(-0.023,offset-0.019,'BYU')

offset += offstep
ph_off = 0.533
ph_kpno = (1 + (((lc_kpno['HJD']-t_zero)/period) % 1)) % 1
plt.scatter(ph_kpno-ph_off,lc_kpno['dmag']+offset,s=1)
t = (ph_mod+ph_off)*period+t_zero
m = -2.5*np.log10(ellc.lc(t, radius_1=r_1, radius_2=r_2,
      incl=incl, sbratio=sbratio_jband,f_s=f_s, f_c=f_c, ld_1='quad', 
      ldc_1 = ldc_jband, period=period, t_zero=t_zero ,
      grid_1='sparse',grid_2='sparse'))
plt.plot(ph_mod,m+offset)
plt.text(-0.023,offset-0.019,'KPNO (Phase$-$0.533)')

plt.xlabel("Phase")
plt.ylabel("Magnitude")
plt.locator_params(axis = 'x', nbins = 6)
plt.tight_layout()

if args.eps:
  fig.savefig("J0113+31.eps")
else:
  plt.show()
plt.clf()
plt.close()

# Plot parameter distributions
fig_p = corner.corner(samples, 
  labels=['$r_1$','$r_2$','$i [^\circ]$','$f_c$','$f_s$','$S_J$',
          '$K_1$ [km/s]'],label_kwargs={'fontsize':'18'},
           truths=p_0)
if args.eps:
  fig_p.savefig("corner.eps")
else:
  plt.show()



