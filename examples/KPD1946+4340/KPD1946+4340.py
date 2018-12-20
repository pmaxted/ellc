#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function,
                            unicode_literals)
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--eps", help="genereate .eps file",action="store_true")
args = parser.parse_args()

if args.eps:
  import matplotlib
  matplotlib.use('Agg')
import ellc
import batman
import numpy as np
import matplotlib.pyplot as plt

lc = np.loadtxt('kplr007975824.dat',
     dtype={'names':('Phase','Flux'), 'formats': ('f8','f8')})

# Reversing the roles of star 1 and 2 cf. Bloemen et al.
period = 0.40375026
K_2 = 164.0
M_1 = 0.59
M_2 = 0.47
incl = 87.14
a_2 = 0.019771142 * K_2 * period #  * np.sqrt(1 - ecc**2)/np.sin(incl*np.pi/180)
q = M_2/M_1
a = (1+q)*a_2
r_1 = 0.013/a
r_2 = 0.214/a
sbratio = 2.86 
gdc_2 = 0.448  # beta_K
ld_1='claret'
ldc_1 = [0.832, -0.681, 0.621, -0.239]
ld_2 = 'claret'
ldc_2 = [0.818, -0.908, 0.755, -0.252] 

bfac_2 = 1.30

flux_r = ellc.lc(lc['Phase']*period, period=period, t_zero=0,
    radius_1=r_1,  radius_2=r_2, incl=incl, \
    sbratio=sbratio, ld_1=ld_1, ldc_1=ldc_1, ld_2=ld_2, ldc_2=ldc_2, \
    gdc_2=gdc_2, bfac_2=bfac_2, a=a, q=q, shape_1='roche',shape_2='roche')

# Roche potential, exact local gravity calculation
flux_e = ellc.lc(lc['Phase']*period, period=period, t_zero=0,
    radius_1=r_1,  radius_2=r_2, incl=incl, \
    sbratio=sbratio, ld_1=ld_1, ldc_1=ldc_1, ld_2=ld_2, ldc_2=ldc_2, \
    gdc_2=gdc_2, bfac_2=bfac_2, a=a, q=q, shape_1='roche',shape_2='roche',
    exact_grav=True)
print("Max. difference between approx. and exact gravity" +
    " calculation = {0:.1f} ppm".format(1e6*np.max(abs(flux_r-flux_e))))

flux_p = ellc.lc(lc['Phase']*period, period=period, t_zero=0,
    radius_1=r_1,  radius_2=r_2, incl=incl, \
    sbratio=sbratio, ld_1=ld_1, ldc_1=ldc_1, ld_2=ld_2, ldc_2=ldc_2, \
    gdc_2=gdc_2, bfac_2=bfac_2, a=a, q=q, shape_1='roche',shape_2='poly1p5')


flux_r = flux_r/np.median(flux_r/lc['Flux'])
flux_p = flux_p/np.median(flux_p/lc['Flux'])

fig=plt.figure(1,figsize=(12,5))
plt.subplot(211)
plt.xlim([-0.4,0.9])
plt.ylim([0.9945,1.0035])
plt.scatter(lc['Phase']   ,lc['Flux'],color='k',s=1)
plt.scatter(lc['Phase']-1,lc['Flux'],color='k',s=1)
plt.plot(lc['Phase']   ,flux_r,color='darkred')
plt.plot(lc['Phase']-1,flux_r,color='darkred')
plt.plot(lc['Phase']   ,flux_p,'--',color='darkblue')
plt.plot(lc['Phase']-1,flux_p,'--',color='darkblue')
plt.xlabel("Phase")
plt.ylabel("Flux")

plt.subplot(223)
plt.xlim([-0.05,0.05])
plt.ylim([0.996,1.001])
plt.scatter(lc['Phase']   ,lc['Flux'],color='k',s=1)
plt.scatter(lc['Phase']-1,lc['Flux'],color='k',s=1)
plt.plot(lc['Phase']   ,flux_r,color='darkred')
plt.plot(lc['Phase']-1,flux_r,color='darkred')
plt.plot(lc['Phase']   ,flux_p,'--',color='darkblue')
plt.plot(lc['Phase']-1,flux_p,'--',color='darkblue')
plt.xlabel("Phase")
plt.ylabel("Flux")
plt.locator_params(axis = 'x', nbins = 4)

plt.subplot(224)
plt.xlim([0.45,0.55])
plt.ylim([0.994,1.00])
plt.scatter(lc['Phase']   ,lc['Flux'],color='k',s=1)
plt.plot(lc['Phase']   ,flux_r,color='darkred')
plt.plot(lc['Phase']   ,flux_p,'--',color='darkblue')
plt.xlabel("Phase")
plt.ylabel("Flux")
plt.locator_params(axis = 'x', nbins = 4)

plt.tight_layout()
if args.eps:
  fig.savefig("KPD1946+4340.eps")
else:
  plt.show()
