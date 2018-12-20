#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--eps", help="generate .eps file",action="store_true")
args = parser.parse_args()

if args.eps:
  import matplotlib
  matplotlib.use('Agg')

import numpy as np
import ellc
import matplotlib.pyplot as plt

lc_dat = np.loadtxt("NightfallCurve.dat")

phase = lc_dat[:,1]
flux = 10**(-0.4*lc_dat[:,4])  # V-band
rv_1 = lc_dat[:,14]
rv_2 = lc_dat[:,15]

# Note that star 1 is the star labelled Secondary in nightfall
t_zero = 0.0
period = 4.0
sbratio = 1.4
r_1 = 0.09992 
r_2 = 0.09726
incl = 89.99
a = 16.826
q = 1/0.8
rotfac_1 = 1
rotfac_2 = 15
shape_1 = 'roche'
shape_2 = 'roche'
gdc_1 = 1
gdc_2 = 1
ld_1 = 'quad'
ldc_1 = [0.283,  0.385] # V-band coeffs from nightfall (using -Dv flag)
ld_2 = 'quad'
ldc_2 = [ 0.417, 0.285] # V-band coeffs from nightfall (using -Dv flag)
ecc = 0  # e>0 does not work in nightfall!
om  = 90
f_c = np.sqrt(ecc)*np.cos(om*np.pi/180.)
f_s = np.sqrt(ecc)*np.sin(om*np.pi/180.)

t = phase*period
lc = ellc.lc(t,t_zero=t_zero, period=period, a=a, q=q,
    ld_1=ld_1,ldc_1=ldc_1,ld_2=ld_2,ldc_2=ldc_2,
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio, 
    rotfac_1=rotfac_1, rotfac_2=rotfac_2, gdc_1=gdc_1, gdc_2=gdc_2,
    shape_1=shape_1,shape_2=shape_2, grid_1='default',grid_2='default') 
rv_ellc_2,rv_ellc_1 = ellc.rv(t,t_zero=t_zero, period=period, a=a, q=q,
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio, 
    rotfac_1=rotfac_1, rotfac_2=rotfac_2, gdc_1=gdc_1, gdc_2=gdc_2,
    shape_1=shape_1,shape_2=shape_2, grid_1='default',grid_2='default') 

fontsize=9
fig=plt.figure(1,figsize=(8,4))
fig=plt.figure(1)
plt.subplot(211)
plt.xlim([-0.25,0.75])
plt.ylim([0.4,1.1])
plt.plot(phase,flux/np.max(flux),color='darkblue',linestyle='--')
plt.plot(phase,lc/np.max(lc),color='darkgreen',linestyle=':')
plt.xlabel("Time [d]",fontsize=fontsize)
plt.ylabel("Flux",fontsize=fontsize)
plt.tick_params(axis='both', labelsize=fontsize)

plt.subplot(223)
plt.xlim([-0.05,0.05])
plt.ylim([-50,50])
plt.plot(phase,rv_1,color='darkblue',linestyle='--')
plt.plot(phase,rv_2,color='darkgreen',linestyle='--')
plt.plot(phase,rv_ellc_1,color='darkblue',linestyle=':' )
plt.plot(phase,rv_ellc_2,color='darkgreen',linestyle=':')
plt.locator_params(axis = 'x', nbins = 4)
plt.tick_params(axis='both', labelsize=fontsize)
plt.xlabel("Time [d]",fontsize=fontsize)
plt.ylabel("Radial velocity [km/s]",fontsize=fontsize)
plt.tick_params(axis='both', labelsize=fontsize)

plt.subplot(224)
plt.xlim([0.45,0.55])
plt.ylim([-400,400])
plt.plot(phase,rv_1,color='darkblue',linestyle='--')
plt.plot(phase,rv_2,color='darkgreen',linestyle='--')
plt.plot(phase,rv_ellc_1,color='darkblue',linestyle=':' )
plt.locator_params(axis = 'x', nbins = 4)
plt.tick_params(axis='both', labelsize=fontsize)
plt.ylabel("Radial velocity [km/s]",fontsize=fontsize)
plt.xlabel("Time [d]",fontsize=fontsize)

plt.tight_layout()
if args.eps:
  fig.savefig("DoublePartial.eps")
else:
  plt.show()

