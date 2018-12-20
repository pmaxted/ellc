#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--eps", help="genereate .eps file", action="store_true")
args = parser.parse_args()

if args.eps:
  import matplotlib
  matplotlib.use('Agg')

import ellc
import batman
import numpy as np
import matplotlib.pyplot as plt

params = batman.TransitParams()
params.t0 = 0.                      #time of inferior conjunction
params.per = 1.                     #orbital period
params.rp = 0.1                     #planet radius (in units of stellar radii)
params.a = 10.                      #semi-major axis (in units of stellar radii)
params.inc = 90.                    #orbital inclination (in degrees)
params.ecc = 0.1                    #eccentricity
params.w = 60.                      #longitude of periastron (in degrees)
params.u = [0.1, 0.3]               #limb darkening coefficients
params.limb_dark = "quadratic"      #limb darkening model

t = np.linspace(-0.05, 0.05, 1000)

m = batman.TransitModel(params, t, max_err=0.01)   #initializes model
flux_batman = m.light_curve(params)  #calculates light curve

r_1 = 1./params.a 
r_2 = r_1*params.rp
incl = params.inc
ldc_1 = params.u 
ld_1 = 'quad'
sbratio = 0 
f_s = np.sqrt(params.ecc)*np.sin(params.w*np.pi/180.)
f_c = np.sqrt(params.ecc)*np.cos(params.w*np.pi/180.)
flux_ellc_d = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
    ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='sphere',
    grid_1='default',grid_2='default', f_s=f_s, f_c=f_c)
flux_ellc_f = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
    ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='sphere',
    grid_1='fine',grid_2='fine', f_s=f_s, f_c=f_c)
flux_ellc_s = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
    ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='sphere',
    grid_1='sparse',grid_2='sparse', f_s=f_s, f_c=f_c)

fig=plt.figure(1,figsize=(8,12))
plt.subplot(411)
plt.xlim([-0.05,0.05])
plt.ylim([0.988,1.002])
plt.plot(t, flux_batman,color='k')
plt.ylabel("Flux",fontsize=16)

plt.subplot(412)
plt.ylim([-20,20])
plt.xlim([-0.05,0.05])
line_d, = plt.plot(t, 1e6*(flux_ellc_d-flux_batman),label='default')
line_f, = plt.plot(t, 1e6*(flux_ellc_f-flux_batman),'--',label='fine')
line_s, = plt.plot(t, 1e6*(flux_ellc_s-flux_batman),':',label='sparse')
plt.ylabel("Flux difference [ppm]",fontsize=16)
plt.legend(handles=[line_d,line_f,line_s],loc='upper left')

# Polytrope and Roche potential
q = 0.001
flux_ellc_r = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
    q=q,ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='roche', 
    f_s=f_s, f_c=f_c)
flux_ellc_p = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
    q=q,ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='poly1p5', 
    f_s=f_s, f_c=f_c)
plt.subplot(413)
plt.ylim([-20,20])
plt.xlim([-0.05,0.05])
line_p, = plt.plot(t, 1e6*(flux_ellc_p-flux_batman),label='Polytrope')
line_r, = plt.plot(t, 1e6*(flux_ellc_r-flux_batman),'--',label='Roche')
plt.legend(handles=[line_p,line_r],loc='upper left')
plt.ylabel("Flux difference [ppm]",fontsize=16)

# Compare for i<90
del params
params = batman.TransitParams()
params.t0 = 0.                      #time of inferior conjunction
params.per = 1.                     #orbital period
params.rp = 0.1                     #planet radius (in units of stellar radii)
params.a = 10.                      #semi-major axis (in units of stellar radii)
params.inc = 87.                    #orbital inclination (in degrees)
params.ecc = 0.1                    #eccentricity
params.w = 60.                      #longitude of periastron (in degrees)
params.u = [0.1, 0.3]               #limb darkening coefficients
params.limb_dark = "quadratic"      #limb darkening model

del m
m = batman.TransitModel(params, t, max_err=0.01)   #initializes model
flux_batman = m.light_curve(params)  #calculates light curve
r_1 = 1./params.a 
r_2 = r_1*params.rp
incl = params.inc
ldc_1 = params.u 
ld_1 = 'quad'
sbratio = 0 
f_s = np.sqrt(params.ecc)*np.sin(params.w*np.pi/180.)
f_c = np.sqrt(params.ecc)*np.cos(params.w*np.pi/180.)
flux_ellc = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
    ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='sphere',
    grid_1='default',grid_2='default', f_s=f_s, f_c=f_c)
plt.subplot(414)
plt.ylim([-80,80])
plt.xlim([-0.05,0.05])
line, = plt.plot(t, 1e6*(flux_ellc-flux_batman),label='$i=87^{\circ}$')
plt.legend(handles=[line],loc='upper left')
plt.ylabel("Flux difference [ppm]",fontsize=16)
plt.xlabel("Time from central transit",fontsize=16)


plt.tight_layout()
if args.eps:
  fig.savefig("batman.eps")
else:
  plt.show()
