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
import matplotlib.pyplot as plt

# Output from arome-1.0.0/examples/f77multiple.f
#
arome_dat = np.loadtxt("arome.dat")

phase  = arome_dat[:,0]
f      = arome_dat[:,4]
v_mean = arome_dat[:,10]

# For testing purposes, calculate R-M effect twice using both star 1 and star 2 in
# the role of the planet.
# For calculation of orbital velocity, assume R_* = 1R_sun

# !/* planet orbital parameters */
# sma = 4.0d0              !/* stellar radii */
# inc = 86.0d0*pi/180.0d0  !/* radian */
a = 4.0
r_star = 1/a
incl = 86.0

# !/* spin-orbit angle */
# lambda = 30.0d0*pi/180.0d0 !/* radian */
lambda_star = 30.0

# !/* limb-darkening */
# c1 = 0.701d0
# c2 = 0.149d0
# c3 = 0.277d0
# c4 =-0.297d0
ld = 'claret'
ldc = [0.701, 0.149, 0.277, -0.297]

# !/* line profile */
# Vsini  = 15.0d0  !/* Vsini */
# Rp     =  0.1d0  !/* radius of the planet */
vsini = 15.0
r_planet = 0.1*r_star

q  = 0.0001  # Nominal mass ratio
# For calculation of period, assume M_* = 1M_sun
period = (a / 4.20944009361)**(3./2.)

t = np.arange(0.18,0.32,0.0005)
phase_1 = (t-0.25)/period + 0.25
f_1 = ellc.lc(t, radius_1 = r_star, radius_2 = r_planet, sbratio = 1e-9, 
    incl=incl, t_zero=0.25, period=period, a=a, q=q, ldc_1=ldc,
    lambda_1=lambda_star, vsini_1=vsini, ld_1='claret',shape_1='sphere')
rv_1,dum = ellc.rv(t, radius_1 = r_star, radius_2 = r_planet, sbratio = 1e-9, 
      incl=incl, t_zero=0.25, period=period, a=a, q=q, ldc_1=ldc,
      lambda_1=lambda_star, vsini_1=vsini, ld_1='claret',shape_1='sphere',
      flux_weighted=True)
rv_orb_1,dum = ellc.rv(t, radius_1 = r_star, radius_2 = r_planet, sbratio = 1e-9, 
      incl=incl, t_zero=0.25, period=period, a=a, q=q, ldc_1=ldc,
      lambda_1=lambda_star, vsini_1=vsini, ld_1='claret',shape_1='sphere',
      flux_weighted=False)
rm_1 = rv_1 - rv_orb_1

# Light travel time correction
solar_radius = 6.957e8
c  = 2.99792458e8
ltt = a*solar_radius/c/86400.0

phase_2 = (t-0.25)/period + 0.25
t_zero = 0.25 + 0.5*period + ltt
f_2 = ellc.lc(t, radius_2 = r_star, radius_1 = r_planet, sbratio = 1e9, 
    incl=incl, t_zero=t_zero, period=period, a=a, q=1/q, ldc_2=ldc,
    lambda_2=lambda_star, vsini_2=vsini, ld_2='claret',shape_2='sphere')
dummy,rv_2 = ellc.rv(t, radius_2 = r_star, radius_1 = r_planet, sbratio = 1e9, 
      incl=incl, t_zero=t_zero, period=period, a=a, q=1/q, ldc_2=ldc,
      lambda_2=lambda_star, vsini_2=vsini, ld_2='claret',shape_2='sphere',
      flux_weighted=True)
dummy,rv_orb_2 = ellc.rv(t, radius_2 = r_star, radius_1 = r_planet, sbratio = 1e9, 
      incl=incl, t_zero=t_zero, period=period, a=a, q=1/q, ldc_2=ldc,
      lambda_2=lambda_star, vsini_2=vsini, ld_2='claret',shape_2='sphere',
      flux_weighted=False)
rm_2 = rv_2 - rv_orb_2

fig=plt.figure(1,figsize=(8,8))

plt.subplot(211)
plt.scatter(phase,1-f,color='k',label='arome')
plt.plot(phase_1,f_1,color='b')
plt.plot(phase_2,f_2,color='g')
plt.xlim([ 0.18,0.32])
plt.xlabel("Phase")
plt.ylabel("Flux")

plt.subplot(212)
plt.scatter(phase,v_mean,color='k',label='arome')
plt.plot(phase_1,rm_1,color='b')
plt.plot(phase_2,rm_2,color='g')
plt.xlim([ 0.18,0.32])
plt.xlabel("Phase")
plt.ylabel("V_RM [km/s]")

plt.tight_layout()

if args.eps:
  fig.savefig("arome.eps")
else:
  plt.show()
