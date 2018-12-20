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


lc_dat = np.loadtxt("lc.dat")


rho = 0.4
period = 1
a=(3.75226985055)*((period)**(2./3.))*((rho)**(1./3.))
print("a =",a)
r_1 = 1/a
r_2 = 0.1*r_1
ld_1 = 'quad'
ldc_1 = [0.4,0.3]
t_zero=0.25
rotfac_1 = 1/23.9
spots_1 = [[180-360*t_zero*rotfac_1,120-360*t_zero*rotfac_1],
           [0,0],
           [5.739170477266787, 8.11641272572196],
           [0.5,0.5]]
      
flux = ellc.lc(lc_dat[:,0],t_zero=t_zero, period=period, \
    radius_1=r_1, radius_2=r_2,incl=90,sbratio=0, rotfac_1 = rotfac_1, \
    ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='sphere',\
    grid_1='sparse',grid_2='sparse',spots_1=spots_1)

fontsize=12
fig=plt.figure(1,figsize=(12,8))
plt.subplot(211)
plt.xlim([0,24])
plt.scatter(lc_dat[:,0],lc_dat[:,1],color='darkgreen',marker='x',s=2)
plt.plot(lc_dat[:,0],flux,linewidth=1,color='darkblue')
plt.xlabel("Time [d]",fontsize=fontsize)
plt.ylabel("Normalized flux",fontsize=fontsize)
plt.tick_params(axis='both', labelsize=fontsize)
plt.subplot(234)
plt.plot(lc_dat[:,0]-6.25,flux,linewidth=1,color='darkblue')
plt.scatter(lc_dat[:,0]-6.25,lc_dat[:,1],color='darkgreen',marker='x',s=2)
plt.locator_params(axis = 'x', nbins = 4)
plt.xlabel("Time-6.25 [d]",fontsize=fontsize)
plt.ylabel("Normalized flux",fontsize=fontsize)
plt.xlim([-0.08,0.08])
plt.ylim([0.975,0.995])
plt.tick_params(axis='both', labelsize=fontsize)
plt.subplot(235)
plt.scatter(lc_dat[:,0]-9.25,lc_dat[:,1],color='darkgreen',marker='x',s=2)
plt.plot(lc_dat[:,0]-9.25,flux,linewidth=1,color='darkblue')
plt.locator_params(axis = 'x', nbins = 4)
plt.xlabel("Time-9.25 [d]",fontsize=fontsize)
plt.xlim([-0.08,0.08])
plt.ylim([0.97,0.99])
plt.tick_params(axis='both', labelsize=fontsize)
plt.subplot(236)
plt.scatter(lc_dat[:,0]-14.25,lc_dat[:,1],color='darkgreen',marker='x',s=2)
plt.plot(lc_dat[:,0]-14.25,flux,linewidth=1,color='darkblue')
plt.locator_params(axis = 'x', nbins = 4)
plt.xlabel("Time-14.25 [d]",fontsize=fontsize)
plt.xlim([-0.08,0.08])
plt.ylim([0.98, 1.0])
plt.tick_params(axis='both', labelsize=fontsize)
plt.tight_layout()

if args.eps:
  fig.savefig("KSint.eps")
else:
  plt.show()


