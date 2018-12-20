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


lc_dat = np.loadtxt("NightfallCurve.dat")

phase = lc_dat[:,1]
lc_q = 10**(-0.4*lc_dat[:,6])  # 
lc_q = lc_q /np.min(lc_q)
rv_q = lc_dat[:,15]  #   RV(Secondary)

period = 0.103
t_zero = 0.5*period
time = period*phase
sbratio = 0.0001
r_1 = 0.023765
r_2 = 0.169751
incl = 29.3
q = 0.234
ld_1 = 'quad'
ldc_1 = [0.028,  0.264]
ld_2 = 'quad'
ldc_2 = [0.272, 0.567]
shape = 'roche'
heat = [0.6,3.5,0.0]
dummy,rv_0 = ellc.rv(time,t_zero=t_zero, q=q, heat_2 = heat, 
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,  
    ld_1=ld_1, ldc_1=ldc_1, ld_2=ld_2, ldc_2=ldc_2, a=0.737,
    period=period, shape_1=shape,shape_2=shape,flux_weighted=False) 

lc_1 = ellc.lc(time,t_zero=0.0, q=1/q, heat_1 = heat, 
    radius_1=r_2, radius_2=r_1,incl=incl,sbratio=1./sbratio,  a=0.737,
    period=period, ld_1=ld_2, ldc_1=ldc_2, ld_2=ld_1, ldc_2=ldc_1, 
    shape_1=shape,shape_2=shape) 
lc_1 = lc_1 /np.min(lc_1)
rv_1,dummy = ellc.rv(time,t_zero=0.0, q=1/q, heat_1 = heat, 
    radius_1=r_2, radius_2=r_1,incl=incl,sbratio=1./sbratio,  a=0.737,
    period=period, ld_1=ld_2, ldc_1=ldc_2, ld_2=ld_1, ldc_2=ldc_1, 
    shape_1=shape,shape_2=shape) 

# Test here is the same results are found if role of star 1 and star 2 are
# swapped (get test for coordinate transforms in stellar.f90/bright()
lc_2 = ellc.lc(time,t_zero=t_zero, q=q, heat_2 = heat, 
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,  a=0.737,
    period=period, ld_1=ld_1, ldc_1=ldc_1, ld_2=ld_2, ldc_2=ldc_2, 
    shape_1=shape,shape_2=shape) 
lc_2 = lc_2 /np.min(lc_2)
dummy,rv_2 = ellc.rv(time,t_zero=t_zero, q=q, heat_2 = heat, 
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,  a=0.737,
    period=period, ld_1=ld_1, ldc_1=ldc_1, ld_2=ld_2, ldc_2=ldc_2, 
    shape_1=shape,shape_2=shape) 

fig=plt.figure(1,figsize=(8,8))

plt.subplot(211)
line_q,=plt.plot(phase,lc_q,color='k',label='Nightfall')
line_1,=plt.plot(phase,lc_1,color='g',linestyle='--',label='ellc 1')
line_2,=plt.plot(phase,lc_2,color='r',linestyle='-.',label='ellc 2')
plt.xlim([-0.25,0.75])
plt.xlabel("Phase")
plt.ylabel("Flux")
plt.legend(handles=[line_1,line_2,line_q],loc='upper right')

plt.subplot(212)
line_q,=plt.plot(phase,rv_q,color='k',label='Nightfall')
line_0,=plt.plot(phase,rv_0,color='b',linestyle=':',label='c-of-m')
line_1,=plt.plot(phase,rv_1,color='g',linestyle='--',label='ellc 1')
line_2,=plt.plot(phase,rv_2,color='r',linestyle='-.',label='ellc 2')
plt.xlim([-0.25,0.75])
plt.ylim([-200,200])
plt.xlabel("Phase")
plt.ylabel("RV [km/s]")
plt.legend(handles=[line_0,line_1,line_2,line_q],loc='lower right')

plt.tight_layout()

if args.eps:
  fig.savefig("GD448_rv.eps")
else:
  plt.show()
