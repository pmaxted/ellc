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
flux_q = 10**(-0.4*lc_dat[:,6])  # 
flux_q = flux_q /np.min(flux_q)

t_zero = 0.5
sbratio = 0.0
r_1 = 0.023765
r_2 = 0.169751
incl = 29.3
q = 0.234
ld_1 = 'quad'
ldc_1 = [0.028,  0.264]
ld_2 = 'quad'
ldc_2 = [0.272, 0.567]
shape_1 = 'roche'
shape_2 = 'roche'
heat_2 = [0.6,3.5,0.0]
lc_1 = ellc.lc(phase,t_zero=t_zero, q=q, heat_2 = heat_2, 
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,  
    ld_1=ld_1, ldc_1=ldc_1, ld_2=ld_2, ldc_2=ldc_2, 
    shape_1=shape_1,shape_2=shape_2) 
heat_2 = [1.0 ,1.5, 0.0]
lc_2 = ellc.lc(phase,t_zero=t_zero, q=q, heat_2 = heat_2, 
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,  
    ld_1=ld_1, ldc_1=ldc_1, ld_2=ld_2, ldc_2=ldc_2, 
    shape_1=shape_1,shape_2=shape_2) 
heat_2 =  2.1
lc_3 = ellc.lc(phase,t_zero=t_zero, q=q, heat_2 = heat_2, 
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,  
    ld_1=ld_1, ldc_1=ldc_1, ld_2=ld_2, ldc_2=ldc_2, 
    shape_1=shape_1,shape_2=shape_2) 

fig=plt.figure(1,figsize=(6,4))
line_q,=plt.plot(phase,flux_q,color='k',label='Nightfall')
line_1,=plt.plot(phase,lc_1,color='g',linestyle='--',label='$0.6,3.5$')
line_2,=plt.plot(phase,lc_2,color='b',linestyle=':',label='$1.0,1.5$')
line_3,=plt.plot(phase,lc_3,color='r',linestyle='-.',label='$\\alpha=2.1$')
plt.xlim([-0.25,0.75])
plt.ylim([0.995,1.065])
plt.xlabel("Phase")
plt.ylabel("Flux")
plt.legend(handles=[line_1,line_2,line_3,line_q],loc='upper right')

plt.tight_layout()

if args.eps:
  fig.savefig("GD448.eps")
else:
  plt.show()
