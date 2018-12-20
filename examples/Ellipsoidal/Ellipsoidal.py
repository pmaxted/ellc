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
from ellc import ldy,lc

import matplotlib.pyplot as plt


lc_dat = np.loadtxt("NightfallCurve.dat")

phase = lc_dat[:,1]
flux_q = 10**(-0.4*lc_dat[:,6])  # 
flux_q = flux_q /np.min(flux_q)

# Star 1 = "Secondary" to that t_zero = 0
sbratio = 1
r_2 = 0.166822
r_1 = 0.007043 
incl = 80.0
q = 0.5
shape_1 = 'sphere'
ldy_ = ldy.LimbGravityDarkeningCoeffs('I')
a1,a2,a3,a4,y = ldy_(20000,4.0,0.0)
# y = 0
print('y =',y)
ld_2 = 'claret'
ldc_2 = [a1,a2,a3,a4]
ecc = 0.15
om  = (180+120) * np.pi/180.
f_c = np.sqrt(ecc)*np.cos(om)
f_s = np.sqrt(ecc)*np.sin(om)
shape_2 = 'roche'
lc_1 = lc(phase,t_zero=0, q=q, 
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,  
    ld_2=ld_2, ldc_2=ldc_2, gdc_2=y,
    f_c=f_c, f_s=f_s,
    shape_1=shape_1,shape_2=shape_2,exact_grav=True) 
lc_2 = lc(phase,t_zero=0, q=q, 
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,  
    ld_2=ld_2, ldc_2=ldc_2, gdc_2=y,
    f_c=f_c, f_s=f_s,
    shape_1=shape_1,shape_2=shape_2,exact_grav=False) 
shape_2 = 'roche_v'
lc_3 = lc(phase,t_zero=0, q=q, 
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,  
    ld_2=ld_2, ldc_2=ldc_2, gdc_2=y,
    f_c=f_c, f_s=f_s,
    shape_1=shape_1,shape_2=shape_2,exact_grav=True) 
lc_4 = lc(phase,t_zero=0, q=q, 
    radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,  
    ld_2=ld_2, ldc_2=ldc_2, gdc_2=y,
    f_c=f_c, f_s=f_s,
    shape_1=shape_1,shape_2=shape_2,exact_grav=False) 

flux_q = flux_q/np.median(flux_q)
lc_1 = lc_1/np.median(lc_1)
lc_2 = lc_2/np.median(lc_2)
lc_3 = lc_3/np.median(lc_3)
lc_4 = lc_4/np.median(lc_4)


fig=plt.figure(1,figsize=(8,6))
line_q,=plt.plot(phase,flux_q,color='k',label='Nightfall')
line_1,=plt.plot(phase,lc_1,color='g',linestyle='--',label='ellc: exact, roche')
line_2,=plt.plot(phase,lc_2,color='b',linestyle=':',label='ellc: poly, roche')
line_3,=plt.plot(phase,lc_3,color='r',linestyle='-.',label='ellc: exact, roche_v')
line_4,=plt.plot(phase,lc_4,color='c',linestyle='-',label='ellc: poly, roche_v')
plt.xlim([-0.25,0.75])
plt.ylim([0.95,1.05])
plt.xlabel("Phase")
plt.ylabel("Flux")
plt.legend(handles=[line_1,line_2,line_3,line_4,line_q],loc='upper right')

plt.tight_layout()

if args.eps:
  fig.savefig("Ellipsoidal.eps")
else:
  plt.show()
