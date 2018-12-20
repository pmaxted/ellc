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

import numpy as np
import matplotlib.pyplot as plt
import ellc

M_1 = 0.43  
M_2 = 0.17
R_1 = 0.0148
R_2 = 0.0214
period = 39.1/60/24
incl = 90.0
vsini_1 = 1148.8
vsini_2 = 601.3
bfac_1 = 0.36
bfac_2 = 1.10
rotfac_1 = vsini_1/(50.57877*R_1/period)
rotfac_2 = vsini_2/(50.57877*R_2/period)
ld_1 = 'lin'
ld_2 = 'lin'
ldc_1 = 0.4
ldc_2 = 0.6
gdc_1 = 0.4
gdc_2 = 0.4

sbratio = 16485/10000  # = T_2/T_1
t_0 = 0.5*period  # .. to get phases consistent with Groot's Fig. 1

a = 4.20944009361 * period**(2./3.) * (M_1+M_2)**(1./3.)
r_1 = R_1/a
r_2 = R_2/a


t1 = np.arange(-0.2,0.8 , 0.001)*period
pro= ellc.lc(t1, radius_1=r_1, radius_2=r_2, sbratio=sbratio, incl=incl,
     q=M_2/M_1, a=a, period=period,t_zero=t_0,bfac_1=bfac_1, bfac_2=bfac_2, 
     ld_1=ld_1,ld_2=ld_2,ldc_1=ldc_1,ldc_2=ldc_2,gdc_1=gdc_1, gdc_2=gdc_2,
     rotfac_1=rotfac_1, rotfac_2=rotfac_2,shape_1='roche',shape_2='roche')
ret= ellc.lc(t1, radius_1=r_1, radius_2=r_2, sbratio=sbratio, incl=incl,
     q=M_2/M_1, a=a, period=period,t_zero=t_0,bfac_1=bfac_1, bfac_2=bfac_2, 
     ld_1=ld_1,ld_2=ld_2,ldc_1=ldc_1,ldc_2=ldc_2,gdc_1=gdc_1, gdc_2=gdc_2,
     rotfac_1=-rotfac_1, rotfac_2=-rotfac_2,shape_1='roche',shape_2='roche')
incl = 89.0
l30 =ellc.lc(t1, radius_1=r_1, radius_2=r_2, sbratio=sbratio, incl=incl,
     q=M_2/M_1, a=a, period=period,t_zero=t_0,bfac_1=bfac_1, bfac_2=bfac_2, 
     ld_1=ld_1,ld_2=ld_2,ldc_1=ldc_1,ldc_2=ldc_2,gdc_1=gdc_1, gdc_2=gdc_2,
     vsini_1=vsini_1, vsini_2=vsini_2, lambda_1=30, lambda_2=30) 
l210=ellc.lc(t1, radius_1=r_1, radius_2=r_2, sbratio=sbratio, incl=incl,
     q=M_2/M_1, a=a, period=period,t_zero=t_0,bfac_1=bfac_1, bfac_2=bfac_2, 
     ld_1=ld_1,ld_2=ld_2,ldc_1=ldc_1,ldc_2=ldc_2,gdc_1=gdc_1, gdc_2=gdc_2,
     vsini_1=vsini_1, vsini_2=vsini_2, lambda_1=210, lambda_2=210) 

fig=plt.figure(1,figsize=(6,5))
plt.subplot(211)
plt.xlim([-0.1,0.6])
plt.ylim([0.55,1.05])
plt.ylabel('Flux')
plt.plot(t1/period,pro,color='darkblue')
plt.plot(t1/period,l30,color='darkgreen',linestyle='--')

plt.subplot(223)
plt.locator_params(axis = 'x', nbins = 4)
plt.xlim([-0.03,0.03])
plt.ylim([-800,1200])
plt.ylabel('R-M effect [PPM]')
plt.xlabel('Phase')
plt.plot(t1/period, 0.5*(pro-ret)*1e6,color='darkblue')
plt.plot(t1/period, 0.5*(l30-l210)*1e6,color='darkgreen',linestyle='--')

plt.subplot(224)
plt.xlim([ 0.47,0.53])
plt.locator_params(axis = 'x', nbins = 4)
plt.ylabel('Flux difference [PPM]')
plt.ylim([-200,200])
plt.xlabel('Phase')
plt.plot(t1/period  , 0.5*(pro-ret)*1e6,color='darkblue')
plt.plot(t1/period, 0.5*(l30-l210)*1e6,color='darkgreen',linestyle='--')

plt.tight_layout()

if args.eps:
  fig.savefig("PhotRM.eps",dpi=400)
else:
  plt.show()

