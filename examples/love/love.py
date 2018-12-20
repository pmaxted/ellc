#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--eps", help="genereate .eps file", action="store_true")
args = parser.parse_args()

if args.eps:
  import matplotlib
  matplotlib.use('Agg')

import ellc
import numpy as np
import matplotlib.pyplot as plt


t = np.linspace(-0.05, 0.05, 1000)

r_1 = 0.15
r_2 = 0.02
incl = 89.5
ldc_1 = [0.1, 0.3]
ld_1 = 'quad'
sbratio = 0 
e = 0.1
om = 60
f_s = np.sqrt(e)*np.sin(om*np.pi/180.)
f_c = np.sqrt(e)*np.cos(om*np.pi/180.)
q = 0.001

flux_0 = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
    ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='sphere',
    grid_1='default',grid_2='default', f_s=f_s, f_c=f_c)

flux_2 = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
    ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='love',hf_2=2.0,q=q,
    grid_1='default',grid_2='default', f_s=f_s, f_c=f_c)

flux_4 = ellc.lc(t,radius_1=r_1, radius_2=r_2,incl=incl,sbratio=sbratio,
    ld_1=ld_1, ldc_1=ldc_1,shape_1='sphere',shape_2='love',hf_2=4.0,q=q,
    grid_1='default',grid_2='default', f_s=f_s, f_c=f_c)

fig=plt.figure(1,figsize=(6,8))
plt.subplot(211)
plt.xlim([-0.05,0.05])
plt.ylim([0.975,1.002])
line_0, = plt.plot(t, flux_0, 'k-',label='sphere')
line_2, = plt.plot(t, flux_2, 'c--',label='h_f = 2')
line_4, = plt.plot(t, flux_4, 'g:', label='h_f = 4')
plt.ylabel("Flux",fontsize=16)
plt.legend(handles=[line_0,line_2,line_4],loc='lower left')

plt.subplot(212)
plt.xlim([-0.05,0.05])
line_2, = plt.plot(t, 1e6*(flux_2-flux_0), 'c--', label='h_f = 2')
line_4, = plt.plot(t, 1e6*(flux_4-flux_0), 'g:', label='h_f = 4')
plt.ylabel("Flux difference [ppm]",fontsize=16)
plt.legend(handles=[line_2,line_4],loc='upper left')

plt.tight_layout()
if args.eps:
  fig.savefig("love.eps")
else:
  plt.show()
