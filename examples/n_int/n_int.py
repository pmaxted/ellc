#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--png", help="genereate .png file", action="store_true")
args = parser.parse_args()

if args.png:
  import matplotlib
  matplotlib.use('Agg')

import numpy as np
import ellc
import matplotlib.pyplot as plt


phase = np.linspace(-0.25,0.75,1001)
lc_0 = ellc.lc(phase,t_zero=0, period=1, q=0.8, 
    radius_1=0.28, radius_2=0.26,incl=89.9,sbratio=0.9,  
    ld_1='quad', ldc_1=[0.3,0.2], ld_2='quad', ldc_2=[0.2,0.1], 
    shape_1='poly1p5',shape_2='poly1p5') 
lc_1 = ellc.lc(phase,t_zero=0, period=1, q=0.8, n_int=4, t_exp=0.05,
    radius_1=0.28, radius_2=0.26,incl=89.9,sbratio=0.9,  
    ld_1='quad', ldc_1=[0.3,0.2], ld_2='quad', ldc_2=[0.2,0.1], 
    shape_1='poly1p5',shape_2='poly1p5') 
n_int_2 = np.zeros(len(phase))
n_int_2[::10] = 1
lc_2 = ellc.lc(phase,t_zero=0, period=1, q=0.8, n_int=n_int_2,
    radius_1=0.28, radius_2=0.26,incl=89.9,sbratio=0.9,  
    ld_1='quad', ldc_1=[0.3,0.2], ld_2='quad', ldc_2=[0.2,0.1], 
    shape_1='poly1p5',shape_2='poly1p5') 


fig=plt.figure(1,figsize=(6,4))
plt.plot(phase,lc_0,label='Original')
plt.plot(phase,lc_1,label='Integrate n_int=4, t_exp=0.05')
plt.plot(phase,lc_2,label='Interpolate 1-in-10')
plt.xlim([-0.25,0.75])
plt.xlabel("Phase")
plt.ylabel("Flux")
plt.legend()

plt.tight_layout()

if args.png:
  fig.savefig("n_int.png")
else:
  plt.show()
