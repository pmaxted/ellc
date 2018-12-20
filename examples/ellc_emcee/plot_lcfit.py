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
import matplotlib.pyplot as plt
from astropy.table import Table, Column

#------------------------------------------------------------------------------

# Handy class for printing floats...
class prettyfloat(float):
  def __repr__(self):
    return "%0.6f" % self

#------------------------------------------------------------------------------
# Generate plots

model = Table.read('model.csv')
model = Table(model, masked=True)

mask = model['flag'] < 0 

fig=plt.figure(1,figsize=(12,8))

ph_obs = np.extract(mask == False,model['phase']) 
m_obs = np.extract(mask == False,model['mag'])
m_fit = np.extract(mask == False,model['fit'])
m_res = m_obs - m_fit
print ('RMS residual = {0:0.3f}'.format(1000*np.std(m_res)))
ph_plt = np.array(model['phase']) 
f_plt = np.array(model['fit'])
i_sort = np.argsort(ph_plt)
ph_plt = ph_plt[i_sort]
f_plt  = f_plt[i_sort]
fontsize=12

plt.subplot(321)
plt.scatter(ph_obs,m_obs,color='darkgreen',marker='x',s=3)
plt.scatter(ph_obs-1,m_obs,color='darkgreen',marker='x',s=3)
plt.plot(ph_plt,f_plt,color='darkblue')
plt.plot(ph_plt-1,f_plt,color='darkblue')
# plt.scatter(ph_obs,m_res-0.01,color='darkgreen',marker='x',s=3)
# plt.scatter(ph_obs-1,m_res-0.01,color='darkgreen',marker='x',s=3)
# plt.plot([-0.08,0.08],[-0.01,-0.01],':',color='darkblue')
plt.ylabel("Kp [mag]",fontsize=fontsize)
plt.xlim([-0.08,0.08])
plt.ylim([0.12,-0.02])

plt.subplot(322)
plt.scatter(ph_obs,m_obs,color='darkgreen',marker='x',s=3)
# plt.scatter(ph_obs,m_res-0.01,color='darkgreen',marker='x',s=3)
plt.plot(ph_plt,f_plt,color='darkblue')
# plt.plot([ 0.42,0.58],[-0.01,-0.01],':',color='darkblue')
plt.xlim([ 0.42,0.58])
plt.ylim([0.12,-0.02])

plt.subplot(312)
plt.scatter(ph_obs,m_obs,color='darkgreen',marker='x',s=3)
plt.scatter(ph_obs-1,m_obs,color='darkgreen',marker='x',s=3)
plt.plot(ph_plt,f_plt,color='darkblue')
plt.plot(ph_plt-1,f_plt,color='darkblue')
plt.ylabel("Kp [mag]",fontsize=fontsize)
plt.xlim([-0.35,0.85])
plt.ylim([0.12,-0.01])


plt.subplot(313)
plt.scatter(ph_obs,m_res,color='darkgreen',marker='x',s=3)
plt.scatter(ph_obs-1,m_res,color='darkgreen',marker='x',s=3)
plt.plot([-0.35,0.85],[0,0],color='darkblue')
plt.ylabel("O-C [mag]",fontsize=fontsize)
plt.xlabel("Phase",fontsize=fontsize)
plt.xlim([-0.35,0.85])
plt.ylim([0.0035,-0.0035])



if args.eps:
  fig.savefig("lcfit.eps")
else:
  plt.show()



