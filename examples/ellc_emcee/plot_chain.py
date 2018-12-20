#!/usr/bin/env python
from __future__ import (absolute_import, division, print_function,
                            unicode_literals)
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--eps", help="genereate .eps file", action="store_true")

parser.add_argument("--burn_in_factor", "-b", 
  default='4',
  help='Number of burn-in steps = BURN_IN_FACTOR x autocorrelation length')

args = parser.parse_args()

if args.eps:
  import matplotlib
  matplotlib.use('Agg')
import numpy as np
import corner
import matplotlib.pyplot as plt
from astropy.table import Table, Column
from emcee.autocorr import function as autocorr_function

#------------------------------------------------------------------------------

# Handy class for printing floats...
class prettyfloat(float):
  def __repr__(self):
    return "%0.6f" % self

def _autocorrelation_length(v):
  """
  Estimate the autocorrelation length of a time series

  Finds the point where the autocorrelation function drops to e^-1
  """
  x = np.arange(len(v))
  acf = autocorr_function(v)
  return np.argmax(x[acf > np.exp(-1)])

#------------------------------------------------------------------------------
# Generate plots


chain = Table.read('chain.csv')
chain = Table(chain, masked=True)

nwalkers = 1+np.max(chain['walker'])
nsteps = 1+np.max(chain['step'])
print('Read chain.csv with {:d} walkers of {:d} steps'
    .format(nwalkers,nsteps))

step = chain['step']
walker = chain['walker']
chain.remove_column('step')
chain.remove_column('walker')
chain.remove_column('loglike')
lrat = chain['sb2']*(chain['r_2']/chain['r_1'])**2
chain.add_column(Column(lrat,name='lrat'))

ndim = len(chain.colnames)
  
acl = np.zeros([ndim,nwalkers],dtype=int)
print('\nCalculating autocorrelation length of chains for each walker.')
print('  For each parameter, results are quoted as median and range.')   
for icol,colname in enumerate(chain.colnames):
  p =  (chain[colname]).reshape([nwalkers,nsteps])
  for ii in range(0,nwalkers):
    acl[icol,ii] = _autocorrelation_length(p[ii,:])
  print(' {:s} : {} ({} - {})'.format(colname,
        np.median(acl[icol,:]),np.min(acl[icol,:]),np.max(acl[icol,:])))

burnin = np.zeros(nwalkers,dtype=int)
for ii in range(0,nwalkers):
  burnin[ii] = np.int(args.burn_in_factor)*np.max(acl[:,ii])
mask = (step == 0)
for ii in range(0,nwalkers):
  mask[walker == ii] = step[walker == ii] < burnin[ii]
nok = len(chain) - len(mask.nonzero()[0])

print ('Calculating results for chains with total of {:d} steps'.format(nok))

print("\n Median and std. dev. of parameter values.")  
for icol,colname in enumerate(chain.colnames):
  p =  np.extract(mask == False,chain[colname])
  if colname == 'T_0' :
    print('{:s} = {:14.6f} +/- {:8.6f}'.format(colname, 
      np.median(p),  np.std(p)))
  elif colname == 'P' :
    print('{:s} = {:12.8f} +/- {:10.8f}'.format(colname, 
      np.median(p),  np.std(p)))
  else:
    print('{:s} = {:10.6g} +/- {:10.6g}'.format(colname, 
      np.median(p),  np.std(p)))

colnames_plot = ['r_1','r_2','incl','sb2','lrat']
labels_plot = ['$r_1$','$r_2$','$i [^{\circ}]$','$s$','$l_2/l_1$']
nplot = len(colnames_plot)
medians = np.zeros(nplot)
samples = np.zeros([nok,nplot])
for ii,colname in enumerate(colnames_plot) :
  p =  np.extract(mask == False,chain[colname])
  samples[:,ii] = p
  medians[ii] = np.median(p)

# Results from Trevor et al., arXiv:1602.01901v1
previous = [0.1450, 0.1262, 78.21, 0.4859, 0.355]
# Plot parameter distributions
fig_p = corner.corner(samples, labels=labels_plot, truths=previous)
if args.eps:
  fig_p.savefig("corner.eps")
else:
  plt.show()



