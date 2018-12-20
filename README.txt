====
ellc
====

Generate light curves or radial velocity curves for a binary
star system with the ellc light curve model. [1]_

Installation of this module requires python >2.7 and numpy >1.10.0

To install this module linux-like systems run the following command
$ python setup.py install 
or to install on your system if you have root access
$ sudo  python setup.py install 

Routines in this module:

  lc(t_obs, radius_1, radius_2, sbratio, incl, ... ) 

  rv(t_obs, radius_1, radius_2, sbratio, incl, t_zero, period, a, q, ...)

  fluxes(t_obs, radius_1, radius_2, sbratio, incl, ... ) 
  
  ldy.LimbGravityDarkeningCoeffs()

Documentation and examples for each routine are included in the file headers
and can be viewed from within python in the normal way, e.g., 
 >>> import ellc
 >>> print(ellc.ldy.__doc__)

The easiest way to fit a single light curve is use ellc_emcee (you will need
to add the appropriate /bin/ directory to your path). Run the command
ellc_emcee and follow the prompts. You can grab the output from this script in
the log file and use this as the basis of a new input file for this script.
Use ./ellc_emcee.py -h to see command-line options.

Please cite Maxted (2016) if you published results based on this software.
See sub-directory examples/ for scripts used to generate the figures in this
paper. 

.. rubric:: References
.. [1] Maxted, P.F.L. 2016. A fast, flexible light curve model for detached
   eclipsing binary stars and transiting exoplanets. A&A 591, A111
