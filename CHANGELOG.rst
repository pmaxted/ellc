0.1.0 (2016-01-07)
~~~~~~~~~~~~~~~~~~
 First version of ellc succesfully installed as a python package

0.2.0 (2016-01-08)
~~~~~~~~~~~~~~~~~~
 Added ellc.tmin() 
 Fixed parameter passing for heat_1, lambda_1, etc. in lc.py and rv.py
 Fixed bub in apsidal motion calculation.

0.3.0 (2016-02-02)
~~~~~~~~~~~~~~~~~~
 Added limb-darkening look-up function ldy.py
 Added examples/ellc_emcee/

0.4.0 (2016-02-02)
~~~~~~~~~~~~~~~~~~
 Minor changes.

0.5.0 (2016-02-06)
~~~~~~~~~~~~~~~~~~
 Minor changes to fortran source files to allow compilation with less tolerant
 version of gfortran.

0.6.0 (2016-02-08)
~~~~~~~~~~~~~~~~~~
 Added lrat as a prior and heating coefficients as variables in ellc_emcee.py.
 Fixed bug in calculation of partial eclipses (uninitialised variables)
 Change contents of examples/ellc_emcee/ to analysis of HD23642.

0.7.0 (2016-02-10)
~~~~~~~~~~~~~~~~~~
 Simplified treatment of reflection.

0.7.1 (2016-02-23)
~~~~~~~~~~~~~~~~~~
 Added ldy to inline doc in __init__.py
 Added astropy requirement to setup.py

0.7.2 (2016-03-22)
~~~~~~~~~~~~~~~~~~
 Added extra output to error message if stars exceed Roche lobes
 Changed ellc_emcee.py so that uniform priors are tested before lc() is called.

0.8.0 (2016-03-28)
~~~~~~~~~~~~~~~~~~
 Fixed bug in printed information re: filter names in ellc_emcee.py

0.9.0 (2016-04-19)
~~~~~~~~~~~~~~~~~~
 Removed tmin function.
 Fixed bug in printed information re: filter names in ellc_emcee.py
 Added exact_grav option to lc().
 Improved root polishing in ell_ell_roots function (added do while loop)

1.0.0 (2016-05-11)
~~~~~~~~~~~~~~~~~~
Version sent to A&A

1.1.0 (2016-05-11)
~~~~~~~~~~~~~~~~~~
 Changed definition of third-light contribution to be more intuitive.

1.2.0 (2016-05-11)
~~~~~~~~~~~~~~~~~~
 Added clean exit if polynomial root finding step fails in ell_ell_intersect,
 and raised warning for cases where results may be inaccurate due to
 root-polishing convergence problems.

1.3.0 (2016-05-11)
~~~~~~~~~~~~~~~~~~
 Added test to catch rare case in module ellipse, subroutine ell_ell_intersect
 where test for tangent/intersection point fails because one of the mid-points 
 between the intersection/tangent points is also a tangent point (or nearly 
 so).

1.4.0 (2016-05-24)
~~~~~~~~~~~~~~~~~~
 Added speed-up to ellc.rv() for flux_weight=False.
 Corrected bug in initialisation of star shapes/fluxes for eccentric orbits.
 Fixed crash in lc.py for some choices of limb darkening law.

1.4.1 (2016-06-04)
~~~~~~~~~~~~~~~~~~
 Avoided stop in ell_ell_intersect when ellipse intersection algorithm fails.

1.4.2 (2016-06-09)
~~~~~~~~~~~~~~~~~~
 Avoided problems with NaNs returned in rare cases where the spot is on the
 limb but no tangent points are defined for the projected ellipse.
 Also removed stop for the case ntouch > 4.
 Corrected first line of doc string for rv() (Errors spotted by Bogumil Pilecki)
 Correct bug in ellc.rv() - calculation of rv for iobs=0 can cause crash

1.4.3 (2016-06-21)
~~~~~~~~~~~~~~~~~~
 Corrected bug in calculation of simplified reflection  
 Corrected bug in calculation of flux-weighted radial velocity for stars with
 detailed reflection.
 Added check for flux_weighted radial velocity combined with simplified
 reflection in rv().

1.5.0 (2016-09-26)
~~~~~~~~~~~~~~~~~~
 Revised coordinate systems used internally.  This is needed for a consistent
 fix to a bug in RM effect in 1.4.3 and a bug in the flux-weighted radial
 velocity when using detailed reflection.

 Added module coords.f90

 Added GD448_rv.py to examples/GD448 to test RV+heating calculations.

 Added notes/figures on coordinate systems in subdirectory doc/.

 Fixed bug initialisation of vsini_2 value in lc.py and rv.py. This changes
 the output of examples/PhotRM/PhotRM.py

 Improved calculation of gravity in the exact_grav case - now done on the
 Roche equipotential at the same angular coordinate as the ellipsoid surface.

 Added option to use an improved calculation of ellipsoid size for Roche
 potential case so that volume of the Roche equipotential surface is constant
 in an eccentric orbit and is equal to the volume of a sphere with radius
 specified by the user. With the star model "roche", the volume is calculated
 either the volume of the approximating ellipsoid. For synchronous rotation
 only, the new star shape model "roche_v" uses equation (2.18) from Kopal
 "Dynamics of Close Binary Systems" (Springer, 1978) to calculate the volume
 enclosed by the equipotential surface. This makes very little difference
 unless the star is very distorted and is slower.

 Made printing of parameters python3 compatible in J0113+31.py
 Added bin/ellc_emcee

1.5.2 (2016-11-21)
~~~~~~~~~~~~~~~~~~
 Added ugriz bands to ldy

1.5.3 (2016-11-21)
~~~~~~~~~~~~~~~~~~
 Tidy up output of ldy.list_bands() - fixes crash in ellc_emcee
 Made use of priors in ellc_emcee python3 safe.
 Added "import setuptools" to setup.py

1.6.0 (2017-05-05)
~~~~~~~~~~~~~~~~~~
 - Added "mugrid" limb darkening law.
 - Updated factor for converting K and P to a.sini to be consistent with nominal
   solar constants from IAU Resolution B3.
   (https://www.iau.org/static/resolutions/IAU2015_English.pdf)
 - Changed polynomial root polishing algorithm in ell_ell_roots to Maehly's
   procedure.
 - Fixed bug with eclipses of spots on the limb (iiswitch values swapped).
 - Updated machine constants to double precision in solve_real_poly.f90 
 - Limit correction for eclipsed spots to less than eclipsed flux.

1.6.1 (2017-05-05)
~~~~~~~~~~~~~~~~~~
 - Bug fix for total eclipses with spots introduced in previous version.
 - Added note to in-line documentation regarding simplified reflection and
   geometric albedo.

1.6.3 (2017-11-16)
~~~~~~~~~~~~~~~~~~
 - Added "power-2" limb darkening law
 - Added tolerance within which to set mu values > 1 to mu=1 in stellar.f90 
 - Fixed problems with MANIFEST.in

1.7.0 (2017-12-01)
~~~~~~~~~~~~~~~~~~
 - Added "love" model for star shape from Correia et al.,  2014A&A...570L...5C

1.7.1 (2017-12-01)
~~~~~~~~~~~~~~~~~~
 - Fixed bug in calling of ellc.rv with flux_weighted=True

1.8.0 (2018-02-21)
~~~~~~~~~~~~~~~~~~
 - Made bright() routine more robust for the case radius_1 = radius_2
 - Added iteration limits to while loops in starshape()
 - Added elegent exit from failed calculation of integration limits.
 - Added tolerance to minimum area for numerical integration to avoid 
   numerical errors for near-tangent ellipse intersections.
 - Caught error in ell_ell_overlap due to four very close near-tangent points
 - Made subroutine bright() more robust when working near the limb.
 - Changed default grid size to "default" in rv() and fluxes()

1.8.1 (2018-08-24)
~~~~~~~~~~~~~~~~~~
 - Fixed bug in lc() that caused a crash NaN appears in the array of
   calculated fluxes and verbose=True
 - Removed "external" declrations in gauss_legendre.f90.
 - Fixed indentation error in fluxes.py
 
1.8.2 (2018-12-19)
~~~~~~~~~~~~~~~~~~
  - improved reliability of bright() function for the case exact_grav=True
  - fix for case of exactly overlapping circles.
