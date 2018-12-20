! This file is part of the ellc binary star model
! Copyright (C) 2016 Pierre Maxted
! 
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module stellar
! Stellar astrophysics functions for ellc
!
! HISTORY
! -------
!  24 Feb 2017
!   Added catch to avoid getting stuck in while loop.
!
! 5 May 2017
!  Added tabulated limb darkening law.
!
! 14 Sep 2016
! Improved definition/handling of coordinates to fix problems with radial
! velocities in bright().
! Improved calculation of local surface gravity - now done at Roche
! equipotential instead of at the ellipsoid surface.
! p.maxted@keele.ac.uk
!
! Apr 2016
! Added exact_grav option to bright
! p.maxted@keele.ac.uk
!
! Jan 2016
! First version 
! p.maxted@keele.ac.uk
!
!
use constants
use utils

implicit none

! Routines for calculating geometry
public starshape     ! Ellipsoidal approximation to shape of a body in a binary
public roche         ! Roche potential from  Wilson, 1979ApJ...234.1054
public roche_l1      ! Position of the inner Lagrangian point
public droche        ! Derivative of Roche potential 
! Stellar surface flux distribution
public limbdark      ! Limb darkening laws
public bright        ! Stellar surface brightness and/or rad. vel.
public ld_quad_match ! Quadratic coefficients to match another limb darking law
public gmodel_coeffs ! Calculate coeffs of surface gravity approximating func.
! Keplerian orbits
public eanom         ! Eccentric anomaly from mean anomaly
public trueanom      ! True anomaly from mean anomaly
public t_ecl_to_peri ! Time of periastron passage prior to a time of eclipse.
public delta_func    ! Helper function for t_ecl_to_peri

contains

function starshape(radius, rsep, spar, ecc, q, model, rtol, &
                   verbose) result(abcd)
implicit none 
! Semi-major axes (A,B,C) and offset towards companion, D, from centre-of-mass
! to geometrical centre of an ellipsoid that approximates the shape of a
! star in a binary system. 
!
! All quantities measured in units of semi-major axis of binary.
! A = semi-major axis on x-axis towards companion
! B = semi-major axis in orbital plane perpendicular to x-axis
! C = semi-major axis on an axis perpendicular to orbital angular momentum
!     vector
!
! Rotation and orbital angular momentum vectors are assumed parallel.
!
! For invalid input all return values are set to bad_dble
!
! Polytropic models use rotational distortion from interpolation in Table 1 of
! James, 1964ApJ...140..552J.
!
! Tidal distortion for polytropes is from Chandrasekhar, 1933MNRAS..93..449C.
!
! Tidal and rotational distortion are assumed to be independent.
!
! The definition of the Roche potential is from Wilson, 1979ApJ...234.1054W
! 
! RADIUS= Radius of a sphere with the same volume as the ellipsoid to be
!         calculated. For the Roche model only, set radius=1 to specify
!         a star that fills its Roche lobe.
! RSEP  = Separation of the stars in units of the semi-major axis.
! SPAR  = angular velocity of star/pseudo-synchronous angular velocity OR
!         for Correia (2014A&A...570L...5C) model, fluid second Love number 
!         for radial displacement, h_f.
! ECC   = eccentricity
! Q     = mass of companion/mass of star
! RTOL  = Tolerance on fractional difference (ABC**(1/3)-radius)/radius.
!         If RTOL <= EPSILON(0.0), e.g., 0.0, then EPSILON(0.0) 
!         (machine precision) is used as the tolerance.
!
! MODEL
!       = -2 - Roche model with precise volume calculation
!       = -1 - Roche model
!       =  0 - spherical star
!       =  1 - polytrope n = 1.5
!       =  2 - polytrope n = 3
!
! VERBOSE : 0 = silent
!           1 = error messages
!           2 = warning messages
!           3 = informational messages
!           4 = debugging
!
!-

double precision, intent(in) :: radius  ! Radius of a sphere with the same volume as the ellipsoid
double precision, intent(in) :: rsep    ! Separation of the stars in units of the semi-major axis
double precision, intent(in) :: spar    ! angular velocity of star/pseudo-synchronous angular velocity
double precision, intent(in) :: ecc     ! eccentricity of the orbit
double precision, intent(in) :: q       ! mass of companion/mass of star
double precision, intent(in) :: rtol    ! Tolerance on fractional difference (ABC**(1/3)-radius)/radius
integer, intent(in)  :: model   ! Model used to calculate star shape (see constants)
integer, intent(in)  :: verbose ! Verbosity level for error reporting
double precision             :: abcd(4) ! semi-major axis on x-axis towards companion

! Local variables
double precision :: rt, t, ttol, tlo, thi
double precision :: a,b,c,d
double precision :: ksi,r_0,vfac,n,n2,r2,r3,r4,r5,r6,q2,q3
double precision :: rsep_l
integer :: it
integer, parameter :: itmax = 100

abcd = (/bad_dble, bad_dble, bad_dble, bad_dble/)

if(verbose >= v_debug) then
  print *,'starshape: radius, q, ecc, rsep, model = ', &
  radius,q,ecc,rsep,model
endif 
if ((radius <= 0.d0).or.(radius.gt.1.d0))then
  if(verbose >= v_error)then
    print *,'starshape: error: invalid input radius: ',radius
  endif
  return
endif

! If rsep is within 4*epsilon(0.d0) of 1-ecc or 1+ecc then assume that this
! is due to round-off error 
! Use rsep_l as a local copy of rsep which is within the allowed range.
tlo = rsep - (1.d0-ecc)
thi = (1.d0+ecc) - rsep
ttol = 4.0d0*epsilon(0.d0)
if ((tlo > 0).and.(thi > 0)) then
  rsep_l = rsep
else if (abs(tlo) < ttol) then
  rsep_l = (1.d0-ecc)
else if (abs(thi) < ttol) then
  rsep_l = (1.d0+ecc)
else
  if(verbose >= v_error)then
    print *,'starshape: error: rsep out of range: ', rsep,ecc
  endif
  return
endif

ttol = max(rtol, epsilon(0.d0))
rt = radius

select case (model)

case(starshape_roche_v)

  t = 1
  it = 0 
  do while (t.gt.ttol)
    it = it + 1
    if (it > itmax) then
      if(verbose >= v_error) then
        print *,'starshape: failed to converge : t  = ',t
      endif
      return
    endif

    if(d.eq.bad_dble) return
    if (radius.eq.1.d0) then
      t = 0
    else
      ksi = rsep_l*roche(x=rt,y=0.0d0,z=0.0d0,d=rsep_l,q=q,f=1.0d0)
      n = (q+1d0)/2d0
      r_0 = 1d0/(ksi-q)
      r2 = r_0**2
      r3 = r2*r_0
      r4 = r3*r_0
      r5 = r4*r_0
      r6 = r5*r_0
      q2 = q**2
      q3 = q2*q
      n2 = n**2
      vfac = (1d0 + r6*(12d0*q2/5d0  &
      + 15d0*q2*r2/7d0 + 2d0*q2*r4 + 22d0*q3*r3/7d0 + 157*q3*r5/7d0 &
      + 2d0*n/r3 + 32d0*n2/5d0 + 176*n**3*r3/7d0 + 8d0*n*q/5d0 &
      + 296*n*q*(2d0*q+n)*r3/35d0 + 26d0*n*q*(q+3d0*n)*r5/35d0 ))
      t = vfac**third*r_0*rsep_l/radius
      rt = rt/t
      t = abs(t-1.d0)
    endif
  end do
  call shape_roche(rt,rsep_l,spar,q,a,b,c,d,verbose-1)

case(starshape_roche)

  t = 1
  it = 0
  do while (t.gt.ttol)
    it = it + 1
    if (it > itmax) then
      if(verbose >= v_error) then
        print *,'starshape: failed to converge : t  = ',t
      endif
      return
    endif

    call shape_roche(rt,rsep_l,spar,q,a,b,c,d,verbose-1)
    if(d.eq.bad_dble) return
    if (radius.eq.1.d0) then
      t = 0
    else
      t = (a*b*c)**third/radius
      rt = rt/t
      t = abs(t-1.d0)
    endif
  end do

case(starshape_sphere)

  a = radius
  b = radius
  c = radius
  d = 0
  t = 0

case(starshape_poly1p5)

  t = 1
  it = 0
  do while (t.gt.ttol)
    it = it + 1
    if (it > itmax) then
      if(verbose >= v_error) then
        print *,'starshape: failed to converge : t, radius  = ',t,radius
      endif
      return
    endif
    call shape_n1p5(rt,rsep_l,spar,ecc,q,a,b,c,d,verbose-1)
    if(d.eq.bad_dble) return
    t = (a*b*c)**third/radius
    rt = rt/t
    t = abs(t-1.d0)
  end do

case(starshape_poly3p0)

  t = 1
  it = 0
  do while (t.gt.ttol)
    it = it + 1
    if (it > itmax) then
      if(verbose >= v_error) then
        print *,'starshape: failed to converge : t  = ',t
      endif
      return
    endif
    call shape_n3p0(rt,rsep_l,spar,ecc,q,a,b,c,d,verbose-1)
    if(d.eq.bad_dble) return
    t = (a*b*c)**third/radius
    rt = rt/t
    t = abs(t-1.d0)
  end do


case(starshape_love)

  t = 1
  it = 0
  do while (t.gt.ttol)
    it = it + 1
    if (it > itmax) then
      if(verbose >= v_error) then
        print *,'starshape: failed to converge : t  = ',t
      endif
      return
    endif
    call shape_love(rt,rsep_l,spar,q,a,b,c,d,verbose-1)
    if(d.eq.bad_dble) return
    t = (a*b*c)**third/radius
    rt = rt/t
    t = abs(t-1.d0)
  end do

case default
  if(verbose >= v_error) then
    print *,'starshape: error: invalid input model'
  endif
  return

end select

if(verbose >= v_debug) then
  print *,'starshape: normal exit : a,b,c,d,t  = ',a,b,c,d,t
endif

abcd = (/a, b, c, d/)
return

end function starshape

!-----------------------------------------------------------------------

subroutine shape_n1p5(radius,rsep,frot,ecc,qmass,a,b,c,d,verbose)
implicit none
double precision, intent(in)  :: radius, rsep, frot, ecc, qmass
double precision, intent(out) :: a, b, c, d
integer, intent(in)           :: verbose

! Local variables
integer, parameter :: npar = 1
double precision   :: ar, par(npar), a1, a2, w
double precision   :: dxip,dxie
double precision   :: dsig0
double precision   :: qnu3,qnu4,qnu5
integer            :: ii
double precision, parameter :: delta2=1.2892d0
double precision, parameter :: delta3=1.1079d0
double precision, parameter :: delta4=1.0562d0
double precision, parameter :: tol = 1.0d-6

! These fractional changes in polar and equatorial radii and 
! fraction change in volume equivalent radius based on columns xi_p, xi_e and
! 10^-2V in Table 1 of James, 1964
!
double precision tdxip(25),tdxie(25)
data tdxip/ &
 0.0000000d0,-0.0042422d0,-0.0084843d0,-0.0127265d0,-0.0169960d0, &
-0.0212929d0,-0.0255898d0,-0.0298867d0,-0.0342383d0,-0.0385900d0, &
-0.0429690d0,-0.0473753d0,-0.0518091d0,-0.0562976d0,-0.0607860d0, &
-0.0653293d0,-0.0699272d0,-0.0745525d0,-0.0792510d0,-0.0840221d0, &
-0.0888390d0,-0.0927801d0,-0.0947780d0,-0.0968033d0,-0.0978160d0/ 
data tdxie/ &
0.00000000d0,0.00747167d0,0.01524440d0,0.02331819d0,0.03174777d0, &
0.04056051d0,0.04981116d0,0.05949970d0,0.06970825d0,0.08046417d0, &
0.09190432d0,0.10408342d0,0.11713832d0,0.13117850d0,0.14642290d0, &
0.16311785d0,0.18161914d0,0.20250151d0,0.22653128d0,0.25521375d0, &
0.29153210d0,0.33124418d0,0.35885927d0,0.40092506d0,0.45073622d0/

par(1) = frot**2 * radius**3 * (1.d0+qmass) * (1.d0-ecc**2) / (1.d0-ecc)**4

a1 = 0.d0
a2 = 1.09d-2
ar =  zbrent(func_n1p5,a1,a2,tol,npar,par,verbose-1)
if(verbose >= v_debug) print *,'shape_n1p5: ar =',ar
if (ar < 0.d0) then
  a = bad_dble
  b = bad_dble
  c = bad_dble
  d = bad_dble
  if(verbose >= v_error) print *,'shape_n1p5: radius out of range.',radius,frot,ecc,qmass
  return
endif

! Interpolate values of fractional changes in radius
if(ar < 1.0d-2) then
  w = 0.2d4*ar
  ii = 1+int(w)
  w = dble(ii)-w
  dxip = w*tdxip(ii) + (1.d0-w)*tdxip(ii+1)
  dxie = w*tdxie(ii) + (1.d0-w)*tdxie(ii+1)
else if(ar < 1.04d-2) then
  w = 0.25d4*(ar-1.0d-2)
  dxie = (1.d0-w)*tdxie(21) + w*tdxie(22)
  dxip = (1.d0-w)*tdxip(21) + w*tdxip(22)
else if(ar <= 1.08d-2) then
  w = 0.5d4*(ar-1.04d-2)
  ii = 22+int(w)
  w = dble(ii-21)-w
  dxie = w*tdxie(ii) + (1.d0-w)*tdxie(ii+1)
  dxip = w*tdxip(ii) + (1.d0-w)*tdxip(ii+1)
else if(ar <= 1.09d-2) then
  w = ar-1.08d-2
  dxie = (1.d0-w)*tdxie(24) + w*tdxie(25)
  dxip = (1.d0-w)*tdxip(24) + w*tdxip(25)
else
  print *,'shape_n1p5: invalid value of ar',ar
  stop
endif

a = radius*(1.d0+dxie)
b = a
c = radius*(1.d0+dxip)

! Tidal distortion
ar = radius/rsep
qnu3 = qmass*ar**3 
qnu4 = qnu3*ar
qnu5 = qnu4*ar
a = a*(1.d0 + 0.5d0*(delta2*qnu3+delta4*qnu5))
dsig0 = 1.d0 - 0.5d0*delta2*qnu3 + 0.375d0*delta4*qnu5
b = b*dsig0
c = c*dsig0
d = delta3*qnu4*radius

return

end subroutine shape_n1p5

!-----------------------------------------------------------------------

double precision function func_n1p5(a,npar,par,verbose)
implicit none
!
! par(1) = f^2 * (1+q) * r_v**3 * (1.d0-ecc**2)/(1.d0-ecc)**4
! 
integer, intent(in)          :: npar, verbose
double precision, intent(in) :: a, par(npar)

! Local variables
double precision array(25)
integer ii
double precision w

! The values in this array are 2*pi*m/3*v where m and v are from Table 1 of
! James, 1964ApJ...140..552J. Return value is linear interpolation in this
! table multiplied by par(1) then subtracted from a
! Out-of range values returned as bad_dble
data array/ &
0.02782105d0,0.02768557d0,0.02754373d0,0.02739694d0,0.02724204d0, &
0.02708024d0,0.02691018d0,0.02673288d0,0.02654358d0,0.02634414d0, &
0.02613295d0,0.02590672d0,0.02566646d0,0.02540777d0,0.02512701d0, &
0.02482182d0,0.02448562d0,0.02410966d0,0.02368468d0,0.02319031d0, &
0.02259155d0,0.02198336d0,0.02160056d0,0.02109994d0,0.02071318d0/

if (verbose > v_debug) print *,'func_n1p5: a =', a

if(a < 0.d0) then
  func_n1p5 = bad_dble
else if(a < 1.0d-2) then
  w = a/0.05d-2
  ii = 1+int(w)
  w = dble(ii)-w
  func_n1p5 = a-par(1)*(w*array(ii) + (1.d0-w)*array(ii+1))
else if(a < 1.04d-2) then
  w = (a-1.0d-2)/0.04d-2
  func_n1p5 = a-par(1)*((1.d0-w)*array(21) + w*array(22))
else if(a <= 1.08d-2) then
  w = (a-1.04d-2)/0.02d-2
  ii = 22+int(w)
  w = dble(ii-21)-w
  func_n1p5 = a-par(1)*(w*array(ii) + (1.d0-w)*array(ii+1))
else if(a <= 1.09d-2) then
  w = a-1.08d-2
  func_n1p5 = a-par(1)*((1.d0-w)*array(24) + w*array(25))
else
  func_n1p5 = bad_dble
endif
return

end function func_n1p5

!-----------------------------------------------------------------------

subroutine shape_n3p0(radius,rsep,frot,ecc,qmass,a,b,c,d,verbose)
implicit none
double precision, intent(in)  :: radius, rsep, frot, ecc, qmass
double precision, intent(out) :: a, b, c, d
integer, intent(in)           :: verbose

! Local variables
integer, parameter :: npar = 1
double precision   :: ar, par(npar), a1, a2, w
double precision   :: dxip,dxie
double precision   :: dsig0
double precision   :: qnu3,qnu4,qnu5
integer            :: ii
double precision, parameter :: delta2=1.0289d0
double precision, parameter :: delta3=1.00736d0
double precision, parameter :: delta4=1.00281d0
double precision, parameter :: tol = 1.0d-6

! These fractional changes in polar and equatorial radii and 
! fraction change in volume equivalent radius based on columns xi_p, xi_e and
! 10^-2v in table 1 of James, 1964
!
double precision tdxip(24),tdxie(24)
data tdxip/ &
 0.0000000d0,-0.0015862d0,-0.0031623d0,-0.0047471d0,-0.0063333d0, &
-0.0079138d0,-0.0094956d0,-0.0110731d0,-0.0126667d0,-0.0142500d0, &
-0.0158232d0,-0.0174036d0,-0.0189855d0,-0.0205674d0,-0.0221536d0, &
-0.0237558d0,-0.0253406d0,-0.0269311d0,-0.0285232d0,-0.0291582d0, &
-0.0297919d0,-0.0301109d0,-0.0304298d0,-0.0307503d0/
data tdxie/ &
0.00000000d0,0.00692635d0,0.01416734d0,0.02175631d0,0.02972951d0, &
0.03812610d0,0.04700117d0,0.05640546d0,0.06641728d0,0.07712796d0, &
0.08861727d0,0.10106498d0,0.11461754d0,0.12951130d0,0.14608118d0, &
0.16478972d0,0.18633434d0,0.21185324d0,0.24351697d0,0.25880946d0, &
0.27637255d0,0.28627707d0,0.29714145d0,0.30924118d0/

par(1)=frot**2 * radius**3 * (1.d0+qmass) * (1.d0-ecc**2)/(1.d0-ecc)**4
a1 = 0.d0
a2 = 9.7d-4
ar =  zbrent(func_n3p0,a1,a2,tol,npar,par,verbose-1)
if(verbose >= v_debug) print *,'shape_n3p0: ar =',ar
if (ar < 0.d0) then
  a = bad_dble
  b = bad_dble
  c = bad_dble
  d = bad_dble
  if(verbose >= 1) print *,'shape_n3p0: radius out of range.',radius,frot,ecc,qmass
  return
endif

! Interpolate values of fractional changes in radius
if(ar.lt.9.0d-4) then
  w = ar/0.5d-4
  ii = 1+int(w)
  w = dble(ii)-w
  dxip = w*tdxip(ii) + (1.d0-w)*tdxip(ii+1)
  dxie = w*tdxie(ii) + (1.d0-w)*tdxie(ii+1)
else if(ar.lt. 9.4d-4) then
  w = (a-9.0d-4)/0.2d-4
  ii = 19+int(w)
  w = dble(ii-18)-w
  dxip = w*tdxip(ii) + (1.d0-w)*tdxip(ii+1)
  dxie = w*tdxie(ii) + (1.d0-w)*tdxie(ii+1)
else if(ar <= 9.7d-4) then
  w = (a-9.4d-4)/0.1d-4
  ii = 21+int(w)
  w = dble(ii-20)-w
  dxip = w*tdxip(ii) + (1.d0-w)*tdxip(ii+1)
  dxie = w*tdxie(ii) + (1.d0-w)*tdxie(ii+1)
else
  print *,'shape_n3p0: invalid value of ar',ar
  stop
endif
a = radius*(1.d0+dxie)
b = a
c = radius*(1.d0+dxip)

! Tidal distortion
ar = radius/rsep
qnu3 = qmass*ar**3 
qnu4 = qnu3*ar
qnu5 = qnu4*ar
a = a*(1.d0 + 0.5d0*(delta2*qnu3+delta4*qnu5))
dsig0 = 1.d0 - 0.5d0*delta2*qnu3 + 0.375d0*delta4*qnu5
b = b*dsig0
c = c*dsig0
d = delta3*qnu4*radius

return

end subroutine shape_n3p0

!-----------------------------------------------------------------------

double precision function func_n3p0(a,npar,par,verbose)
implicit none
!
! par(1) = f^2 * (1+q) * r_v**3 * (1.d0-ecc**2)/(1.d0-ecc)**4
! 
integer, intent(in)          :: npar, verbose
double precision, intent(in) :: a, par(npar)

! Local variables
double precision array(24)
integer ii
double precision w
! The values in this array are 2*pi*m/3*v where m and v are from table 1 of
! james, 1964apj...140..552j. Return value is linear interpolation in this
! table multiplied by par(1) then subtracted from a
! out-of range values returned as bad_dble
data array/ &
0.00307603d0,0.00304253d0,0.00300812d0,0.00297269d0,0.00293619d0, &
0.00289851d0,0.00285958d0,0.00281926d0,0.00277739d0,0.00273384d0, &
0.00268838d0,0.00264078d0,0.00259069d0,0.00253778d0,0.00248146d0, &
0.00242103d0,0.00235546d0,0.00228318d0,0.00220144d0,0.00216504d0, &
0.00212572d0,0.00210468d0,0.00208254d0,0.00205903d0/

if (verbose > v_debug) print *,'func_n1p5: a =', a

if(a < 0.d0) then
  func_n3p0 = bad_dble
else if(a < 9.0d-4) then
  w = a/0.5d-4
  ii = 1+int(w)
  w = dble(ii)-w
  func_n3p0 = a-par(1)*(w*array(ii) + (1.d0-w)*array(ii+1))
else if(a < 9.4d-4) then
  w = (a-9.0d-4)/0.2d-4
  ii = 19+int(w)
  w = dble(ii-18)-w
  func_n3p0 = a-par(1)*(w*array(ii) + (1.d0-w)*array(ii+1))
else if(a <=  9.7d-4) then
  w = (a-9.4d-4)/0.1d-4
  ii = 21+int(w)
  w = dble(ii-20)-w
  func_n3p0 = a-par(1)*(w*array(ii) + (1.d0-w)*array(ii+1))
else
  func_n3p0 = bad_dble
endif
return
end function func_n3p0

!-----------------------------------------------------------------------

subroutine shape_roche(radius,rsep,frot,qmass,a,b,c,d,verbose)
implicit none

double precision, intent(in)  :: radius, rsep, frot, qmass
double precision, intent(out) :: a, b, c, d
integer, intent(in)           ::  verbose

! Local variables
double precision   :: pot, rl1, xf, xb
integer, parameter :: npar = 4
double precision   :: par(npar)
double precision, parameter :: tol=1.0d-6

rl1 = roche_l1(qmass,frot)
if (radius == 1.d0) then
  xf = rl1
else
  if (radius > rl1) then
    if(verbose >= v_error) then
      print *,' shape_roche: star exceeds roche lobe'
      print *,' shape_roche: radius,rl1,q = ',radius,rl1,qmass
    endif
    a = bad_dble
    b = bad_dble
    c = bad_dble
    d = bad_dble
    return
  endif
  xf = radius
endif
pot = roche(xf,0.d0,0.d0,rsep,qmass,frot)
par(1) = pot
par(2) = rsep
par(3) = qmass
par(4) = frot
b = zbrent(func_b,radius*0.5d0,radius*1.1d0,tol,npar,par,verbose-1)
c = zbrent(func_c,radius*0.5d0,radius*1.1d0,tol,npar,par,verbose-1)
xb = zbrent(func_a,-radius*1.5d0,-radius*0.5d0,tol,npar,par,verbose-1)
d = 0.5*(xf+xb)
a = xf - d

end subroutine shape_roche

!-----------------------------------------------------------------------

subroutine shape_love(radius,rsep,h_f,qmass,a,b,c,d,verbose)
implicit none

double precision, intent(in)  :: radius, rsep, h_f, qmass
double precision, intent(out) :: a, b, c, d
integer, intent(in)           :: verbose
double precision              :: qr

! RADIUS = Radius of a spherical planet in unit of stellar radii.
! RSEP   = Separation of between bodies in unit of stellar radii.
! qmass  = mass of star/mass of planet
! h_f    = fluid second love number as used in Correia()
! a,b,c  = the semi-major axes of the ellipse in units of stellar radii
! d      = d is an offset towards the companion from the centre-of-mass of the
!           planet to the geometric centre of the ellipsoid
!
! Local variables.
!qr      = planet asymmetry as defined in correia


b=radius
qr=0.5d0*h_f*qmass*(b/rsep)**3       !from eqn(12) in correia(2014)
if (qr >= 1) then
    if (verbose > v_silent) print *,'shape_love: qr >1; b, qr  =', b, qr
    a = bad_dble
    b = bad_dble
    c = bad_dble
    d = bad_dble
    return
endif

a = b*(1.d0+3.d0*qr)          !eqn(10)
c = b*(1.d0-qr)               !eqn(11)
d = radius*qmass*radius**4    !From Chandrasekhar using delta_3 = 1

if (verbose > v_debug) print *,'shape_love: a,c,d,qr =', a,c,d,qr
end subroutine shape_love

!-----------------------------------------------------------------------

double precision function func_a(a,npar,par,verbose)
implicit none
integer, intent(in)          :: npar, verbose
double precision, intent(in) :: a, par(npar)
if (verbose > v_debug) print *,'func_a: ', a
func_a=par(1)-roche(a,0.d0,0.d0,par(2),par(3),par(4))
return 
end function func_a

!-----------------------------------------------------------------------

double precision function func_b(b,npar,par,verbose)
implicit none
integer, intent(in)          :: npar, verbose
double precision, intent(in) :: b, par(npar)
if (verbose > v_debug) print *,'func_b: ', b
func_b = par(1)-roche(0.d0,b,0.d0,par(2),par(3),par(4))
return 
end function func_b

!-----------------------------------------------------------------------

double precision function func_c(c,npar,par,verbose)
implicit none
integer, intent(in)          :: npar, verbose
double precision, intent(in) :: c, par(npar)
if (verbose > v_debug) print *,'func_c: ', c
func_c=par(1)-roche(0.d0,0.d0,c,par(2),par(3),par(4))
return 
end function func_c

!-----------------------------------------------------------------------

double precision function func_r(r,npar,par,verbose)
implicit none
integer, intent(in)          :: npar, verbose
double precision, intent(in) :: r, par(npar)
if (verbose > v_debug) print *,'func_r: ', r
func_r=par(1)-roche(par(2)*r,par(3)*r,par(4)*r,par(5),par(6),par(7))
!print *,'func_r (',real(r),') = ',real(func_r),', target =',real(par(1))
!print *,'func_r',real(r),real(par(1)-func_r),real(par(1))
return 
end function func_r

!-----------------------------------------------------------------------

double precision function roche(x, y, z, d, q, f, iscomp)
!+
!
! Roche potential according to the definition of Wilson, 1979ApJ...234.1054W
!
! Input position is (X, Y, Z)
!  - origin at centre-of-mass of the star.
!  - x-axis towards companion
!  - y-axis in orbtial plane opposite to direction of orbital motion
!  - z-axis parallel to orbital angular momentum vector
!  - all lengths relative to the semi-major axis of the binary.
! 
! Q is mass ratio = mass companion/mass of star.
!
! D is instantaneous separation of the centres.
!
! F is the synchronization index of the star, i.e. the ratio of sidereal 
! rotational and orbital frequency
!
! To evaluate the potential for the companion to a star with consistent
! units for two stars:
!  - set "ISCOMP" to .TRUE.
!  - Use the value of F appropriate for the companion
!  - Measure position from the centre of the companion with x-axis
!    towards the star
!  - Use mass ratio Q = mass of companion/mass star, as before.
!  
! Invalid input returns the values bad_dble
!
! p.maxted@keele.ac.uk, Jun 2015 
!
!-
implicit none
double precision, intent(in)  :: q       ! mass ratio = mass of companion/mass of star
double precision, intent(in)  :: x, y, z ! position (x, y, z) at which to evaluate the function
double precision, intent(in)  :: d       ! instantaneous separation of the stars
double precision, intent(in)  :: f       ! synchronization index of the star
logical, intent(in), optional :: iscomp  ! Set .true. to calculate potential due to companion star

! Local variables
double precision qq, r, r2, nu, lam, d2
logical iscomp_l

if (present(iscomp)) then
  iscomp_l = iscomp
else
  iscomp_l = .false.
endif

if (q <= 0.d0) then
  roche = bad_dble
  return
endif

if (iscomp_l) then
  qq = 1.d0/q
else
  qq = q
endif

r2 = x**2 + y**2 + z**2
if (r2 == 0.d0) then
  roche = bad_dble
  return
endif

r = sqrt(r2)
nu = z/r
lam = x/r
d2 = d**2
roche = 1.d0/r + qq*( 1.d0/sqrt(d2+r2-2.0d0*r*lam*d) - r*lam/d2 ) &
+ 0.5d0 * f**2 * (1.d0+qq) * r2 * (1.d0-nu**2) 
if (iscomp_l) roche = roche/qq + 0.5d0*(qq-1.d0)/qq

return
end function roche

!-----------------------------------------------------------------------

double precision function roche_l1(q, f)
implicit none
double precision, intent(in) ::  q   ! (mass of companion) / (mass of star)  
double precision, intent(in),optional ::  f   ! Rotation factor
!+
!
! Position of the inner Lagrangian point in terms of the separation of 
! two stars which have mass ratio, q.
! Works by solving for root of a quintic polynomial, by Newton-Raphson
! iteration
! Q = (mass of companion) / (mass of star)
! L1 = (distance from centre of star to inner Lagrangian point)/separation
!
! Optional input factor f is the asynchronous rotation factor as defined by
! Wilson, 1979ApJ...234.1054W
!
! For invalid input return value is bad_dble
!
!-
double complex   ::  az(6) = (/ &
   dcmplx(-1.d0,0.0d0),dcmplx(2.0d0,0.0d0),dcmplx(-1.d0,0.0d0), &
   dcmplx(0.d0,0.0d0),dcmplx(0.0d0,0.0d0),dcmplx(0.d0,0.0d0) /)
double complex   ::  z
integer ::  m
double precision :: fac

if(q <= 0.d0) then
  roche_l1 = bad_dble
  return
endif

if (present(f)) then 
  fac = f**2*(1.0d0+q)
  az(4) =dcmplx( fac + 2.0d0*q,0.0d0)
  az(5) =dcmplx(-2.0d0*fac - q,0.0d0)
  az(6) =dcmplx( fac          ,0.0d0)
else
  az(4) =dcmplx( 1.0d0 + 3.0d0*q,0.0d0)
  az(5) =dcmplx(-2.0d0 - 3.0d0*q,0.0d0)
  az(6) =dcmplx( 1.0d0 + q      ,0.0d0)
endif
z     =dcmplx( 1.0d0/(1.0d0+q),0.0d0)
m=5
roche_l1 = dble(laguerre(az,m,z))
return

end function roche_l1

!-----------------------------------------------------------------------

double precision function droche(q, x, y, z, d, f)
!+
!  Evaluates magnitude of the derivative of the dimensionless Roche potential
! at a point X,Y,Z measured from the centre of mass of the star in units of
! the semi-major axis of the binary.
!
! D = instantaneous separation of the stars. 
!
! Q = mass ratio = mass of companion/mass of star.
!
! F = synchronization index of the star, i.e. the ratio of sidereal 
!     rotational and orbital frequency
! 
! Roche potential according to the definition of Wilson, 1979ApJ...234.1054W
!
! Invalid input returns the values bad_dble
!
! p.maxted@keele.ac.uk, Jun 2015 
!
!-
implicit none
double precision, intent(in) :: q       ! mass ratio = mass of companion/mass of star
double precision, intent(in) :: x, y, z ! position (x, y, z) at which to evaluate the function
double precision, intent(in) :: d       ! instantaneous separation of the stars
double precision, intent(in) :: f       ! synchronization index of the star

! Local variables
double precision domx,domy,domz
double precision t1, t2, f21q

if (d <= 0.d0) then
  droche = bad_dble
  return
endif

t1 = x**2 + y**2 + z**2
if (t1 == 0.d0) then
  droche = bad_dble
  return
endif
t1 = t1**(-1.5d0)

t2 = (d-x)**2 + y**2 + z**2
if (t2 == 0.d0) then
  droche = bad_dble
  return
endif
t2 = q*t2**(-1.5d0)

f21q = f**2*(1.d0+q)

domx = -x*t1 + (d-x)*t2 + f21q*x - q/d**2
domy = -y*(t1 + t2 - f21q)
domz = -z*(t1 + t2)
droche = sqrt(domx**2 + domy**2 + domz**2)
return
end function droche

!-------------------------------------------------------------------------------

double precision function limbdark(mu, law, ldc)
! Limb darkening
implicit none
double precision,intent(in) :: mu     ! Cosine of the viewing angle
double precision,intent(in) :: ldc(:) ! Limb darkening coefficients/data
integer, intent(in) :: law    ! Limb darkening law - see constants.f90
! If law < -1 then ldc(:) is assumed to be a grid of specific intensity values
! tabulated on a uniform grid of abs(law) mu values, for 0 to 1. The specific
! intensity for the input value of mu is then determined  by linear 
! interpolation.

integer :: i, n
double precision :: w

if (law  < -1 ) then
  n =  abs(law)
  w = 1 + mu * (n-1)
  i = int(w)
  w = w - dble(i)
  if (i < n) then 
    limbdark =  (1.0d0-w)*ldc(i) + w*ldc(i+1)
  else
    limbdark = ldc(n)
  endif
  return
endif

select case (law)

case(ld_none) ! no limb darkening
  limbdark = 1.d0
  return

case(ld_linear) ! linear
  limbdark = 1.d0 - ldc(1)*(1.d0-mu)
  return

case(ld_quadratic) ! quadratic
  limbdark = 1.d0 - ldc(1)*(1.d0-mu) - ldc(2)*(1.d0-mu)**2
  return

case(ld_sing) ! sing 3-parameter
  limbdark = 1.d0 - ldc(1)*(1.d0-mu) - ldc(2)*(1.d0-mu**1.5d0) &
           - ldc(3)*(1.d0-mu**2)
  return

case(ld_claret) ! claret 4-parameter
  limbdark = 1.d0 - ldc(1)*(1.d0-sqrt(mu)) - ldc(2)*(1.d0-mu) &
           - ldc(3)*(1.d0-mu**1.5d0)- ldc(4)*(1.d0-mu**2)
  return

case(ld_logarithmic) ! logarithmic
  limbdark = 1.d0 - ldc(1)*(1.d0-mu) - ldc(2)*mu*log(mu)
  return

case(ld_square_root) ! square-root
  limbdark = 1.d0 - ldc(1)*(1.d0-mu) - ldc(2)*(1.d0-sqrt(mu))
  return

case(ld_exponential) ! exponential
  limbdark = 1.d0 - ldc(1)*(1.d0-mu) - ldc(2)/(1.d0-exp(mu))
  return

case(ld_power2) ! power-2 
  limbdark = 1.d0 - ldc(1)*(1.d0-mu**ldc(2))
  return

case default ! 
  limbdark = bad_dble
end select

end function limbdark

!-------------------------------------------------------------------------------

function ld_quad_match(law, ldc) result(ldc_q) 
! Set coefficients of a quadratic limb darkening law so that the intensity
! profile matches at mu= 0, 0.5, 1.0
! N.B.  these are the coefficients on the quadratic limb darkening law as used
! in eker, i.e.,  I_0[1 - u_1.mu + u_2.mu^2], so u_2 is the negative of the
! normal quadratic limb darkening coefficient.
implicit none
integer, intent(in) :: law    ! Limb darkening law - see constants.f90
double precision,intent(in) :: ldc(:)    ! Limb darkening coefficients
double precision            :: ldc_q(2) ! Quadratic limb darkening coefficients

double precision :: x0, x1

select case (law)

case(ld_none) ! no limb darkening
  ldc_q = (/ 0.0d0, 0.0d0 /)
  return

case(ld_linear) ! linear
  ldc_q = (/ ldc(1), 0.0d0 /)
  return

case(ld_quadratic) ! quadratic
  ldc_q = (/ ldc(1), -ldc(2) /)
  return

case default ! other
  x0 = limbdark(0.0d0, law, ldc)
  x1 = limbdark(0.5d0, law, ldc)
  ldc_q = (/ 3.0d0 - 4.0d0*x1 + x0, -4.0d0*x1 + 2.0d0*x0 + 2.0d0 /)
  return
end select

end function ld_quad_match
        
!-------------------------------------------------------------------------------
function gmodel_coeffs(abcd, rsep, frot, qmass,  verbose) result(gmodel) 
implicit none
double precision,intent(in) :: abcd(4)
double precision,intent(in) :: rsep  
double precision,intent(in) :: frot  
double precision,intent(in) :: qmass
integer, intent(in)         :: verbose
double precision :: gmodel(3)

! Local variables
double precision,parameter  :: rtol = 0.0d0
double precision  :: gpole,gside,gback,gfront
double precision  :: ax(3),xarr(3),garr(3)

gpole = droche(qmass, 0.0d0, 0.0d0, abcd(3), rsep, frot)
gside = droche(qmass, 0.0d0, abcd(2), 0.0d0, rsep, frot)
gback = droche(qmass, abcd(4)-abcd(1), 0.0d0, 0.0d0, rsep, frot)
gfront= droche(qmass, abcd(4)+abcd(1), 0.0d0, 0.0d0, rsep, frot)

if (verbose >= v_debug) then
  print *,' gmodel_coeff:  '
  print *,' abcd = ',abcd
  print *,' rsep = ',rsep
  print *,' frot = ',frot
  print *,' g_pole = ',gpole
  print *,' g_side = ',gside
  print *,' g_back = ',gback
  print *,' g_front= ',gfront
endif

xarr = [abcd(4)+abcd(1), abcd(4), abcd(4)-abcd(1)]
garr = [gfront/gpole-1.0d0,gside/gpole-1.0d0,gback/gpole-1.0d0]
ax = parfit(xarr,garr)
if (verbose >= v_debug) then
  print *,'xfit ',xarr
  print *,'gfit ',garr
  print *,'ax = ',ax
endif
gmodel(1:3) = ax

end function gmodel_coeffs
!-------------------------------------------------------------------------------
double precision function bright(f, g, npar, par, verbose)
!
!  Surface brightness at a point on an ellipsoidal star. 
!
!  To apply Doppler boosting set vscale to the rotation speed at the point
! (x,y,z) = (0,b,0) in km/s and set kboost to the Doppler boosting factor.
! 
!  To return the flux-weighted radial velocity set rvflag/=0 and vscale.
!
! For spherical stars only, sky-projected angle lambda between orbital and
! stellar rotation axes, can be non-zero.
!
! par(1) = Intensity at the centre of the stellar disc
! par(2) = a
! par(3) = b
! par(4) = c 
! par(5) = d 
! par(6) = inc
! par(7) = phi
! par(8) = rsep (separation between centres-of-mass)
!
! int(par(9)) = ldlaw = code for limb darkening law. 
! If ldlaw > -1 ...
! par(10) = Limb darkening coefficient 1
! par(11) = Limb darkening coefficient 2
! par(12) = Limb darkening coefficient 3
! par(13) = Limb darkening coefficient 4
!
! par(14:16) = Surface gravity model coefficients (if .not.exact_grav)
!              OR
! par(14) = q, mass ratio = mass of companion/mass of star (if exact_grav)
! par(15) = g0, local surface gravity at the pole (if exact_grav)
! par(16) = frot, asynchronous rotation factor
!
! par(17) = Gravity darkening coefficient (-v in Wood's notation)
! par(18) = Incident flux from companion
! par(19) = Heating/reflection model coefficient
! par(20) = Heating/reflection model exponent
! par(21) = Heating/reflection linear limb darkening coefficient
! par(22) = Companion radius/semi-major axis
! par(23) = lambda = sky-projected angle between orbital and rotation axes.
! par(24) = vscale = rotation speed at (x,y,z)=(0,b,0) in km/s
! par(25) = kboost = Doppler boosting factor
! par(26) = rvflag /= 0 to calculate flux-weighted radial velocity.
! par(27) = { 0 to flag .not.exact_grav
!           { potential at sub-stellar point for exact_grav case
! par(28-31) = Coordinate transformation 
! par(32)    = 0 to ignore coordinate transformation
! par(33:npar) = specific intensity values from mu=0 to mu=1 (if ldlaw < -1).
!
! Position (s,t) in projected star coordinates s given by the following
! transformation of the  input coordinates  (if npar >= 31)
!   s = par(28) + par(30)*f + par(31)*g
 !  t = par(29) - par(31)*f + par(30)*g
! otherwise
!   s = f
!   t = g
! 
! **  If you use this transformation, remember to multiply the result of your
! integration by the Jacobian = par(30)**2 + par(31)**2
!
!
implicit none
integer, intent(in)  :: npar, verbose
double precision, intent(in) :: f,g,par(npar)

! Local variables
double precision :: I_0, F_0,s,t,a,b,c,d,inc,phi,rsep,ldc(4)
double precision :: gmodel(3),vgrav,H_0,H_1,uheat,rcomp, q, g0, frot
double precision :: ap,bp,cp,dp
double precision, save :: cosi,sini,sinsqi,sinsqphi,cossqi,cossqphi
double precision, save :: sin2i,sin2phi,sinphi,cosphi
double precision, save :: phi_save = not_set_dble, inc_save = not_set_dble
double precision :: p,x,xp,y,yp,z,zp,s2,t2,st,a2,b2,c2,mu,lambda,rv,kboost
double precision :: ldfac,gfac,heat,radpt,r,rsq,cosgam,dcosgam
double precision :: omsini, omcosi, sinl, cosl, pot
double precision, parameter :: c_kms = iau_c*1.0d-3 ! Speed of light in km/s
integer,parameter :: nrpar=7
double precision :: rpar(nrpar)
double precision, parameter :: rtol=1.0d-4,mutol=1.0d-12
integer :: ldlaw
logical :: rvflag,exact_grav

I_0   = par(1)
a     = par(2) 
b     = par(3) 
c     = par(4) 
d     = par(5) 
inc   = par(6) 
phi   = par(7) 
rsep  = par(8) 
ldlaw = int(par(9))
ldc(1:4)    = par(10:13) 
pot = par(27)
exact_grav = pot.ne.0.0d0
if (exact_grav) then
  q  = par(14) 
  g0 = par(15) 
  frot = par(16) 
else
  gmodel(1:3) = par(14:16)
  g0 = not_set_dble
endif
vgrav  = par(17)
F_0    = par(18) 
H_0    = par(19)
H_1    = par(20)
uheat  = par(21)
rcomp  = par(22)
lambda = par(23)
omsini = par(24)/b
kboost = par(25)
rvflag = par(26).ne.0.0d0
if (par(32) /= 0) then
  s = par(28) + par(30)*f + par(31)*g
  t = par(29) - par(31)*f + par(30)*g
else
  s = f
  t = g
endif

if ((omsini/=0.0d0).and.(lambda/=0.0d0).and.(abs(a-c)>epsilon(0.0d0))) then
  if(verbose >= v_error) then
    print *,'Error in bright - lambda /= 0 for non-spherical star.'
    print *,lambda,omsini
    print *,a,b,c,inc,phi,s,t
  endif
  bright = bad_dble
  return
endif

! Dark star, no heating ...
heat = F_0*H_0
if ((I_0 == 0.0d0).and.(heat == 0.0d0)) then
  bright = 0
  return
endif

if (inc /= inc_save) then
  inc_save = inc
  cosi = cos(inc)
  sini = sin(inc)
  cossqi = cosi**2
  sinsqi = sini**2
  sin2i = 2.d0*cosi*sini
endif
if (phi /= phi_save) then
  phi_save = phi
  cosphi = cos(phi)
  sinphi = sin(phi)
  cossqphi = cosphi**2
  sinsqphi = sinphi**2
  sin2phi = 2.d0*cosphi*sinphi
endif

a2 = a**2
b2 = b**2
c2 = c**2
s2 = s**2
t2 = t**2
st = s*t

ap = sinsqi*cossqphi/a2 + sinsqi*sinsqphi/b2 + cossqi/c2
bp = (sini*sin2phi*s - sin2i*cossqphi*t)/a2   &
   - (sini*sin2phi*s + sin2i*sinsqphi*t)/b2   &
   + sin2i*t/c2
cp = (sinsqphi*s2 + cossqi*cossqphi*t2 - cosi*sin2phi*st)/a2 &
   + (cossqphi*s2 + cossqi*sinsqphi*t2 + cosi*sin2phi*st)/b2 &
   + sinsqi*t2/c2 - 1.d0
dp = bp**2 - 4.d0*ap*cp

if (dp < 0.d0) then
  ! The is usually a round-off error problem. 
  ! Assume this is the case and do the calculation for the point on the limb
  if(verbose >= v_debug) then
    print *,'bright: No such point on surface.'
    print *,hypot(s,t)-a,hypot(s,t)-b,hypot(s,t)-c
    print *,a,b,c,inc,phi,dp
    print *,s,t,npar
    print *,f,g,par(28:31)
    print *,ap,bp,cp
  endif
  dp = 0
endif

p = (-bp + sqrt(dp))/(2.d0*ap)

! Coordinates of input point (s,t) in (x', y', z') coordinates centred on
! geometric centre of ellipsoid
xp = sinphi*s - cosi*cosphi*t + sini*cosphi*p
yp = cosphi*s + cosi*sinphi*t - sini*sinphi*p 
zp =            sini*t        + cosi*p        

! Coordinates of input point (s,t) in (x,y,z) coordinates centred on
! centre-of-mass of the star
x = xp+d
y = yp
z = zp

! cosine of the angle between the line of sight and the surface normal to a
! triaxial ellipsoid.
mu = (xp*cosphi*sini/a2 -yp*sinphi*sini/b2 + zp*cosi/c2) &
   /       sqrt((xp/a2)**2 + (yp/b2)**2 +(zp/c2)**2)   

if(verbose >= v_debug) then
  print *,'bright: mu = ',mu
endif

if (abs(mu) < mutol) then
    mu = mutol
else if (abs(mu-1.d0) < mutol) then
    mu = 1-mutol
else if ((mu < 0.d0).or.(mu > 1.d0)) then
    print *,'Error in bright - mu out of range'
    print *,a,b,c,inc,phi,dp
    print *,s,t,npar
    print *,xp,y,z,a2,b2,c2
    print *,sini,cosi,sinphi,cosphi,mu
    bright = bad_dble
    return
endif

! Limb darkening
if (ldlaw < -1) then
  ldfac = limbdark(mu, ldlaw, par(33:npar))
else
  ldfac = limbdark(mu, ldlaw, ldc)
endif
if (ldfac == bad_dble) then
  if(verbose >= v_error) then
    print *,'Error in bright - invalid parameters for limbdark.'
    print *,ldlaw,ldc
    print *,par(33),par(npar)
  endif
  bright = bad_dble
  return
endif

if (vgrav /= 0.0d0) then
  if (exact_grav) then
    rpar(1) = pot
    rpar(2) = x
    rpar(3) = y
    rpar(4) = z
    rpar(5) = rsep
    rpar(6) = q
    rpar(7) = frot
    r = zbrent(func_r,0.9d0,1.1d0,rtol,nrpar,rpar,verbose-1)
    gfac = (droche(q, r*x, r*y, r*z, rsep, frot)/g0)**vgrav
  else
    gfac =  (poly(gmodel(1:3),x)*(1.0d0-z**2/c2) + 1.0d0)**vgrav
  endif
else
  gfac = 1
endif

! Heating/reflection
if (heat /=  0.d0) then
  radpt   = sqrt(xp**2 + y**2 + z**2)
  rsq =  radpt**2 - 2.d0*xp*(rsep-d) + (rsep-d)**2
  r =  sqrt(rsq)
  cosgam = (xp*(rsep-d)/radpt - radpt)/r
  dcosgam = (radpt-rcomp)/r
  if (cosgam < -dcosgam) then
    heat = 0
  else if (cosgam > dcosgam) then
    heat = heat*(1.0d0/rsq)**H_1*(1.d0-uheat*(1.d0-mu)) 
  else 
    heat = heat*(0.5d0*(cosgam+dcosgam)/dcosgam/rsq)**H_1*(1.d0-uheat*(1.d0-mu))
  endif
endif

bright = I_0*ldfac*gfac+heat
if (omsini /= 0.0d0) then
  if (lambda == 0.0d0) then
    rv = omsini*s
  else
    cosl = cos(lambda)
    sinl = sin(lambda)
    if (sini /= 0.0d0) then
      omcosi = omsini/sini*cosi
    else
      omcosi = bad_dble
    endif
    rv = omsini*(cosl*s - sinl*t) 
  endif

  if (rvflag) then
    bright = rv*bright*(1.0d0 - kboost*rv/c_kms)
  else
    bright = bright*(1.0d0 - kboost*rv/c_kms)
  endif
else
  rv = 0
endif

if(verbose >= v_debug) then
  print '(a,11G13.4)','bright,s,t,p,xp,y,z,rv,inc,phi,dp = ',bright,s,t,p &
   ,xp,y,z,rv,inc ,dp
   ! phi,(xp/a)**2+(y/b)**2+(z/c)**2 - 1.0d0
endif

return

end function bright

!-------------------------------------------------------------------------------
double precision function eanom(m, e) 
!  Calculate the eccentric anomaly of a Keplerian orbit with
! eccentricity e, from the mean anomaly, m.
!
!  Solves Kepler''s equation using Newton-Raphson iteration using formula for
! initial estimate from Heintz DW, 'Double stars' (Reidel, 1978).
! 
!  Input:
!   m - Mean anomaly in radians.
!   e - Eccentrcity, 0 <= E < 1.
! 
!  Output:
!    Eccentric anomaly in the range 0 < eanom < 2*pi.
!    If e is out-of-range return bad_dble
!
implicit none
double precision,intent(in) :: m, e
integer, parameter  :: itmax = 9999
double precision, parameter :: etol=1.0d-9

! Local variables
double precision :: e0, e1, test
integer  :: it

if ( (e < 0.0d0).or.(e >= 1.0d0)) then
 print *, 'Invalid eccentricity value in function eanom'
 print *, e
 eanom = bad_dble
 return
endif

it = 0
e1 = mod(m,twopi) + e*sin(m) + e*e*sin(2.0d0*m)/2.0d0
test = 1.0d0
do while (test > etol)
 it = it + 1
 if (it .gt. itmax) then
  print *,'function eanom failed to converge'
  print *,m,e,e0,e1,test
  stop
 endif
 e0 = e1
 e1 = e0 + (m-(e0 - e*sin(e0)))/(1.0d0 - e*cos(e0))
 test = abs(e1 - e0)
end do
e1 = mod(e1, twopi)
if (e1 < 0.0d0) then
 e1 = e1 + twopi
endif 
eanom = e1
return
end function eanom

!-----------------------------------------------------------------------------

double precision function trueanom(m, e) 
!  Calculate the true anomaly of a Keplerian orbit with eccentricity e,
! from the mean anomaly, m.
!
! Uses: eanom
!
!  Input:
!   m - Mean anomaly in radians.
!   e - Eccentrcity, 0 <= e < 1.
!
! Output:
!   True anomaly in the range 0 to 2*PI.
!   If e is out-of-range return bad_dble
!  
implicit none
double precision,intent(in) :: m, e

! Local variables
double precision :: ee

if ( (e < 0.0d0).or.(e >= 1.0d0)) then
 print *, 'invalid eccentricity value in function trueanom'
 print *, e
 trueanom = bad_dble
 return
endif

ee = eanom(m,e)
trueanom = 2.0d0*atan(sqrt((1.0d0 + e)/(1.0d0 - e))*tan(ee/2.0d0))
return
end function trueanom

!-----------------------------------------------------------------------------

double precision function radvel(t, t0, p, v0, k, e, omrad) 
! Calculate radial velocity for a Keplerian orbit
!
! Uses: TRUEANOM
!
! Input:
!  T      - Time of observation
!  T0     - Time of periastron
!  P      - Orbital period
!  V0     - Systemic velocity
!  K      - Velocity semi-amplitude in the same units as V0
!  E      - Eccentricity of the orbit
!  OMRAD  - Longitude of periastron in radians.
!
!  Output:
!   Radial velocity in the same units as V0.
!
implicit none
double precision, intent(in) :: t, t0, p, v0, k, e, omrad

! Local variables
double precision :: m

m = twopi*mod((t-t0)/p,1.0d0)
if (e == 0.0d0) then
  radvel = v0 +k*cos(m+omrad)
  return
endif
if ( (e < 0.0d0).or.(e >= 1.0d0)) then
  print *, 'Invalid eccentricity value in function radvel'
  print *, e
  radvel = bad_dble
  return
endif

radvel = v0 + k*( e*cos(omrad) + cos(trueanom(m,e) + omrad))         
return
end function radvel

!------------------------------------------------------------------------------
double precision function t_ecl_to_peri(t_ecl, ecc, omega, incl, p_sid, verbose)
!
! Calculate the time of periastron passage immediately prior to a give time of
! eclipse. Equation numbers from Hilditch, "An Introduction to Close Binary
! Stars"
!
double precision, intent (in) :: t_ecl ! Time of eclipse
double precision, intent (in) :: ecc   ! Orbital eccentricity
double precision, intent (in) :: omega ! Longitude of periastron, radians
double precision, intent (in) :: incl  ! Orbital inclination, radians
double precision, intent (in) :: p_sid ! Siderial period
integer, intent (in) :: verbose
! Local variables
double precision :: theta, theta_0, delta_t, sin2i, efac, ee, eta, d
double precision, parameter :: tol = 1.0d-8
integer :: verbose1
integer, parameter :: npar = 4
double precision :: par(npar)

verbose1 = verbose_for_calls(verbose)
efac  = (1.0d0 - ecc**2)
sin2i = sin(incl)**2
! Value of theta for i=90 degrees
theta_0 = halfpi - omega ! True anomaly at superior conjunction
if (verbose >= v_debug) print *,'t_ecl_to_peri: theta_0 = ',theta_0
if (incl /= halfpi) then
 par = (/ efac, sin2i, omega, ecc /)
 d =  brent(theta_0-halfpi,theta_0,theta_0+halfpi, delta_func, npar, par, tol, &
            theta, verbose1)
else
  theta = theta_0
endif
if (verbose >= v_debug) print *,'t_ecl_to_peri: theta = ',theta

! (4.10)
if (theta == pi) then
  ee = pi
else
 ee = 2.0d0 * atan(sqrt((1.0d0-ecc)/(1.0d0+ecc)) * tan(theta/2.0d0))
endif
eta = ee - ecc*sin(ee)
delta_t = eta*p_sid/twopi
t_ecl_to_peri = t_ecl  - delta_t

end function t_ecl_to_peri

!------------------------------------------------------------------------------

double precision function delta_func(theta, npar, par, verbose)
implicit none
integer, intent(in)          :: npar, verbose
double precision, intent(in) :: theta, par(npar)
! Local variables
double precision :: efac, sin2i, omega, ecc

efac  = par(1)
sin2i = par(2)
omega = par(3)
ecc   = par(4)
delta_func = efac*sqrt(1.0d0 - sin2i*sin(theta+omega)**2)/(1.0d0+ecc*cos(theta))
if (verbose >= v_debug) print *,par,delta_func
end function delta_func

end module stellar
