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

MODULE  constants
!
! Program constants for ellc and an interface routine to make these accesible
! from python
!
! HISTORY
! -------
! Jan 2016
! First version 
! p.maxted@keele.ac.uk
!
!

implicit none

public  !  All constants are here to be used. 

! Return error values for functions.
! Integer error return value
integer,  parameter :: bad_int  = -2147483647     
! Double precision error return value
double precision, parameter :: bad_dble = -1.7976931348623157d+308

! Dummy value for quantities not yet defined
! Undefined double precision
double precision, parameter :: not_set_dble = 8.9884656743115785d+307
! Undefined integer
integer,  parameter :: not_set_int  = 2147483646

! Number of parameters in the light curve model
integer, parameter :: n_par = 39   ! Number of double-precision parameters
integer, parameter :: n_ipar =10   ! Number of integer parameters

! Star shape models
integer, parameter :: starshape_roche_v   = -2  ! Roche model, exact volume
integer, parameter :: starshape_roche     = -1  ! Roche model
integer, parameter :: starshape_sphere    =  0  ! Sphere
integer, parameter :: starshape_poly1p5   =  1  ! polytrope n = 1.5
integer, parameter :: starshape_poly3p0   =  2  ! polytrope n = 3
integer, parameter :: starshape_love      =  3  ! Correia (2014) model
character*7, parameter, dimension(-2:3) :: &
 starshape_name = (/ 'roche_v', 'roche  ', 'sphere ', 'poly1p5', &
 'poly3p0', 'love   ' /)

! Indices for spot parameters in a spot data array
integer, parameter :: n_spot_par = 4
integer, parameter :: i_spot_lon  = 1 ! longitude, alpha
integer, parameter :: i_spot_lat  = 2 ! latitude, beta
integer, parameter :: i_spot_gam  = 3 ! opening angle/size, gamma
integer, parameter :: i_spot_fac  = 4 ! contrast factor, f

! Verbosity levels
! Tests against these values should always be of the form
! "IF (verbose >= v_error) ... " since verbosity levels can be outside this
! range of values. 
integer, parameter :: v_silent = 0  ! strictly no output
integer, parameter :: v_error  = 1  ! print error messages only
integer, parameter :: v_warn   = 2  ! print error messages and warnings
integer, parameter :: v_user   = 3  ! print user information
integer, parameter :: v_debug  = 4  ! debugging level

! Limb darkening laws
! N.B. values < -1 are the negative of the number of points in an array of 
! specific intensity values over a regular grid of mu values from 0 to 1.
integer, parameter :: ld_mugrid       =  -1
integer, parameter :: ld_none         =  0
integer, parameter :: ld_linear       =  1
integer, parameter :: ld_quadratic    =  2
integer, parameter :: ld_sing         =  3
integer, parameter :: ld_claret       =  4
integer, parameter :: ld_logarithmic  =  5
integer, parameter :: ld_square_root  =  6
integer, parameter :: ld_exponential  =  7
integer, parameter :: ld_power2       =  8

! Bit flags for eclipse types
integer, parameter :: bit_eclipse          =  0 ! There is an eclipse
integer, parameter :: bit_eclipse_of_star1 =  1 ! Star 1 is eclipsed by star 2
integer, parameter :: bit_eclipse_of_star2 =  2 ! Star 2 is eclipsed by star 1
integer, parameter :: bit_eclipse_is_total =  3 ! Eclipse is total
integer, parameter :: bit_eclipse_is_a_transit =  4 ! Eclipse is a transit 
! Eclipse is "double-partial" (4 intersection points).
integer, parameter :: bit_eclipse_is_double_partial  =  5

! Useful mathematical constants

double precision, parameter :: pi =  3.141592653589793238d0 ! pi
double precision, parameter :: rtod = 57.295779513082320876d0 ! rad to deg
double precision, parameter :: dtor =  0.017453292519943295d0 ! deg to rad
double precision, parameter :: twopi   = 6.283185307179586476d0 ! 2*pi
double precision, parameter :: threepi = 9.424777960769379715d0 ! 3*pi
double precision, parameter :: halfpi  = 1.570796326794896619d0 ! pi/2
double precision, parameter :: third   = 0.333333333333333333d0 ! 1/3
double precision, parameter :: fourpi_3 = 4.188790204786391d0    ! 4pi/3
  
! Astronomical and physical constants from 
! The IAU 2009 system of astronomical constants: the report of the IAU working
! group on numerical standards for Fundamental Astronomy
! Luzum  et al., ! 2011CeMDA.110..293L
! Speed of light, SI units
double precision, parameter :: iau_c  = 2.99792458d8     
! Astronomical units, SI units
double precision, parameter :: iau_au = 1.49597870700d11 
! Gravitational constant, SI units
! 2014 CODATA value from Mamajek et al. arXiv:1510.07674
double precision, parameter :: iau_G = 6.67408d-11
  
! IAU nominal values of solar constants in SI units and related constants and
! functions in consistent units.
! From www.iau.org/static/resolutions/IAU2015_English.pdf
!
! Solar radius, SI units
double precision, parameter :: solar_radius = 6.957d8
! Solar luminosity, SI units
double precision, parameter :: solar_luminosity = 3.828d26
! Solar mass * gravitational constant, SI units
double precision, parameter :: solar_gm = 1.3271244d20
! Solar mass, SI units
double precision, parameter :: solar_mass = solar_gm/iau_G
! Solar effective temperature, SI units
double precision, parameter :: solar_teff = 5772.0d0

! a.sin(i) in solar radii for K in km/s, P in mean solar days
double precision, parameter :: solar_asini_kms_d = 0.019765685498705825d0

! Mean solar day in seconds - taken from
! http://tycho.usno.navy.mil/leapsec.html, Apr 2015.
double precision, parameter :: seconds_per_mean_solar_day = 86400.002d0

contains 

integer function get_integer(variable_name)
implicit none
character(*), intent(in) :: variable_name

if (index(variable_name,'bad_int') == 1) then
  get_integer = bad_int
else if (index(variable_name,'n_par') == 1) then
  get_integer = n_par
else if (index(variable_name,'n_ipar') == 1) then
  get_integer = n_ipar
else if (index(variable_name,'starshape_roche_v') == 1) then
  get_integer = starshape_roche_v
else if (index(variable_name,'starshape_roche') == 1) then
  get_integer = starshape_roche
else if (index(variable_name,'starshape_poly1p5') == 1) then
  get_integer = starshape_poly1p5
else if (index(variable_name,'starshape_poly3p0') == 1) then
  get_integer = starshape_poly3p0
else if (index(variable_name,'n_spot_par') == 1) then
  get_integer = n_spot_par
else if (index(variable_name,'i_spot_lon') == 1) then
  get_integer = i_spot_lon
else if (index(variable_name,'i_spot_lat') == 1) then
  get_integer = i_spot_lat
else if (index(variable_name,'i_spot_gam') == 1) then
  get_integer = i_spot_gam
else if (index(variable_name,'i_spot_fac') == 1) then
  get_integer = i_spot_fac
else if (index(variable_name,'ld_none') == 1) then
  get_integer = ld_none
else if (index(variable_name,'ld_linear') == 1) then
  get_integer = ld_linear
else if (index(variable_name,'ld_quadratic') == 1) then
  get_integer = ld_quadratic
else if (index(variable_name,'ld_sing') == 1) then
  get_integer = ld_sing
else if (index(variable_name,'ld_claret') == 1) then
  get_integer = ld_claret
else if (index(variable_name,'ld_logarithmic') == 1) then
  get_integer = ld_logarithmic
else if (index(variable_name,'ld_square_root') == 1) then
  get_integer = ld_square_root
else if (index(variable_name,'ld_exponential') == 1) then
  get_integer = ld_exponential
else if (index(variable_name,'ld_power2') == 1) then
  get_integer = ld_power2
else if (index(variable_name,'ld_mugrid') == 1) then
  get_integer = ld_mugrid
else if (index(variable_name,'bit_eclipse') == 1) then
  get_integer = bit_eclipse
else if (index(variable_name,'bit_eclipse_of_star1') == 1) then
  get_integer = bit_eclipse_of_star1
else if (index(variable_name,'bit_eclipse_of_star2') == 1) then
  get_integer = bit_eclipse_of_star2
else if (index(variable_name,'bit_eclipse_is_total') == 1) then
  get_integer = bit_eclipse_is_total
else if (index(variable_name,'bit_eclipse_is_a_transit') == 1) then
  get_integer = bit_eclipse_is_a_transit
else if (index(variable_name,'bit_eclipse_is_double_partial') == 1) then
  get_integer = bit_eclipse_is_double_partial
else
  get_integer = bad_int
endif
return

end function get_integer

end module constants

