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

module ellipse
! Module for storage and manipulation of ellipses 
!
! Ellipse data stored in one double precision array with n_ell_par parameters
!
! Ellipses can be represented as a quadratic curve
!   ax**2 + 2bxy + cy**2 + 2dx + 2fy + g = 0
!
! Ellipses can also be described by the following 5 parameters
!  a_p = semi-major axis
!  b_p = semi-minor axis
!  x_c, y_c = coordinates of centre
!  phi  = angle measured ccw from x-axis to the major axis in radians 
!
! The parametric equation of the ellipse is
!   x = a_p*cos(phi)*cos(t) - b_p*sin(phi)*sin(t) + x_c
!   y = a_p*sin(phi)*cos(t) + b_p*cos(phi)*sin(t) + y_c
!
! The ellipse parameters and properties can be accessed as follows
!  ellipse(i_ell_qcoeff) = [a, b, c, d, f, g]
!  ellipse(i_ell_ellpar) = [a_p, b_p, x_c, y_c, phi]
!  ellipse(i_ell_a_p) = a_p
!  ellipse(i_ell_b_p) = b_p
!  ellipse(i_ell_x_c) = x_c
!  ellipse(i_ell_y_c) = y_c
!  ellipse(i_ell_phi) = phi
!  ellipse(i_ell_area) = area
!  ellipse(i_ell_cosphi) = x-axis projection factor cos(phi)
!  ellipse(i_ell_sinphi) = y-axis projection factor sin(phi)
!
! Affine transformations are stored in an array tr(2,3)
!   x_t = tr(1,3) + tr(1,1)*x + tr(1,2)*y
!   y_t = tr(2,3) + tr(2,1)*x + tr(2,2)*y
! 
! HISTORY
! Apr 2019
!  Added check for very-nearly identical ellipses
!
! Jun 2017
!  Changed root polishing in ell_ell_roots to Maehly's procedure
!
! May 2016
!  Added a warning if root polishing failed to converge in ell_ell_roots
!
! Apr 2016
! Added do while loop to root polishing in ell_ell_roots
! p.maxted@keele.ac.uk
!
! Jan 2016
! First version 
! p.maxted@keele.ac.uk
!

use constants
use utils
use solve_real_poly

private

! Functions to create ellipses
public ell_init_from_par     ! New ellipse from parameters
public ell_init_from_qcoeff  ! New ellipse from coefficients of quadratic curve
public ell_project_ellipsoid ! Ellipse that is the projection of an ellipsoid
! Functions to transform ellipses
public ell_affine            ! Apply an affine transformation to an ellipse
public ell_move              ! Translation ellipse
public ell_rotate            ! Rotate ellipse around the origin
public ell_shrink            ! Scaling along x- and y- axes
public ell_ell_normalize     ! Transform to unit circle + horizontal ellipse
! Functions to find a point on an ellipse
public ell_xy_to_t           ! Parametric angle of a point on an ellipse
public ell_t_to_xy           ! Point on ellipse for a given parametric angle
! Functions to calculate properties of ellipses
public ell_evaluate_qcoeff   ! Evaluate quadratic curve polynomial at a point.
public ell_point_is_inside   ! Test if point is inside ellipse
public ell_line_intersect    ! Intersection of a line with ellipse
public ell_nearfar           ! Nearest/furthest points from the origin
public ell_ell_overlap       ! Area bounded by overlap of two ellipses
public ell_ell_intersect     ! Intersection points of two ellipses (non-tangent)
public ell_ucircle_overlap   ! Ellipse/unit circle overlap area
public ell_sector_area       ! Area of an elliptical sector

!------------------------------------------------------------------------------

! Indices for ellipse parameters
integer, parameter, public :: n_ell_par = 14
integer, parameter, public :: i_ell_qcoeff(6) = [1, 2, 3, 4, 5, 6]
integer, parameter, public :: i_ell_ellpar(5) = [7, 8, 9, 10, 11]
integer, parameter, public :: i_ell_a_p = 7
integer, parameter, public :: i_ell_b_p = 8
integer, parameter, public :: i_ell_x_c = 9
integer, parameter, public :: i_ell_y_c = 10
integer, parameter, public :: i_ell_phi = 11
integer, parameter, public :: i_ell_area = 12
integer, parameter, public :: i_ell_cosphi = 13
integer, parameter, public :: i_ell_sinphi = 14

! Bit flags for ell_ell_overlap, ell_ell_intersect, ell_ucircle_overlap 
integer, parameter, public :: b_ell_no_overlap         = 0
integer, parameter, public :: b_ell_1_inside_2         = 1
integer, parameter, public :: b_ell_2_inside_1         = 2
integer, parameter, public :: b_ell_two_intersects     = 3
integer, parameter, public :: b_ell_four_intersects    = 4
integer, parameter, public :: b_ell_identical          = 5
integer, parameter, public :: b_ell_invalid            = 7
integer, parameter, public :: b_ell_error              = 8
integer, parameter, public :: b_ell_warn_inaccurate    = 9
integer, parameter, public :: b_ell_ellipse_inside_circle = 1
integer, parameter, public :: b_ell_circle_inside_ellipse = 2

contains

function ell_init_from_par(ellpar) result(ellipse)
! Set ellipse data from parameters of the ellipse.
!  ellpar = [a_p, b_p, x_c, y_c, phi]
!  
!  If a_p <= 0 or b_p <= then all elements of ellipse are set to bad_dble
implicit none
double precision,intent(in) :: ellpar(5)
!f2py integer, parameter :: n_ell_par = 14
double precision            :: ellipse(n_ell_par) 

! Local variables
double precision :: a_p, b_p, x_c, y_c, phi
double precision :: cosphi,cosphi2,sinphi,sinphi2,asq,bsq
double precision :: tmp0,tmp1,anorm, qcoeff(6)

a_p = ellpar(1)
b_p = ellpar(2)
x_c = ellpar(3)
y_c = ellpar(4)
phi = ellpar(5)
if ((a_p <= 0.0d0).or.(b_p<=0.0d0)) then
  ellipse(:) = bad_dble
  return
endif
if (a_p < b_p) then
  call swap(a_p,b_p)
  phi = mod(phi+halfpi,twopi)
endif

! Calculate coefficients of the quadratic curve
cosphi = cos(phi)
sinphi = sin(phi)
cosphi2 = cosphi**2
sinphi2 = 1.0d0 - cosphi2
asq = a_p**2
bsq = b_p**2
anorm = min(asq,bsq)
tmp0 = (cosphi*x_c + sinphi*y_c)/asq
tmp1 = (sinphi*x_c - cosphi*y_c)/bsq
qcoeff(1) = (cosphi2/asq + sinphi2/bsq) * anorm
qcoeff(2) = cosphi*sinphi*(1.d0/asq - 1.d0/bsq) * anorm
qcoeff(3) = (sinphi2/asq + cosphi2/bsq) * anorm
qcoeff(4) = (-cosphi*tmp0 - sinphi*tmp1) * anorm
qcoeff(5) = (-sinphi*tmp0 + cosphi*tmp1) * anorm
qcoeff(6) = (tmp0**2*asq + tmp1**2*bsq - 1.d0) * anorm
ellipse(i_ell_qcoeff) = qcoeff

! Copy ellipse parameters to ellipse array
ellipse(i_ell_a_p) = a_p
ellipse(i_ell_b_p) = b_p
ellipse(i_ell_x_c) = x_c
ellipse(i_ell_y_c) = y_c
ellipse(i_ell_phi) = phi

! Calculate and store area and projection factors
ellipse(i_ell_area) = pi*a_p*b_p
ellipse(i_ell_cosphi) = cosphi 
ellipse(i_ell_sinphi) = sinphi 

end function ell_init_from_par

!-------------------------------------------------------------------------------

function ell_init_from_qcoeff(qcoeff) result(ellipse)
! Set ellipse data from coefficients of the quadratic curve.
!  ax**2 + 2bxy + cy**2 + 2dx + 2fy + g = 0
!  qcoeff = [a, b, c, d, f, g]
!
! If these coefficients do not describe an ellipse then all elements of
! ellipse are set to bad_dble
!
! Equations for a_p, b_p and phi adapted from 
!  mathworld.wolfram.com/Ellipse.html
! (equation numbers given in comments).
! Modified calculation so that a_p > b_p in all cases.
!
      
implicit none
double precision,intent(in) :: qcoeff(6)
!f2py integer, parameter :: n_ell_par = 14
double precision            :: ellipse(n_ell_par) 

! Local variables
double precision :: a, b, c, d, f, g
double precision :: a_p, b_p, x_c, y_c, phi
double precision :: num, ss, t1,t2
double precision :: b2, d2, f2, dd, jj, ii

! Initialise to error flag value so we can return on error
ellipse(:) = bad_dble

a = qcoeff(1)
b = qcoeff(2)
c = qcoeff(3)
d = qcoeff(4)
f = qcoeff(5)
g = qcoeff(6)

! Calculate ellipse parameters
b2 = b**2
d2 = d**2
f2 = f**2

! (16)
dd = a*(c*g-f2) - b*(b*g-f*d) + d*(b*f-c*d) 
if (dd == 0.d0) return
! (17)
jj = a*c - b2  
if (jj <= 0.d0) return
! (18)
ii = a + c     
if (ii == 0.d0) return
if (dd/ii >= 0.d0) return

! (19), (20)
x_c = (b*f - c*d)/jj 
y_c = (b*d - a*f)/jj

num =2.d0*(a*f2 + c*d2 + g*b2 - 2.d0*b*d*f - a*c*g) 
ss = sqrt((a-c)**2 + 4.d0*b2)

! (21), (22) - added logic here so that a_p > b_p
t1 = sqrt(num/(jj*(ii - ss)))  
t2 = sqrt(num/(jj*(ii + ss)))  
if (t1 > t2) then
  a_p = t1 
  b_p = t2
else
  a_p = t2 
  b_p = t1
endif

! (23)
if (a == c) then
  if (t1 > t2) then
    phi = 1.5d0*halfpi
  else
    phi = 0.5d0*halfpi
  endif
else if ((a < c).eqv.(t1 > t2)) then
  if (b == 0.d0) then
    phi = 0.d0
  else
    phi =  acot((a-c)/(2.d0*b))/2.d0
  endif
else
  if (b == 0.d0) then
    phi = halfpi
  else
    phi = halfpi+ acot((a-c)/(2.d0*b))/2.d0
  endif
endif

! Copy quadratic curve function into ellipse array
ellipse(i_ell_qcoeff) = qcoeff

! Copy  ellipse parameters into into ellipse array
ellipse(i_ell_a_p) = a_p
ellipse(i_ell_b_p) = b_p
ellipse(i_ell_x_c) = x_c
ellipse(i_ell_y_c) = y_c
ellipse(i_ell_phi) = phi

! Calculate and store area and projection factors
ellipse(i_ell_area) = pi*a_p*b_p
ellipse(i_ell_cosphi) = cos(phi)
ellipse(i_ell_sinphi) = sin(phi)

end function ell_init_from_qcoeff

!------------------------------------------------------------------------------

double precision function ell_evaluate_qcoeff(xy,ellipse)
! Evaluate  ax**2 + 2bxy + cy**2 + 2dx + 2fy + g 
implicit none
double precision, intent(in) :: xy(2)       ! (x, y)
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: ellipse(n_ell_par) ! Ellipse

! Local variable
double precision :: x, y
x = xy(1)
y = xy(2)
ell_evaluate_qcoeff = ellipse(1)*x**2 + 2.d0*ellipse(2)*x*y  &
                    + ellipse(3)*y**2 + 2.d0*ellipse(4)*x  &
                    + 2.d0*ellipse(5)*y + ellipse(6)
return
end function ell_evaluate_qcoeff

!-----------------------------------------------------------------------------

function ell_move(u, v, ellipse) result(ellipse_tr)
! Apply a translation (u, v) to ellipse
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: u, v, ellipse(n_ell_par)
double precision :: ellipse_tr(n_ell_par)

!Local variables
double precision :: a, b, c, d, f, g

ellipse_tr(:) = ellipse(:)

ellipse_tr(i_ell_x_c) = ellipse(i_ell_x_c) + u
ellipse_tr(i_ell_y_c) = ellipse(i_ell_y_c) + v

a = ellipse(1)
b = ellipse(2)
c = ellipse(3)
d = ellipse(4)
f = ellipse(5)
g = ellipse(6)
ellipse_tr(4) = d - a*u - b*v 
ellipse_tr(5) = f - b*u - c*v
ellipse_tr(6) = g + a*u*u  + c*v*v + 2.d0*(b*u*v - d*u - f*v)
      
end function ell_move

!------------------------------------------------------------------------------

function ell_rotate(theta, ellipse) result(ellipse_r)
! Rotate an ellipse theta radians counter-clockwise around the origin
implicit none
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: theta, ellipse(n_ell_par)
double precision :: ellipse_r(n_ell_par)

!Local variables
double precision :: ellpar(5)

ellpar(1) = ellipse(i_ell_a_p)
ellpar(2) = ellipse(i_ell_b_p)
ellpar(3) = cos(theta)*ellipse(i_ell_x_c) - sin(theta)*ellipse(i_ell_y_c)
ellpar(4) = sin(theta)*ellipse(i_ell_x_c) + cos(theta)*ellipse(i_ell_y_c)
ellpar(5) = mod(ellipse(i_ell_phi)+twopi+theta,twopi)

ellipse_r = ell_init_from_par(ellpar)

end function ell_rotate

!------------------------------------------------------------------------------

function ell_shrink(xscale, yscale, ellipse) result(ellipse_s)
! Apply change of scale x' = xscale*x and y' = yscale*y to ellipse
implicit none
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: xscale, yscale, ellipse(n_ell_par)
double precision :: ellipse_s(n_ell_par)

!Local variables
double precision :: qcoeff(6)

qcoeff(1) = xscale**2 * ellipse(1)
qcoeff(2) = xscale*yscale * ellipse(2)
qcoeff(3) = yscale**2 * ellipse(3)
qcoeff(4) = xscale * ellipse(4)
qcoeff(5) = yscale * ellipse(5)
qcoeff(6) = ellipse(6)

ellipse_s = ell_init_from_qcoeff(qcoeff)

end function ell_shrink

!-----------------------------------------------------------------------------

logical function ell_point_is_inside(xy,ellipse,q_value)
! Return true if the given point is inside the given ellipse
! Points on the circumference of ellipse return false.
implicit none
double precision, intent(in) :: xy(2)       ! (x, y)
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: ellipse(n_ell_par) ! 
double precision, intent(out), optional :: q_value 

! Local variables
double precision :: q
q = ell_evaluate_qcoeff(xy,ellipse)
ell_point_is_inside = (q.ne.0.d0) .and. &
                      ((q < 0.d0).eqv.(ellipse(i_ell_qcoeff(1)) > 0.d0))
if (present(q_value)) q_value = q

end function ell_point_is_inside

!-----------------------------------------------------------------------------

function ell_line_intersect(ellipse, qline) result(t)
!
! Intersection of a line with ellipse.
!
! The line is defined by the parametric equation
!  (x, y) = (x_0, y_0) + t*(dx, dy) 
!
! qline = [x_0, y_0, dx, dy]
!
! Values of t for which the line intersects the ellipse are returned in t(2)
!
! If there is only one intersection then t(1) = t(2)
! 
! If there are no intersections then t(1) = -huge(0.d0) and t(2) = huge(0.d0)
!
implicit none
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: ellipse(n_ell_par) ! Ellipse
double precision, intent(in) :: qline(4)    ! Parametric equation of the line
double precision :: t(2)                    ! Line parameter values at intercept

! Local variables
double precision :: a,b,c,d,f,g,x0,y0,dx,dy
double precision :: aa,bb,cc,dd

a = ellipse(1)
b = ellipse(2)
c = ellipse(3)
d = ellipse(4)
f = ellipse(5)
g = ellipse(6)
x0 = qline(1)
y0 = qline(2)
dx = qline(3)
dy = qline(4)
aa = a*dx**2 + 2.d0*b*dx*dy + c*dy**2
bb = 2.d0*(a*x0*dx + b*(x0*dy+y0*dx) + c*y0*dy + d*dx + f*dy)
cc = a*x0**2 + c*y0**2 + 2.d0*(d*x0 + b*x0*y0 +f*y0) + g
dd = bb**2 - 4.d0*aa*cc
if (dd < 0.d0) then
  t(1) = -huge(0.d0)
  t(2) =  huge(0.d0)
else
  if (bb >= 0.d0) then
    t(1) = (-bb - sqrt(dd))/(2.d0*aa)
    t(2) = 2.d0*cc/(-bb-sqrt(dd))
  else
    t(1) = 2.d0*cc/(-bb+sqrt(dd))
    t(2) = (-bb + sqrt(dd))/(2.d0*aa)
  endif
endif
return
end function ell_line_intersect

!-----------------------------------------------------------------------------


double precision function ell_xy_to_t(xy, ellipse)
! Parametric angle in the range 0 to 2.pi for a point on an ellipse.
! In fact this works for any point in space apart from the origin of the ellipse
! in which case it is the parametric angle of the point of intersection of the
! line through the centre and the input point with the ellipse. 
implicit none
double precision, intent(in) :: xy(2)  ! Point on ellipse (xy(1), xy(2))
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: ellipse(n_ell_par) ! Ellipse

! Local variables
double precision :: a_p, b_p, x_c, y_c, phi

a_p = ellipse( 7)
b_p = ellipse( 8)
x_c = ellipse( 9)
y_c = ellipse(10)
phi = ellipse(11)
ell_xy_to_t = atan2( ((xy(2)-y_c)*cos(phi)-(xy(1)-x_c)*sin(phi))/b_p, &
                     ((xy(1)-x_c)*cos(phi)+(xy(2)-y_c)*sin(phi))/a_p )
if (ell_xy_to_t < 0.0d0) ell_xy_to_t = ell_xy_to_t + twopi
end function ell_xy_to_t
                                                                        
!-------------------------------------------------------------------------------

function ell_t_to_xy(t, ellipse) result(xy)
! Point on an ellipse corresponding to a given parametric angle.
implicit none
double precision, intent(in) :: t                  ! Parametric angle
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: ellipse(n_ell_par) ! Ellipse
double precision             :: xy(2)              ! (x, y)

! Local variables
double precision :: a_p, b_p, x_c, y_c, cosphi, sinphi, cost, sint

a_p = ellipse(i_ell_a_p)
b_p = ellipse(i_ell_b_p)
x_c = ellipse(i_ell_x_c)
y_c = ellipse(i_ell_y_c)
cosphi = ellipse(i_ell_cosphi)
sinphi = ellipse(i_ell_sinphi)
cost = cos(t)
sint = sin(t)
xy(1) = a_p*cosphi*cost - b_p*sinphi*sint + x_c
xy(2) = a_p*sinphi*cost + b_p*cosphi*sint + y_c

end function ell_t_to_xy

!-------------------------------------------------------------------------------

function ell_project_ellipsoid(abc, phi, incl) result(ellipse)
! 
! Ellipse that is the projection in the sky plane (s,t) of the triaxial 
! ellipsoid with semi-axes a, b, c and centred on the origin after the
! following rotation matrix has been applied.
!
! [s]   [  sin(phi)            cos(phi)           0         ][x']
! [t] = [ -cos(phi)cos(incl)   sin(phi)cos(incl)  sin(incl) ][y']
! [p]   [  cos(phi)sin(incl)  -sin(phi)sin(incl)  cos(incl) ][z']
!
! Input:
!  abc(3) - semi-axes on x'-axis, y'-axis and z'-axis, resp.
!  phi    - angle between the line-of-sight (p-axis, towards observer)
!           projected in the x'-y' plane and the  x'-axis
!  incl   - angle between the p-axis and the z'-axis
!
! Output:
!  Projected ellipse.
!  All elements set to bad_dble for invalid input.
!
implicit none
double precision, intent(in)  :: abc(3), phi, incl
!f2py integer, parameter :: n_ell_par = 14
double precision              :: ellipse(n_ell_par)

! Local variables
double precision :: t1, t2, t3, t4, t5, t6, t7, a2, b2, c2
double precision :: sini, cosi, sinphi, cosphi
double precision :: sin2phi, cos2phi, sin2i, cos2i
double precision :: qcoeff(6), aa, cc, anorm

if ( minval(abc) <= 0.d0) then

  ellipse(:) = bad_dble

else

  sini = sin(incl)
  cosi = cos(incl)
  sinphi = sin(phi)
  cosphi = cos(phi)
  sin2phi = sinphi**2
  cos2phi = cosphi**2
  sin2i = sini**2
  cos2i = cosi**2
  a2 = 1.d0/abc(1)**2
  b2 = 1.d0/abc(2)**2
  c2 = 1.d0/abc(3)**2
  t1 = sinphi*cosphi*(a2-b2)
  t2 = sini*cosi*(c2 - cos2phi*a2 - sin2phi*b2)
  t3 = sin2i*(cos2phi*a2 + sin2phi*b2) + cos2i*c2
  t4 = sin2phi*a2 + cos2phi*b2
  t5 = cos2phi*cos2i*a2 + sin2phi*cos2i*b2 + sin2i*c2 
  t6 = t1*sini
  t7 = -t1*cosi

  aa = t6**2/t3 - t4
  cc = t2**2/t3 - t5
  anorm = max(abs(aa),abs(cc))
  qcoeff(1) = aa/anorm
  qcoeff(2) = (t2*t6/t3 - t7)/anorm
  qcoeff(3) = cc/anorm
  qcoeff(4) = 0.d0
  qcoeff(5) = 0.d0
  qcoeff(6) = 1.d0/anorm

  ellipse = ell_init_from_qcoeff(qcoeff)
endif

return 

end function ell_project_ellipsoid

!-------------------------------------------------------------------------------

function ell_affine(tr, ellipse, inverse) result(ellipse_tr)
! Apply an affine transformation to an ellipse
!  The result of applything the affine transformation to a point (x, y) is
!   x = tr(1,3) + tr(1,1)*x + tr(1,2)*y
!   y = tr(2,3) + tr(2,1)*x + tr(2,2)*y
!
! If the array tr does not contain an invertable transformation then the
! all values in the return value are set to bad_dble
!
! Set inverse=.true, to apply the inverse transformation
!
implicit none
double precision, intent(in) :: tr(2,3)
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: ellipse(n_ell_par)
double precision             :: ellipse_tr(n_ell_par)
logical, optional,intent(in) :: inverse 
logical  :: inverse_l

!Local variables
double precision :: a, b, c, d, f, g
double precision :: p, q, r, s, u, v
double precision :: at, bt, ct, dt, ft, gt, anorm
double precision :: qcoeff(6), det

if (present(inverse)) then
  inverse_l = inverse
else
  inverse_l = .false.
endif

a = ellipse(1)
b = ellipse(2)
c = ellipse(3)
d = ellipse(4)
f = ellipse(5)
g = ellipse(6)

if (inverse_l) then

  p = tr(1,1)
  q = tr(1,2)
  r = tr(2,1)
  s = tr(2,2)
  u = tr(1,3)
  v = tr(2,3)

  ! First the translation
  gt = g + a*u**2 + c*v**2 + 2.d0*(b*u*v + d*u + f*v)
  anorm = max(abs(a),abs(c))
  at  = a/anorm
  bt  = b/anorm
  ct  = c/anorm
  dt  = (d + a*u + b*v)/anorm
  ft  = (f + b*u + c*v)/anorm
  gt  = gt/anorm
  
  ! Then the rotation/scaling/shearing matrix
  qcoeff(1) = at*p**2 + 2.d0*bt*p*r + ct*r**2
  qcoeff(2) = at*p*q + bt*(s*p + r*q) + ct*r*s
  qcoeff(3) = at*q**2 + 2.d0*bt*s*q + ct*s**2
  qcoeff(4) = dt*p + ft*r
  qcoeff(5) = dt*q + ft*s
  qcoeff(6) = gt

else

  det = tr(1,1)*tr(2,2) - tr(1,2)*tr(2,1)
  if (det == 0.d0) then
    ellipse_tr(:) = bad_dble
    return
  endif
  p =  tr(2,2)/det
  q = -tr(1,2)/det
  r = -tr(2,1)/det
  s =  tr(1,1)/det
  u = -tr(1,3)
  v = -tr(2,3)

  ! First the rotation/scaling/shearing matrix
  at = a*p**2 + 2.d0*b*p*r + c*r**2
  bt = a*p*q + b*(s*p + r*q) + c*r*s
  ct = a*q**2 + 2.d0*b*s*q + c*s**2
  dt = d*p + f*r
  ft = d*q + f*s

  ! Now apply the translation
  anorm = max(abs(at),abs(ct))
  qcoeff(1) = at/anorm
  qcoeff(2) = bt/anorm
  qcoeff(3) = ct/anorm
  qcoeff(4) = (dt + at*u + bt*v)/anorm
  qcoeff(5) = (ft + bt*u + ct*v)/anorm
  qcoeff(6) = (g + at*u**2 + ct*v**2 + 2.d0*(bt*u*v + dt*u + ft*v))/anorm

endif

ellipse_tr = ell_init_from_qcoeff(qcoeff)

end function ell_affine

!------------------------------------------------------------------------------

function ell_ell_normalize(ellipse_1, ellipse_2) result (tr)
! Calculate the affine transformation that transforms ellipse_1 to a unit circle
! centred at the origin and ellipse_2 to a horizontal ellipse  centred in the
! upper quadrant.
!
! The affine transformation is returned as an array tr(2,3) that can be applied 
! to a point (x,y) as follows
!  u = tr(1,3) + tr(1,1)*x + tr(1,2)*y
!  v = tr(2,3) + tr(2,1)*x + tr(2,2)*y
! 
! In this coordinate system the parametric equation of the ellipse is
!   u = u_c + a_u cos(s)
!   v = v_c + b_u sin(s)
!  with
!   u_c >= 0
!   v_c >= 0
!   a_u >= b_u
!
!
implicit none
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: ellipse_1(n_ell_par)
double precision, intent(in) :: ellipse_2(n_ell_par)
double precision             :: tr(2,3)

! Local variables
double precision             :: s(2,3)
double precision :: uv(2), cosdphi,sindphi
double precision :: ellipse_t(n_ell_par)

! Affine transformation for rotation by -phi_1 and scaling by axis lengths.
s(1,1) =  ellipse_1(i_ell_cosphi)/ellipse_1(i_ell_a_p)
s(1,2) =  ellipse_1(i_ell_sinphi)/ellipse_1(i_ell_a_p)
s(2,1) = -ellipse_1(i_ell_sinphi)/ellipse_1(i_ell_b_p)
s(2,2) =  ellipse_1(i_ell_cosphi)/ellipse_1(i_ell_b_p)
s(1,3) = 0.d0
s(2,3) = 0.d0

! Apply this affine transformation to the origin of ellipse_1
uv(1) = ellipse_1(i_ell_x_c)
uv(2) = ellipse_1(i_ell_y_c)
uv = affine(uv,s)

! Apply the same affine transformation ellipse_2 to find its new orientation
ellipse_t = ell_affine(s, ellipse_2)
cosdphi = ellipse_t(i_ell_cosphi)
sindphi = ellipse_t(i_ell_sinphi)

tr(1,1) =  cosdphi*s(1,1) + sindphi*s(2,1) 
tr(1,2) =  cosdphi*s(1,2) + sindphi*s(2,2) 
tr(2,1) = -sindphi*s(1,1) + cosdphi*s(2,1) 
tr(2,2) = -sindphi*s(1,2) + cosdphi*s(2,2) 
tr(1,3) = -cosdphi*uv(1) - sindphi*uv(2) 
tr(2,3) =  sindphi*uv(1) - cosdphi*uv(2) 

return

end function ell_ell_normalize

!------------------------------------------------------------------------------

function ell_ell_overlap(ellipse_1, ellipse_2, verbose, &
                         xy_intersect, t_int_1, t_int_2) result (area_flags)
!
! Area in the region bounded by the overlap of two ellipses and type of overlap
! The result is returned as a 2-element array
!
! area_flags(1) = area (or bad_dble for invalid input)
!
! area_flags(2) = bit mask with the following flags
!  b_ell_no_overlap       = zero area of overlap (*)
!  b_ell_1_inside_2       = ellipse_1 inside ellipse_2 (*)
!  b_ell_2_inside_1       = ellipse_2 inside ellipse_1 (*)
!  b_ell_two_intersects   = there are two intersection points  (*)
!  b_ell_four_intersects  = there are four intersection points 
!  b_ell_identical        = ellipses are identical
!  b_ell_invalid          = invalid input
!  b_ell_error            = error during calculation
!  b_ell_warn_inaccurate  = results may be inaccurate
!
! Notes
! (*) - the ellipses may (also) share one or two tangent points in these cases.
!  To test these flags, use btest(int(area_flags(2)),<bit flag name>)
!
implicit none
double precision             :: area_flags(2)
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: ellipse_1(n_ell_par)
double precision, intent(in) :: ellipse_2(n_ell_par)
integer, intent(in)  :: verbose 
double precision, optional, intent(out)  :: xy_intersect(2,4)
double precision, optional, intent(out)  :: t_int_1(4), t_int_2(4)

! Fraction limit on areas to avoid spurious 4-point intersections due to
! near-tangent points.
double precision, parameter :: atol = 1.0d-5

! Local variables
double precision :: area_tol, separation
double precision :: xyc_1(2), xyc_2(2)
double precision :: area, xyqc(2), xy_int(2,4), qline(4), s_t_a, t_1, t_2, tq(4)
double precision :: xy1(2), xy2(2), theta_1(4), theta_2(4)
integer  :: flags, k(4), j0, j1, i, intersect_flags
logical  :: l_1, l_2, whichell

! Initialize
flags = 0
! Minimum area for integration
area_tol = atol*min(ellipse_1(i_ell_area),ellipse_2(i_ell_area))

call ell_ell_intersect(ellipse_1, ellipse_2, verbose-1, intersect_flags, &
                       xy_int, theta_1, theta_2)

if (present(xy_intersect)) xy_intersect(:,:) = xy_int(:,:)
if (present(t_int_1)) t_int_1(:) = theta_1(:)
if (present(t_int_2)) t_int_2(:) = theta_2(:)

if (btest(intersect_flags, b_ell_identical)) then
  if (verbose >= v_debug) then
    print *,'ell_ell_overlap: identical ellipses.'
  endif
  flags = ibset(flags,b_ell_identical)
  area_flags(1) = ellipse_1(i_ell_area)
  area_flags(2) = flags
  return
endif

if (btest(intersect_flags, b_ell_1_inside_2)) then
    flags = ibset(flags,b_ell_1_inside_2)
    area_flags(1) = ellipse_1(i_ell_area)
    area_flags(2) = flags
    return
endif

if (btest(intersect_flags, b_ell_2_inside_1)) then
  flags = ibset(flags,b_ell_2_inside_1)
  area_flags(1) = ellipse_2(i_ell_area)
  area_flags(2) = flags
  return
endif

if (btest(intersect_flags, b_ell_no_overlap)) then
  flags = ibset(flags, b_ell_no_overlap)
  area_flags(1) = 0.d0
  area_flags(2) = flags
  return
endif

if (btest(intersect_flags, b_ell_two_intersects)) then

  flags = ibset(flags,b_ell_two_intersects)
  xy1 = ell_t_to_xy(midangle(theta_1(1),theta_1(2)), ellipse_1)
  l_1 = ell_point_is_inside(xy1,ellipse_2)
  xyc_1(1:2) = [ellipse_1(i_ell_x_c),ellipse_1(i_ell_y_c)]
  if (l_1) then
    area =  ell_sector_area(theta_1(1),theta_1(2),ellipse_1)  &
         - triangle_area(xyc_1,xy_int(1:2,1),xy_int(1:2,2))
  else
    area = ellipse_1(i_ell_area) &
         - ell_sector_area(theta_1(1),theta_1(2),ellipse_1)  &
         + triangle_area(xyc_1,xy_int(1:2,1),xy_int(1:2,2))
  endif
  xy2 = ell_t_to_xy(midangle(theta_2(1),theta_2(2)), ellipse_2)
  l_2 = ell_point_is_inside(xy2,ellipse_1)
  xyc_2(1:2) = [ellipse_2(i_ell_x_c),ellipse_2(i_ell_y_c)]
  if (l_2) then
    area = area + ell_sector_area(theta_2(1),theta_2(2),ellipse_2)  &
         - triangle_area(xyc_2,xy_int(1:2,1),xy_int(1:2,2))
  else
    area = area + ellipse_2(i_ell_area) &
         - ell_sector_area(theta_2(1),theta_2(2),ellipse_2)  &
         + triangle_area(xyc_2,xy_int(1:2,1),xy_int(1:2,2))
  endif
  if (verbose >= v_debug) then
    print *,'ell_ell_overlap, two intersections, area = ',area
  endif
  area_flags(1) = area
  area_flags(2) = flags
  return

endif

if (btest(intersect_flags, b_ell_four_intersects)) then
  flags = ibset(flags,b_ell_four_intersects)
  
  ! Deal with intersection points in clockwise order, sorted by the angle to the
  ! x-axis measured from their centroid.
  xyqc = sum(xy_int,dim=2)/4.d0
  tq = atan2(xy_int(2,1:4)-xyqc(2),xy_int(1,1:4)-xyqc(1))
  call heapsort(tq,k)

  ! Area of quadrilateral
  area = 0.5d0*abs( &
          (xy_int(1,k(3))-xy_int(1,k(1)))*(xy_int(2,k(4))-xy_int(2,k(2))) &
        - (xy_int(1,k(4))-xy_int(1,k(2)))*(xy_int(2,k(3))-xy_int(2,k(1))) )
  ! Can catch a problem here with a near-tangent point producing four
  ! intersections all very close together
  if (area < area_tol) then
    separation = hypot(ellipse_1(i_ell_x_c)-ellipse_2(i_ell_x_c), &
                       ellipse_1(i_ell_y_c)-ellipse_2(i_ell_y_c) )
    if (separation >= (ellipse_1(i_ell_a_p)+ellipse_2(i_ell_a_p))) then
      flags = ibset(0, b_ell_no_overlap)
      area = 0
    else 
      if (ellipse_1(i_ell_area) < ellipse_2(i_ell_area)) then
        flags = ibset(0,b_ell_1_inside_2)
        area = ellipse_1(i_ell_area)
      else
        flags = ibset(0,b_ell_2_inside_1)
        area = ellipse_2(i_ell_area)
      endif
    endif
    if (verbose >= v_debug) then
      print *,'ell_ell_overlap: resolving near-tangent'
      print *,'ell_ell_overlap: ',ellipse_1
      print *,'ell_ell_overlap: ',ellipse_2
      print *,'ell_ell_overlap: ',flags
      print *,'ell_ell_overlap: ',area
    endif
    area_flags(1) = area
    area_flags(2) = flags
    return
  endif

  !  For each pair of adjacent intersection points we should either use the 
  ! area of the elliptical arc for  ellipse_1 or ellipse_2.
  !  Determine which ellipse to use for the first pair of intersection points
  ! by looking which ellipse is intersected first by a line from the centroid
  ! of the quadrilateral through the mid-point of the quadrilateral edge.
  !  Store the result of this test in whichell, then use this
  ! to  alternate the use of ellipse_1/elliptical_2 arc area as the
  ! calculation  proceeds around intersection points.

  qline(1:2) = xyqc
  qline(3) = 0.5d0*(xy_int(1,k(1))+xy_int(1,k(2))) - xyqc(1)
  qline(4) = 0.5d0*(xy_int(2,k(1))+xy_int(2,k(2))) - xyqc(2)
  t_1 = maxval(ell_line_intersect(ellipse_1, qline))
  t_2 = maxval(ell_line_intersect(ellipse_2, qline))
  whichell = (t_1 < t_2)
  xyc_1(1:2) = [ellipse_1(i_ell_x_c),ellipse_1(i_ell_y_c)]
  xyc_2(1:2) = [ellipse_2(i_ell_x_c),ellipse_2(i_ell_y_c)]
  do i=1,4
    j0 = k(i)
    j1 = k(mod(i,4)+1)
    if (whichell) then 
      !  Note cunning using of triangle area calculation with sign here to
      !  determine whether to use the larger or smaller of the two possible
      !  circle/ellipse sectors defined by the intersection points.
      s_t_a=triangle_area(xyc_1,xy_int(1:2,j0),xy_int(1:2,j1),signed=.true.)
      if(s_t_a > 0.d0) then
        area = area + ell_sector_area(theta_1(j0),theta_1(j1),ellipse_1) - s_t_a
      else
        area = area + ellipse_1(i_ell_area) &
        - ell_sector_area(theta_1(j0),theta_1(j1),ellipse_1) - s_t_a
      endif
    else 
      s_t_a=triangle_area(xyc_2,xy_int(1:2,j0),xy_int(1:2,j1),signed=.true.)
      if(s_t_a > 0.d0) then
        area = area + ell_sector_area(theta_2(j0),theta_2(j1),ellipse_2) - s_t_a
      else
        area = area + ellipse_2(i_ell_area) &
        - ell_sector_area(theta_2(j0),theta_2(j1),ellipse_2) - s_t_a
      endif

    endif
    whichell = .not.whichell
  end do

  area_flags(1) =  area
  area_flags(2) =  flags
  return
endif

if (verbose >= v_error) then
  print *,'ell_ell_overlap : unknown error'
  print *,intersect_flags
  print *,ellipse_1(i_ell_ellpar)
  print *,ellipse_2(i_ell_ellpar)
endif
area_flags(1) =  bad_dble
area_flags(2) =  b_ell_error
return

end function ell_ell_overlap

!------------------------------------------------------------------------------

subroutine ell_ell_intersect(ellipse_1, ellipse_2, verbose, & 
                             flags, xy_intersect, &
                             theta_1, theta_2)

! Non-tangent intersection points for two ellipses.
implicit none
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: ellipse_1(n_ell_par) 
double precision, intent(in) :: ellipse_2(n_ell_par)
integer, intent(in)  :: verbose
integer, intent(out) :: flags
double precision, intent(out):: xy_intersect(2,4)
! Parametric angle to intersection points for ellipse_1 and ellipse_2
double precision, intent(out), optional :: theta_1(4)
double precision, intent(out), optional :: theta_2(4)
! flags is bit mask with the following flags
!  b_ell_no_overlap       = zero area of overlap (*)
!  b_ell_1_inside_2       = ellipse_1 inside ellipse_2 (*)
!  b_ell_2_inside_1       = ellipse_2 inside ellipse_1 (*)
!  b_ell_two_intersects   = there are two intersection points  (+)
!  b_ell_four_intersects  = there are four intersection points 
!  b_ell_identical        = ellipses are identical
!  b_ell_invalid          = invalid input
!  b_ell_error            = error during calculation
!
! Notes
! (*) - the ellipses may (also) share one or two tangent points in these cases.
! (+) - the ellipses may also share an additional tangent point
!  To test these flags, use btest(int(area_flags(2)),<bit flag name>)

! Local variables
double precision :: xroots(4), yroots(4)
double precision :: xy(2,16,3), t_1(4), t_2(4), t_left, t_right
double precision :: xy_left(2), xy_right(2), xyc_1(2), xyc_2(2), separation
double precision :: ellipse_1t(n_ell_par),ellipse_2t(n_ell_par)
double precision, parameter :: roothalf = 0.70710678118654757d0
double precision :: tr1(2,3) = reshape([ roothalf,roothalf, &
                                       -roothalf,roothalf, &
                                           0.0d0,   0.0d0], [2,3])
double precision, parameter :: sin60 = 0.86602540378443860d0
double precision :: tr2(2,3) = reshape([ 0.5d0,sin60, &
                                        -sin60,0.5d0, &
                                           0.0d0,   0.0d0], [2,3])
!double precision :: tr3(2,3) = reshape([ sin60,0.5d0, &
!                                        -0.5d0,sin60, &
!                                           0.0d0,   0.0d0], [2,3])
double precision :: tr3(2,3) = reshape([ 1.0d0,0.0d0, &
                                         0.0d0,1.0d0, &
                                           0.0d0,   0.0d0], [2,3])
integer :: ix, iy, nxroots, nyroots
integer :: ii, jj, kk, ll
integer :: npair(3), ntouch, key(3)
double precision :: xy_touch(2,4)
logical :: fail, new, l_left, l_right, lx1, lx2

! Tolerance parameters

! Maximum separation between points on ellipses to be considered as 
! intersection as a fraction of the minimum
! Also used to remove duplicate intersections.
double precision, parameter :: xytol = 1.0d-9
! Tolerance on ellipse parameters to identify identical ellipses
double precision, parameter :: stol = 1d-12

flags = 0  ! Initialize/clear bit mask flags

! Deal with a few easy cases ...
separation = hypot(ellipse_1(i_ell_x_c)-ellipse_2(i_ell_x_c), &
                   ellipse_1(i_ell_y_c)-ellipse_2(i_ell_y_c) )
if (separation >= (ellipse_1(i_ell_a_p)+ellipse_2(i_ell_a_p))) then
  flags = ibset(flags, b_ell_no_overlap)
  return
endif

if (separation <= (ellipse_1(i_ell_b_p)-ellipse_2(i_ell_a_p))) then
  flags = ibset(flags, b_ell_2_inside_1)
  return
endif

if (separation <= (ellipse_2(i_ell_b_p)-ellipse_1(i_ell_a_p))) then
  flags = ibset(flags, b_ell_1_inside_2)
  return
endif

if (all(abs(ellipse_2(i_ell_ellpar)-ellipse_1(i_ell_ellpar)) < stol)) then
  flags = ibset(flags, b_ell_identical)
  return
endif


do ii = 1,3
  ! Find x-/y-values of intersection points
  select case (ii)
  case (1)
    ellipse_1t = ell_affine(tr1, ellipse_1)
    ellipse_2t = ell_affine(tr1, ellipse_2)
    call ell_ell_roots(ellipse_1t, ellipse_2t, xroots,nxroots,yroots=.false.)
    call ell_ell_roots(ellipse_1t, ellipse_2t, yroots,nyroots,yroots=.true.)
  case(2)
    ellipse_1t = ell_affine(tr2, ellipse_1)
    ellipse_2t = ell_affine(tr2, ellipse_2)
    call ell_ell_roots(ellipse_1t, ellipse_2t, xroots,nxroots,yroots=.false.)
    call ell_ell_roots(ellipse_1t, ellipse_2t, yroots,nyroots,yroots=.true.)
  case(3)
    ellipse_1t = ell_affine(tr3, ellipse_1)
    ellipse_2t = ell_affine(tr3, ellipse_2)
    call ell_ell_roots(ellipse_1t, ellipse_2t, xroots,nxroots,yroots=.false.)
    call ell_ell_roots(ellipse_1t, ellipse_2t, yroots,nyroots,yroots=.true.)
  end select

  if (verbose >= v_debug) then
    print *,'ell_ell_intersect: nxroots, nyroots = ',nxroots, nyroots
  endif


  if ((nxroots == 5).or.(nyroots == 5)) then ! Identical ellipses
    flags = ibset(flags,b_ell_identical)
    return
  endif

  if ((nxroots == -5).or.(nyroots == -5)) then ! root finding failed.
    flags = ibset(flags,b_ell_error)
  endif

  if ((nxroots < 0 ).or.(nyroots < 0 )) then ! roots may be inaccurate
    nxroots = abs(nxroots)
    nyroots = abs(nyroots)
    flags = ibset(flags,b_ell_warn_inaccurate)
  endif

  ! Store every pair-wise combination of x-,y- roots
  npair(ii) = 0 
  do ix=1,nxroots
    do iy=1,nyroots
      npair(ii) = npair(ii) + 1
      select case (ii)
      case (1)
        xy(1:2,npair(ii),1) = affine([xroots(ix), yroots(iy)],tr1,&
          inverse=.true.)
      case (2)
        xy(1:2,npair(ii),2) = affine([xroots(ix), yroots(iy)],tr2,&
          inverse=.true.)
      case (3)
        xy(1:2,npair(ii),3) = affine([xroots(ix), yroots(iy)],tr3,&
          inverse=.true.)
      end select
    end do
  end do

end do

! Search for (x,y) pairs that produce an intersection for at least 2 of the
! ellipse orientations (within distance xytol).
ntouch = 0 
do ii = 1,2
do jj = ii+1,3
  do kk = 1,npair(ii)
  do ll = 1,npair(jj)

    if (distance(xy(1:2,kk,ii),xy(1:2,ll,jj)) < xytol) then
      new = .true.
      do ix=1,ntouch
        if((distance(xy(1:2,kk,ii),xy_touch(1:2,ix)) < xytol) &
        .or.(distance(xy(1:2,ll,jj),xy_touch(1:2,ix)) < xytol))then
          new = .false.
          exit
        endif
      end do
      if (new) then
        ntouch = ntouch + 1
        if (ntouch > 4) then
          print *,'ell_ell_intersect: ntouch > 4'
          print *, ellipse_1
          print *, ellipse_2
          stop
        endif
        xy_touch(1:2,ntouch) = 0.5d0*(xy(1:2,kk,ii)+xy(1:2,ll,jj))
      endif
    endif
  end do
  end do
end do
end do
if (verbose >= v_debug) then
  print *,'ell_ell_intersect: ntouch= ',ntouch
endif


select case (ntouch)
  case (0,1)

  xyc_1(1:2) = [ellipse_1(i_ell_x_c), ellipse_1(i_ell_y_c)]
  if (ell_point_is_inside(xyc_1, ellipse_2).and. &
     (ellipse_1(i_ell_area) < ellipse_2(i_ell_area))) then
    flags = ibset(flags,b_ell_1_inside_2)
    return
  endif

  xyc_2(1:2) = [ellipse_2(i_ell_x_c),ellipse_2(i_ell_y_c)]
  if (ell_point_is_inside(xyc_2, ellipse_1).and. &
      (ellipse_2(i_ell_area) < ellipse_1(i_ell_area))) then
    flags = ibset(flags,b_ell_2_inside_1)
    return
  endif

  flags = ibset(flags, b_ell_no_overlap)
  return

case (2)

  ! Calculate parametric angle in ellipse 1 for the points where the ellipses
  ! touch
  do ii = 1,2
    t_1(ii) = ell_xy_to_t(xy_touch(1:2,ii),ellipse_1)
    t_2(ii) = ell_xy_to_t(xy_touch(1:2,ii),ellipse_2)
  end do
  
  ! Test for tangent points/intersection points by testing whether the points
  ! mid-way between the intersection points are inside/outside ellipse_2.
  t_right = 0.5d0*(t_1(1)+t_1(2))
  t_left = t_right + pi
  xy_right = ell_t_to_xy(t_left ,ellipse_1)
  xy_left  = ell_t_to_xy(t_right,ellipse_1)
  l_right = ell_point_is_inside(xy_right,ellipse_2)
  l_left  = ell_point_is_inside(xy_left ,ellipse_2)
  if (l_right.neqv.l_left) then
    xy_intersect(1:2,1:2) = xy_touch(1:2,1:2)
    if (present(theta_1)) theta_1(1:2) = t_1(1:2)
    if (present(theta_2)) theta_2(1:2) = t_2(1:2)
    flags = ibset(flags, b_ell_two_intersects)
    return
  endif

  ! Need to catch a rare case here where one of the mid-points between the
  ! intersection/tangent points is also a tangent point (or nearly so), in which
  ! case the test above fails. Double-check here whether the xy_touch(1:2,1:2)
  ! are tangent points by testing whether points just either side of these
  ! points are both inside/outside the other ellipse. 
  t_right = t_1(1)+1.0D-3
  t_left  = t_1(1)-1.0D-3
  xy_right = ell_t_to_xy(t_left ,ellipse_1)
  xy_left  = ell_t_to_xy(t_right,ellipse_1)
  l_right = ell_point_is_inside(xy_right,ellipse_2)
  l_left  = ell_point_is_inside(xy_left ,ellipse_2)
  lx1 = (l_right.neqv.l_left) ! True if xy_touch(1:2,1) is not tangent 

  t_right = t_1(2)+1.0D-3
  t_left  = t_1(2)-1.0D-3
  xy_right = ell_t_to_xy(t_left ,ellipse_1)
  xy_left  = ell_t_to_xy(t_right,ellipse_1)
  l_right = ell_point_is_inside(xy_right,ellipse_2)
  l_left  = ell_point_is_inside(xy_left ,ellipse_2)
  lx2 = (l_right.neqv.l_left) ! True if xy_touch(1:2,2) is not tangent 

  if (lx1.or.lx2) then ! At least one point is a non-tangent intersection
    ! Both points should be non-tangent, raise warning if not
    if (lx1.neqv.lx2) then
      flags = ibset(flags,b_ell_warn_inaccurate)
    endif
    xy_intersect(1:2,1:2) = xy_touch(1:2,1:2)
    if (present(theta_1)) theta_1(1:2) = t_1(1:2)
    if (present(theta_2)) theta_2(1:2) = t_2(1:2)
    flags = ibset(flags, b_ell_two_intersects)
    return
  endif


  xyc_1(1:2) = [ellipse_1(i_ell_x_c), ellipse_1(i_ell_y_c)]
  if (ell_point_is_inside(xyc_1, ellipse_2).and. &
  (ellipse_1(i_ell_area) < ellipse_2(i_ell_area))) then
    flags = ibset(flags,b_ell_1_inside_2)
    return
  endif

  xyc_2(1:2) = [ellipse_2(i_ell_x_c),ellipse_2(i_ell_y_c)]
  if (ell_point_is_inside(xyc_2, ellipse_1).and. &
  (ellipse_2(i_ell_area) < ellipse_1(i_ell_area))) then
    flags = ibset(flags,b_ell_2_inside_1)
    return
  endif

  flags = ibset(flags, b_ell_no_overlap)
  return

case (3) ! One of the points where the ellipses touch is a tangent point

  flags = ibset(flags, b_ell_two_intersects)
  do ii = 1,3
    t_1(ii) = ell_xy_to_t(xy_touch(1:2,ii),ellipse_1)
    t_2(ii) = ell_xy_to_t(xy_touch(1:2,ii),ellipse_2)
  end do
  call heapsort(t_1(1:3),key)

  t_right = 0.5d0*(t_1(key(2))+t_1(key(1)))
  t_left  = 0.5d0*(t_1(key(2))+t_1(key(3))) 
  xy_right = ell_t_to_xy(t_left ,ellipse_1)
  xy_left  = ell_t_to_xy(t_right,ellipse_1)
  l_right = ell_point_is_inside(xy_right,ellipse_2)
  l_left  = ell_point_is_inside(xy_left ,ellipse_2)
  if (l_right.eqv.l_left) then
    xy_intersect(1:2,1) = xy_touch(1:2,key(1))
    xy_intersect(1:2,2) = xy_touch(1:2,key(3))
    if (present(theta_1)) theta_1(1:2) = [t_1(key(1)),t_1(key(3))]
    if (present(theta_2)) theta_2(1:2) = [t_2(key(1)),t_2(key(3))]
    return
  endif

  t_right = 0.5d0*(t_1(key(1))+t_1(key(3)))+pi
  t_left  = 0.5d0*(t_1(key(1))+t_1(key(2))) 
  xy_right = ell_t_to_xy(t_left ,ellipse_1)
  xy_left  = ell_t_to_xy(t_right,ellipse_1)
  l_right = ell_point_is_inside(xy_right,ellipse_2)
  l_left  = ell_point_is_inside(xy_left ,ellipse_2)
  if (l_right.eqv.l_left) then
    xy_intersect(1:2,1) = xy_touch(1:2,key(2))
    xy_intersect(1:2,2) = xy_touch(1:2,key(3))
    if (present(theta_1)) theta_1(1:2) = [t_1(key(2)),t_1(key(3))]
    if (present(theta_2)) theta_2(1:2) = [t_2(key(2)),t_2(key(3))]
    return
  endif

  t_right = 0.5d0*(t_1(key(3))+t_1(key(1)))+pi
  t_left  = 0.5d0*(t_1(key(3))+t_1(key(2))) 
  xy_right = ell_t_to_xy(t_left ,ellipse_1)
  xy_left  = ell_t_to_xy(t_right,ellipse_1)
  l_right = ell_point_is_inside(xy_right,ellipse_2)
  l_left  = ell_point_is_inside(xy_left ,ellipse_2)
  if (l_right.eqv.l_left) then
    xy_intersect(1:2,1) = xy_touch(1:2,key(2))
    xy_intersect(1:2,2) = xy_touch(1:2,key(2))
    if (present(theta_1)) theta_1(1:2) = [t_1(key(1)),t_1(key(2))]
    if (present(theta_2)) theta_2(1:2) = [t_2(key(1)),t_2(key(2))]
    return
  endif

  print *,'ell_ell_interect: case(3) error - no tangent found.'
  print *,ellipse_1(i_ell_ellpar)
  print *,ellipse_2(i_ell_ellpar)
  stop

case (4)

  xy_intersect(1:2,1:4) = xy_touch(1:2,1:4)
  flags = ibset(flags, b_ell_four_intersects)
  if (present(theta_1)) then
    do ii = 1,4
    theta_1(ii) = ell_xy_to_t(xy_touch(1:2,ii),ellipse_1)
    end do
  endif
  if (present(theta_2)) then
    do ii = 1,4
      theta_2(ii) = ell_xy_to_t(xy_touch(1:2,ii),ellipse_2)
    end do
  endif
  return

case default

  print *,'ell_ell_intersect, ntouch = ', ntouch
  print *,ellipse_1
  print *,ellipse_2
  stop

end select



contains
  subroutine ell_ell_roots(ellipse_1, ellipse_2, roots, nroot, yroots)
  ! Roots of the polynomial in x or y that correspond to the intersections
  ! of two ellipses.
  ! Return nroot = -(# of roots) if root-polishing exceeds iteration limit
  ! Return nroot = -5 of algorithm fails
  double precision, intent(in) :: ellipse_1(n_ell_par) 
  double precision, intent(in) :: ellipse_2(n_ell_par)
  logical, intent(in), optional :: yroots
  double precision, intent(out) :: roots(4)
  integer, intent(out) :: nroot

  ! Test those potential real roots for which the imaginary part is < ztol. 
  double precision, parameter :: ztol = 1.0d-15

  ! Converge tolerance on x/y root position
  double precision, parameter :: xytol = 1.0d-9
  
  ! Fractional tolerance on sum of successive corrections in Maehly's 
  ! procedure to identify cases where the algorithm is bouncing around 
  ! either side of the root.
  double precision, parameter :: dtol = 1.0d-6

  double precision :: s0,s1,s2,s3,s4,s5
  double precision :: t0,t1,t2,t3,t4,t5
  double precision :: d31,d10,d30,d34,d40,d12,d32,d42,d51,d35,d45
  double precision :: e52,e20,e50,e54,e40,e21,e51,e41,e32,e53,e43
  double precision :: a(5), xreal(4), ximag(4)
  double precision :: xcand, p, p1, x1, x2, s, d1, d2
  integer :: degree, i, j, iter
  integer, parameter :: itmax = 64
  logical :: yroots_l, newroot, root_polishing_failed

  if (present(yroots)) then
    yroots_l = yroots
  else
    yroots_l = .false.
  endif
  s0 = ellipse_1(i_ell_qcoeff(6))
  s1 = 2.0D0*ellipse_1(i_ell_qcoeff(4))
  s2 = 2.0D0*ellipse_1(i_ell_qcoeff(5))
  s3 = ellipse_1(i_ell_qcoeff(1))
  s4 = 2.0D0*ellipse_1(i_ell_qcoeff(2))
  s5 = ellipse_1(i_ell_qcoeff(3))

  t0 = ellipse_2(i_ell_qcoeff(6))
  t1 = 2.0D0*ellipse_2(i_ell_qcoeff(4))
  t2 = 2.0D0*ellipse_2(i_ell_qcoeff(5))
  t3 = ellipse_2(i_ell_qcoeff(1))
  t4 = 2.0D0*ellipse_2(i_ell_qcoeff(2))
  t5 = ellipse_2(i_ell_qcoeff(3))

  if (yroots_l) then
    ! Coefficients of quartic whose roots are y values of intersections.
    d10 = s1*t0 - s0*t1
    d12 = s1*t2 - s2*t1
    d30 = s3*t0 - s0*t3
    d31 = s3*t1 - s1*t3
    d32 = s3*t2 - s2*t3
    d34 = s3*t4 - s4*t3
    d35 = s3*t5 - s5*t3
    d40 = s4*t0 - s0*t4
    d42 = s4*t2 - s2*t4
    d45 = s4*t5 - s5*t4
    d51 = s5*t1 - s1*t5
    a(1) = d34*d45-d35**2
    a(2) = d34*(d42-d51) + d31*d45 - 2.0d0*d35*d32
    a(3) = d34*(d40+d12) + d31*(d42-d51) - d32**2 - 2.0d0*d35*d30 
    a(4) = d34*d10 + d31*(d40+d12) - 2.0d0*d32*d30
    a(5) = d31*d10 - d30**2
  else
    ! Coefficients of quartic whose roots are x values of intersections.
    e20 = s2*t0 - s0*t2
    e21 = s2*t1 - s1*t2
    e50 = s5*t0 - s0*t5
    e52 = s5*t2 - s2*t5
    e51 = s5*t1 - s1*t5
    e54 = s5*t4 - s4*t5
    e53 = s5*t3 - s3*t5
    e40 = s4*t0 - s0*t4
    e41 = s4*t1 - s1*t4
    e43 = s4*t3 - s3*t4
    e32 = s3*t2 - s2*t3

    a(1) = e54*e43-e53**2
    a(2) = e54*(e41-e32) + e52*e43 - 2.0d0*e53*e51
    a(3) = e54*(e40+e21) + e52*(e41-e32) - e51**2 - 2.0d0*e53*e50 
    a(4) = e54*e20 + e52*(e40+e21) - 2.0d0*e51*e50
    a(5) = e52*e20 - e50**2
  endif
  ! Find unique real (or nearly-real) roots of these quartics
  degree = 4
  do i = 1,5
    if (a(i) == 0.0d0) then
      degree = degree -1
    else 
      exit
    endif
  end do
  if (degree == 0) then ! Identical ellipses
    nroot = 5 
    return
  endif

  ! Use Jenkins & Traub algorithm to get candidate roots
  call rpoly(a(5-degree:5),degree,xreal, ximag, fail)
  if (fail) then 
    print *,'ell_ell_roots: rpoly fail, a',fail,a
    print *,ellipse_1
    print *,ellipse_2
    print *,a
    print *,degree
    nroot = -5
    return
  endif

  nroot = 0 
  root_polishing_failed = .false.
  do i = 1, degree
    if(abs(ximag(i)) < ztol) then
      newroot = .true.
      xcand = xreal(i)
      do j = 1,nroot
        if (abs(xcand - roots(j)) <  epsilon(0.0d0)) then
          newroot = .false.
          exit
        endif
      end do
      if (newroot) then
        x1 = xcand
        x2 = xcand + 1.0d0
        iter = 0 
        d1 = not_set_dble

        do while ((abs(x1-x2) > xytol).and.(iter<itmax))
          p = x1*(x1*(x1*(x1*a(1)+a(2))+a(3))+a(4))+a(5)
          if (p == 0.0d0) then 
            x2 = x1
          else! Polish root using Maehly's method
            iter = iter + 1
            if (iter == itmax) then
              root_polishing_failed = .true.
              x2 = x1
            else
              p1 = x1*(x1*(4.0d0*x1*a(1)+3.0d0*a(2))+2.0d0*a(3))+a(4)
              s = 0
              do j = 1, nroot
                s = s + 1.0d0/(x1-roots(j))
              end do
              d2 = d1
              d1 = p/(p1 - p*s)
              if (abs(d1+d2) < dtol*abs(d1)) then
                x2 = x1
              else
                x2 = x1 - d1
              endif
              call swap(x1,x2)
            endif
          endif
        end do
        nroot = nroot + 1
        roots(nroot) = x1
      endif
    endif
  end do
  if (root_polishing_failed) nroot = -nroot

  end subroutine ell_ell_roots

!---

end subroutine ell_ell_intersect

!---------------------------------------------------------------------------

function ell_ucircle_overlap(ellipse, verbose, xy_intersect, t_ell,  &
                             t_circ) result (area_flags)
!
! Area in the region bounded by the overlap of an ellipse and a unit circle
! and type of overlap. The result is returned as a 2-element array
!
! area_flags(1) = area (or bad_dble for invalid input)
!
! area_flags(2) = bit mask with the following flags
!  b_ell_no_overlap            = zero area of overlap (*)
!  b_ell_ellipse_inside_circle = ellipse inside unit circle (*)
!  b_ell_circle_inside_ellipse = unit circle inside ellipse (*)
!  b_ell_two_intersects        = there are two intersection points  (*)
!  b_ell_four_intersects       = there are four intersection points 
!  b_ell_identical             = ellipses are identical
!  b_ell_invalid               = invalid input
!  b_ell_error                 = error during calculation
!  b_ell_warn_inaccurate       = results may be inaccurate
!
! Notes
! (*) - the ellipse and circle may (also) share one or two tangent points in 
!       these cases.
!
!  To test these flags, use btest(int(area_flags(2)),<bit flag name>)
!
implicit none
double precision             :: area_flags(2)
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: ellipse(n_ell_par)
integer, intent(in)  :: verbose 
double precision, optional, intent(out) :: xy_intersect(2,4)
double precision, optional, intent(out) :: t_ell(4), t_circ(4)

! Local variables
double precision :: xyc(2)
double precision :: area, xyqc(2), xy_int(2,4), qline(4), s_t_a, t_1, t_2, tq(4)
double precision :: xy1(2), xy2(2), theta_e(4), theta_c(4)
integer  :: flags, k(4), j0, j1, i, intersect_flags
logical  :: l_1, l_2, whichell
integer :: verbose1
double precision :: circle(n_ell_par) = [1.0d0,0.0d0,1.0d0,0.0d0,0.0d0,-1.0d0,&
                                         1.0d0,1.0d0,0.0d0,0.0d0,0.0d0, 0.0d0,&
                                         1.0d0,0.0d0]
! Initialize
flags = 0
verbose1 = verbose_for_calls(verbose)

call ell_ell_intersect(ellipse, circle, verbose1, intersect_flags, &
                           xy_int, theta_e, theta_c)

if (present(xy_intersect)) xy_intersect(:,:) = xy_int(:,:)
if (present(t_ell)) t_ell(:) = theta_e(:)
if (present(t_circ)) t_circ(:) = theta_c(:)

if (btest(intersect_flags, b_ell_identical)) then
  if (verbose >= v_debug) then
    print *,'ell_ucircle_overlap: ellipse identical to unit circle.'
  endif
  flags = ibset(flags,b_ell_identical)
  area_flags(1) = pi
  area_flags(2) = flags
  return
endif

if (btest(intersect_flags, b_ell_circle_inside_ellipse)) then
  if (verbose >= v_debug) then
    print *,'ell_ucircle_overlap: unit circle inside ellipse.'
  endif
  flags = ibset(flags,b_ell_circle_inside_ellipse)
  area_flags(1) = pi
  area_flags(2) = flags
  return
endif

if (btest(intersect_flags, b_ell_ellipse_inside_circle)) then
  if (verbose >= v_debug) then
    print *,'ell_ucircle_overlap: unit circle inside ellipse.'
  endif
  flags = ibset(flags,b_ell_ellipse_inside_circle)
  area_flags(1) = ellipse(i_ell_area)
  area_flags(2) = flags
  return
endif

if (btest(intersect_flags, b_ell_no_overlap)) then
  if (verbose >= v_debug) then
    print *,'ell_ucircle_overlap: no overlap.'
  endif
  flags = ibset(flags, b_ell_no_overlap)
  area_flags(1) = 0.d0
  area_flags(2) = flags
  return
endif

if (btest(intersect_flags, b_ell_two_intersects)) then

  flags = ibset(flags,b_ell_two_intersects)
  xy1 = ell_t_to_xy(midangle(theta_e(1),theta_e(2)), ellipse)
  l_1 = (xy1(1)**2 + xy1(2)**2) < 1.0d0 
  xyc(1:2) = [ellipse(i_ell_x_c),ellipse(i_ell_y_c)]
  if (l_1) then
    area =  ell_sector_area(theta_e(1),theta_e(2),ellipse)  &
         - triangle_area(xyc,xy_int(1:2,1),xy_int(1:2,2))
  else
    area = ellipse(i_ell_area) &
         - ell_sector_area(theta_e(1),theta_e(2),ellipse)  &
         + triangle_area(xyc,xy_int(1:2,1),xy_int(1:2,2))
  endif
  t_2 = midangle(theta_c(1),theta_c(2))
  xy2 = [ cos(t_2), sin(t_2) ]
  l_2 = ell_point_is_inside(xy2,ellipse)
  if (l_2) then
    area = area + arc_area(xy_int(1:2,1),xy_int(1:2,2))
  else
    area = area + pi - arc_area(xy_int(1:2,1),xy_int(1:2,2))
  endif
  if (verbose >= v_debug) then
    print *,'ell_ucircle_overlap, two intersections, area = ',area
  endif
  area_flags(1) = area
  area_flags(2) = flags
  return

endif

if (btest(intersect_flags, b_ell_four_intersects)) then
  flags = ibset(flags,b_ell_four_intersects)
  
  ! Deal with intersection points in clockwise order, sorted by the angle to the
  ! x-axis measured from their centroid.
  xyqc = sum(xy_int,dim=2)/4.d0
  tq = atan2(xy_int(2,1:4)-xyqc(2),xy_int(1,1:4)-xyqc(1))
  call heapsort(tq,k)

  ! Area of quadrilateral
  area = 0.5d0*abs( &
          (xy_int(1,k(3))-xy_int(1,k(1)))*(xy_int(2,k(4))-xy_int(2,k(2))) &
        - (xy_int(1,k(4))-xy_int(1,k(2)))*(xy_int(2,k(3))-xy_int(2,k(1))) )
  !  For each pair of adjacent intersection points we should either use the 
  ! area of the elliptical arc or the circular arc. Determine which to use for
  ! the first pair of intersection points by looking whether the ellipse or the
  ! circle is intersected first by a line from the centroid of the 
  ! quadrilateral through the mid-point of the quadrilateral edge.
  !  Store the result of this test in whichell, then use this
  ! to  alternate the use of ellipse/circle arc area as the
  ! calculation  proceeds around intersection points.

  qline(1:2) = xyqc
  qline(3) = 0.5d0*(xy_int(1,k(1))+xy_int(1,k(2))) - xyqc(1)
  qline(4) = 0.5d0*(xy_int(2,k(1))+xy_int(2,k(2))) - xyqc(2)
  t_1 = maxval(ell_line_intersect(ellipse, qline))
  t_2 = maxval(ucircle_line_intersect(qline))
  whichell = (t_1 < t_2)
  xyc(1:2) = [ellipse(i_ell_x_c),ellipse(i_ell_y_c)]
  do i=1,4
    j0 = k(i)
    j1 = k(mod(i,4)+1)
    if (whichell) then 
      !  Note cunning using of triangle area calculation with sign here to
      !  determine whether to use the larger or smaller of the two possible
      !  circle/ellipse sectors defined by the intersection points.
      s_t_a=triangle_area(xyc,xy_int(1:2,j0),xy_int(1:2,j1),signed=.true.)
      if(s_t_a > 0.d0) then
        area = area + ell_sector_area(theta_e(j0),theta_e(j1),ellipse) - s_t_a
      else
        area = area + ellipse(i_ell_area) &
        - ell_sector_area(theta_e(j0),theta_e(j1),ellipse) - s_t_a
      endif
    else 
      s_t_a = xy_int(1,j0)*xy_int(2,j1) - xy_int(1,j1)*xy_int(2,j0)
      if(s_t_a > 0.d0) then
        area = area + arc_area(xy_int(1:2,j0),xy_int(1:2,j1))
      else
        area = area + pi - arc_area(xy_int(1:2,j0),xy_int(1:2,j1))
      endif

    endif
    whichell = .not.whichell
  end do

  area_flags(1) =  area
  area_flags(2) =  flags
  return
endif

if (verbose >= v_error) then
  print *,'ell_ucircle_overlap : unknown error'
  print *,intersect_flags
  print *,ellipse(i_ell_ellpar)
endif
area_flags(1) =  bad_dble
area_flags(2) =  b_ell_error
return

end function ell_ucircle_overlap

!---------------------------------------------------------------------------

double precision function ell_sector_area(t_0,t_1,ellipse)
! Area of the elliptical sector defined by the parametric angles t0 and t1
! N.B. This is always the smallest possible sector defined by t0 and t1
! irrespective of whether t0 > t1 or t0 < t1.
implicit none
!f2py integer, parameter :: n_ell_par = 14
double precision,intent(in) :: ellipse(n_ell_par) 
double precision,intent(in) :: t_0, t_1

! Local variables
double precision :: dt

dt = mod(abs(t_1 - t_0), twopi)
if (dt > pi) dt = twopi - dt
ell_sector_area = dt*ellipse(i_ell_area)/twopi
return
end function ell_sector_area

!------------------------------------------------------------------------------

end module ellipse
