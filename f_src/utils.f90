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

module utils
!
! Utility functions and general numerical procedures for ellc
!
! HISTORY
! -------
! Jan 2016
! First version 
! p.maxted@keele.ac.uk
!
!

use constants

implicit none

private  

public swap              ! swap double precision values 
public swap_int          ! swap integer values 
public slowsort          ! sort dble values by direct insertion (slow)
public heapsort          ! calculate sort key for array of dble values 
public poly              ! evaluate a polynomial
public eval_poly_deriv   ! evaluate a polynomial and its derivative
public zbrent            ! Brent's algorithm for bracketed root of a function
public brent             ! Minimisation of a bracketed root
public rtsafe            ! Find bracketed root,  Newton's method + bisection
public laguerre          ! Find a root of a polynomial using Laguerre's method
public affine            ! Apply affine transformation to an array of points
public determinant       ! Determinant of a 2x2 or 3x3 matrix
public deflate_poly      ! Coefficients of a polynomial divided by (x-c)
public midangle          ! Bisector of two angles
public triangle_area     ! Area of a triangle defined by three points
public arc_area          ! Area of an arc on a unit circle
public ucircle_line_intersect ! Intersection of a line with unit circle
public circles_intersect  ! Intersection points of two circles on a plane
public sph_circ_intersect ! Intersection points of two circles on a sphere
public distance          ! Separation between two points in 2-d Euclidean space
public angular_distance  ! Great circle distance between points on a sphere.
public parfit            ! Parabolic fit to three points
public unitfunc          ! =1
public verbose_for_calls ! Set verbose value for function calls.
public acot              ! Inverse cotangent function
public cross3            ! Cross product of 3-element vectors

! Bit flags for circles_intersect and sph_circ_intersect
integer, parameter, public :: b_circles_no_overlap         = 0
integer, parameter, public :: b_circles_1_inside_2         = 1
integer, parameter, public :: b_circles_2_inside_1         = 2
integer, parameter, public :: b_circles_tangent            = 3
integer, parameter, public :: b_circles_intersect          = 4
integer, parameter, public :: b_circles_identical          = 5
integer, parameter, public :: b_circles_invalid            = 6

contains

!-------------------------------------------------------------------------------

subroutine swap(a,b)
  implicit none
  double precision, intent(inout) :: a, b
  !f2py intent(in,out) :: a, b
  double precision :: t
  t=b
  b=a
  a=t
  return
end subroutine swap

!-------------------------------------------------------------------------------

subroutine swap_int(a,b)
  implicit none
  integer, intent(inout) :: a, b
  integer :: t
  t=b
  b=a
  a=t
  return
end subroutine swap_int

!-------------------------------------------------------------------------------

subroutine slowsort(a,reverse)

double precision, intent(inout) :: a(:)    ! Array of values to be sorted
logical, intent(in), optional :: reverse ! Set true to sort in descending order

! Local variables
integer ::  j, k, l
logical :: reverse_l

if (present(reverse)) then
  reverse_l = reverse
else
  reverse_l = .false.
endif

if (reverse_l) then
  do j=1,size(a)-1
    k=maxloc(a(j:),1)
    l = j-1+k
    if (j /= l) call swap(a(l),a(j))
  end do
else
  do j=1,size(a)-1
    k=minloc(a(j:),1)
    l = j-1+k
    if (j /= l) call swap(a(l),a(j))
  end do
endif
return
end subroutine slowsort

!-------------------------------------------------------------------------------

subroutine heapsort(xdata, key) 
  implicit none 
  double precision, intent(in) :: xdata(:)
  integer, intent(out) :: key(:) 

! Local variables
  double precision  :: x 
  integer   :: ndata, i, l, ir, ik, j 

  ndata = size(xdata)
  if (size(key) < ndata) then
   key(1) = bad_int
   return
  endif

  do i = 1, ndata 
    key(i) = i 
  end do 
  if(ndata == 1) return 
!
  l = ndata/2+1 
  ir = ndata 
  do
    if(l > 1) then 
      l = l - 1 
      ik = key(l) 
      x = xdata(ik) 
    else 
      ik = key(ir) 
      x = xdata(ik) 
      key(ir) = key(1) 
      ir = ir - 1 
      if(ir == 1) then 
        key(1) = ik 
        return 
      endif 
    endif 
    i = l 
    j = l + l 
    do while (j <= ir) 
      if(j < ir) then 
        if(xdata(key(j)) < xdata(key(j+1))) j = j + 1 
      endif 
      if(x < xdata(key(j))) then 
        key(i) = key(j) 
        i = j 
        j = j + j 
      else 
        j = ir + 1 
      endif 
    end do
    key(i) = ik 
  end do 
end subroutine heapsort
                                                                        
!-------------------------------------------------------------------------------

double precision function poly(a,x,reverse)
implicit none
double precision, intent(in) :: a(:) ! Array of coefficients
double precision, intent(in) :: x    ! Value at which to evaluate polynomial
logical, intent(in), optional :: reverse ! Coeffs in descending order if true

! Local variables
logical :: reverse_l
integer :: i, n

if (present(reverse)) then
  reverse_l = reverse
else
  reverse_l = .false.
endif

n = size(a)

if (reverse_l) then

  select case (n)
  case(:0)
    poly=bad_dble
  case(1) 
    poly=a(1)
  case(2) 
    poly=x*a(1)+a(2)
  case(3) 
    poly=x*(x*a(1)+a(2))+a(3)
  case(4) 
    poly=x*(x*(x*a(1)+a(2))+a(3))+a(4)
  case(5) 
    poly=x*(x*(x*(x*a(1)+a(2))+a(3))+a(4))+a(5)
  case(6) 
    poly=x*(x*(x*(x*(x*a(1)+a(2))+a(3))+a(4))+a(5))+a(6)
  case(7) 
    poly=x*(x*(x*(x*(x*(x*a(1)+a(2))+a(3))+a(4))+a(5))+a(6))+a(7)
  case(8) 
    poly=x*(x*(x*(x*(x*(x*(x*a(1)+a(2))+a(3))+a(4))+a(5))+a(6))+a(7))+a(8)
  case(9) 
    poly=x*(x*(x*(x*(x*(x*(x*(x*a(1)+a(2))+a(3))+a(4))+a(5)) &
        +a(6))+a(7))+a(8))+a(9)
  case(10) 
    poly=x*(x*(x*(x*(x*(x*(x*(x*(x*a(1)+a(2))+a(3))+a(4))+a(5)) &
        +a(6))+a(7))+a(8))+a(9))+a(10)
  case(11) 
    poly=x*(x*(x*(x*(x*(x*(x*(x*(x*(x*a(1)+a(2))+a(3))+a(4)) &
        +a(5))+a(6))+a(7))+a(8))+a(9))+a(10))+a(11)
  case default
    poly=a(1)
    do i=2,n
      poly = poly*x + a(i)
    end do
  end select

  else

    select case (n)
    case(:0)
      poly=bad_dble
    case(1) 
      poly=a(1)
    case(2) 
      poly=x*a(2)+a(1)
    case(3) 
      poly=x*(x*a(3)+a(2))+a(1)
    case(4) 
      poly=x*(x*(x*a(4)+a(3))+a(2))+a(1)
    case(5) 
      poly=x*(x*(x*(x*a(5)+a(4))+a(3))+a(2))+a(1)
    case(6) 
      poly=x*(x*(x*(x*(x*a(6)+a(5))+a(4))+a(3))+a(2))+a(1)
    case(7) 
      poly=x*(x*(x*(x*(x*(x*a(7)+a(6))+a(5))+a(4))+a(3))+a(2))+a(1)
    case(8) 
      poly=x*(x*(x*(x*(x*(x*(x*a(8)+a(7))+a(6))+a(5))+a(4))+a(3)) &
          +a(2))+a(1)
    case(9) 
      poly=x*(x*(x*(x*(x*(x*(x*(x*a(9)+a(8))+a(7))+a(6))+a(5)) &
          +a(4))+a(3))+a(2))+a(1)
    case(10) 
      poly=x*(x*(x*(x*(x*(x*(x*(x*(x*a(10)+a(9))+a(8))+a(7)) &
          +a(6))+a(5))+a(4))+a(3))+a(2))+a(1)
    case(11) 
      poly=x*(x*(x*(x*(x*(x*(x*(x*(x*(x*a(11)+a(10))+a(9)) &
          +a(8))+a(7))+a(6))+a(5))+a(4))+a(3))+a(2))+a(1)
    case default
      poly=a(n)
      do i=n-1,1,-1
        poly = poly*x + a(i)
      end do
    end select

    endif 

    return
    end function poly

!-----------------------------------------------------------------------------

double precision function zbrent(func,x1,x2,tol,npar,par,verbose)
!
! Find root of function func(x,npar,par,verbose-1) using Brent's method
!
implicit none 
integer, intent(in)  :: npar, verbose
double precision, intent(in) :: x1, x2, tol, par(npar)
interface
  function func(fx,fnpar, fpar, fverbose)
  double precision             :: func
  double precision, intent(in) :: fx
  integer, intent(in)  :: fnpar
  double precision, intent(in) :: fpar(fnpar)
  integer, intent(in)  :: fverbose
  end function func
end interface
integer :: iter
double precision :: fa, fb, fc, a, b, c, d, e, tol1, xm, p, q, r, s
integer, parameter :: itmax = 100
double precision, parameter :: eps = epsilon(0.0d0)

if (verbose >= v_debug) print *,'Enter zbrent'
zbrent = bad_dble

a = x1
b = x2
c = 0.d0
d = 0.d0
e = 0.d0
fa = func(a,npar,par,verbose-1)
fb = func(b,npar,par,verbose-1)
if ((fb == bad_dble).or.(fa == bad_dble)) return

if (fa*fb > 0.0d0) then
  if (verbose >= v_error) then
    print *, 'zbrent: input values do not bracket root'
    print *, a,b
    print *, fa,fb
  endif
  return
endif

fc = fb
do iter=1, itmax
  if (fb*fc > 0.0d0) then
    c = a
    fc = fa
    d = b-a
    e=d
  endif
  if (abs(fc) < abs(fb)) then
    a = b
    b = c
    c = a
    fa = fb
    fb = fc
    fc = fa
  endif
  tol1 = 2.0d0*eps*abs(b)+0.5d0*tol
  xm = (c-b)/2.0d0
  if (abs(xm) <  tol1 .or. fb == 0.0d0) then
    zbrent = b
    return
  endif
  if (abs(e) > tol1 .and. abs(fa) >  abs(fb)) then
    s = fb/fa
    if (a == c) then
      p = 2.0d0*xm*s
      q = 1.0d0-s
    else
      q = fa/fc
      r = fb/fc
      p = s*(2.0d0*xm*q*(q-r)-(b-a)*(r-1.0d0))
      q = (q-1.0d0)*(r-1.0d0)*(s-1.0d0)
    endif
    if (p > 0.0d0) q = - q
    p = abs(p)
    if (2.0d0*p < min(3.0d0*xm*q-abs(tol1*q),abs(e*q))) then
      e = d
      d = p/q
    else
      d = xm
      e = d
    endif
  else
    d = xm
    e = d
  endif
  a = b
  fa = fb
  if( abs(d) > tol1) then
    b = b + d
  else
    b = b + sign(tol1, xm)
  endif
  fb = func(b,npar,par,verbose-1)
  if (fb == bad_dble) return
end do
if (verbose >= v_error) then
  print *, 'zbrent exceeded maximum iterations.'
  print *, a,b
  print *, fa,fb
endif
return
end function zbrent

!-----------------------------------------------------------------------------

double precision function rtsafe(func,dfunc,x1,x2,tol,npar,par,verbose)
!
! Find root of function func(x,npar,par,verbose-1) with derivative
! dfunc(x,npar,par,verbose-1) Newton-Raphson iteration combined
! with bisection.
!
! Follows the algorithm described in Numerical Recipes in Fortran 90
! (Press et al., http://www.nr.com)
!
! If the algorithm fails to converge after maxit iterations return bad_dble
! 
! If input values of x1 and x2 do not bracket the root, print warning and
! return with best estimate of root.
!
implicit none 
integer, intent(in)  :: npar, verbose
double precision, intent(in) :: x1, x2, tol, par(npar)
interface
  function func(fx,fnpar, fpar, fverbose)
  double precision             :: func
  double precision, intent(in) :: fx
  integer, intent(in)  :: fnpar
  double precision, intent(in) :: fpar(fnpar)
  integer, intent(in)  :: fverbose
  end function func
  function dfunc(fx,fnpar, fpar, fverbose)
  double precision             :: dfunc
  double precision, intent(in) :: fx
  integer, intent(in)  :: fnpar
  double precision, intent(in) :: fpar(fnpar)
  integer, intent(in)  :: fverbose
  end function dfunc
end interface

integer, parameter :: maxit = 1024  ! Maximum number of iterations
double precision, parameter :: ftol = 4.d0*epsilon(0.d0)
integer :: j
double precision :: df, dx, dxold, f, fh, fl, temp, xh, xl

rtsafe = bad_dble

fl = func(x1, npar, par, verbose-1)
fh = func(x2, npar, par, verbose-1)

if (((fl > 0.0d0).and.(fh > 0.0d0)).or.((fl < 0.0d0).and.(fh < 0.0d0))) then
if (verbose >= v_error) then
    print *, 'rtsafe: input values do not bracket root'
    print *, x1, x2
    print *, fl,fh
  endif
  return
endif

if (fl == 0.0d0) then
  rtsafe = x1
  return
else if (fh == 0.0d0) then
  rtsafe = x2
  return
else if (fl < 0.0d0) then
  xl = x1
  xh = x2
else
  xh = x1
  xl = x2
endif

rtsafe = 0.5d0*(x1+x2)
dxold  =abs(x2-x1)
dx=dxold
f= func(rtsafe,npar, par, verbose-1)
df = dfunc(rtsafe,npar, par, verbose-1)
do j=1,maxit

  if( ((rtsafe-xh)*df-f)*((rtsafe-x1)*df-f) > 0.0d0 .or.  &
      abs(2.0d0*f) > abs(dxold*df) ) then
    ! Newton out-of-range or not decreasing fast enough - bisect
    dxold = dx
    dx = 0.5d0*(xh-xl)
    rtsafe = xl+dx
    if(xl == rtsafe) return ! Change in root is negligible.
  else
    dxold = dx
    dx=f/df
    temp=rtsafe
    rtsafe=rtsafe-dx
    if(temp == rtsafe) return
  endif
  if (abs(dx) < tol) return  ! Convergence criterion.
  f= func(rtsafe,npar, par, verbose-1)
  ! Added a test here to return if abs(f) < ftol
  if (abs(f) < ftol) return  
  df = dfunc(rtsafe,npar, par, verbose-1)
  if( f < 0.d0) then
    xl = rtsafe
  else 
    xh=rtsafe
  endif
end do
if (verbose >= v_error) then
  print *, 'rtsafe: Failed to converge'
  print *, xl, xh
  print *, x1, x2
  print *, f, df,ftol
endif
return
end function rtsafe

!-----------------------------------------------------------------------------

double complex function laguerre(a,m,x)
! Given a degree m and the complex coefficients a(1:m+1) of the polynomial
! sum(i=1,m+1) a(i) * x**(i-1), and given a complex value x, this routine
! improves x by Laguerre's method until it converges, within the achievable
! round-off limit, to a root of the given polynomial. 
!  If the algorithm fails to converge return value is (bad_dble, bad_dble)
! 
! Follows the algorithm described in Numerical Recipes in Fortran 90
! (Press et al., http://www.nr.com)
!
integer, intent(in)        :: m      ! Degree of the polynomial
double complex, intent(in) :: a(m+1) ! Polynomial coefficients (ascending power)
double complex, intent(in) :: x      ! Initial/final estimate of root

! Parameters
integer, parameter  :: mr = 8         ! Number of fractional values to test
integer, parameter  :: mt = 10        ! Step frequency for limit cycle tests
integer, parameter  :: maxit = mr*mt  ! Maximum number of iterations

! Local variables
integer     :: iter, j
double precision :: abx, abp, abm, err, frac(mr)
double complex :: dx, x1, b, d, f, g, h, sq, gp, gm, g2
save frac
! Fractions used to break limit cycle
data frac /0.5d0, 0.25d0, 0.75d0, 0.13d0, 0.38d0, 0.62d0, 0.88d0, 1.0d0/

laguerre = x

do iter=1,maxit
  b=a(m+1)
  err=abs(b)
  d=dcmplx(0.0d0,0.0d0)
  f=dcmplx(0.0d0,0.0d0)
  abx=abs(laguerre)
  do j=m,1,-1
    f = laguerre*f + d
    d = laguerre*d + b
    b = laguerre*b + a(j)
    err = abs(b) + abx*err
  end do
  err = epsilon(0.0d0)*err
  if (abs(b) <= err) then
    return
  else
    g = d/b
    g2=g**2
    h = g2-2.0d0*f/b
    sq = sqrt((m-1)*(m*h-g2))
    gp = g + sq
    gm = g - sq
    abp = abs(gp)
    abm = abs(gm)
    if (abp < abm) gp = gm
    if (max(abp,abm) > 0.0d0) then
      dx = m/gp
    else
      dx = exp(dcmplx(log(1.0d0+abx),dble(iter)))
    endif
  endif
  x1 = laguerre - dx
  if (laguerre == x1) return
  if (mod(iter,mt) /= 0) then
    laguerre = x1
  else
    laguerre = laguerre - dx*frac(iter/mt)
  endif
end do
laguerre = dcmplx(bad_dble,bad_dble)
return
end function laguerre


!-----------------------------------------------------------------------------

function affine(x,tr,inverse) result(xt)
implicit none
! Apply the affine transformation tr(2,3) to an array of points x(2).
!  The result is
!   xt(1,i) = tr(1,3) + tr(1,1)*x(1,i) + tr(1,2)*x(2,i)
!   xt(2,i) = tr(2,3) + tr(2,1)*x(1,i) + tr(2,2)*x(2,i)
!
! If the argument inverse is present and .true. then the results is
!
!   (  tr(2,2)*(x(1,i)-tr(1,3)) - tr(1,2)*(x(2,i)-tr(2,3)) )/D
!   ( -tr(2,1)*(x(1,i)-tr(1,3)) + tr(1,1)*(x(2,i)-tr(2,3)) )/D
! where D = tr(1,1)*tr(2,2) -  tr(1,2)*tr(2,1) 
! N.B. if D == 0, all values are returned as bad_dble.
!
double precision, intent(in)         :: x(2)
double precision, intent(in)         :: tr(2,3)
double precision                     :: xt(2)
logical, optional,intent(in) :: inverse 

! Local variables
double precision :: d
logical  :: inverse_l

if (present(inverse)) then
  inverse_l = inverse
else
  inverse_l = .false.
endif

if (inverse_l) then
  d = tr(1,1)*tr(2,2) -  tr(1,2)*tr(2,1) 
  if (d == 0.0d0) then
    xt(:) = bad_dble
    return
  endif
  xt(1) = ( tr(2,2)*(x(1)-tr(1,3)) - tr(1,2)*(x(2)-tr(2,3)))/d
  xt(2) = (-tr(2,1)*(x(1)-tr(1,3)) + tr(1,1)*(x(2)-tr(2,3)))/d
else
  xt(1) = tr(1,3) + tr(1,1)*x(1) + tr(1,2)*x(2)
  xt(2) = tr(2,3) + tr(2,1)*x(1) + tr(2,2)*x(2)
endif
end function affine

!-----------------------------------------------------------------------------

double precision function determinant(a,n) 
! Determinant of a 2x2 or 3x3 matrix
integer,  intent(in) :: n
double precision,intent(in) :: a(n,n)

select case (n)
case(0)
  determinant = 0.d0
case(1)
  determinant = a(1,1)
case(2)
  determinant = a(1,1)*a(2,2) - a(1,2)*a(2,1)
case(3) 
  determinant = ( (a(3,1)*a(2,2)-a(2,1)*a(3,2))   &
                 *(a(2,1)*a(1,3)-a(1,1)*a(2,3)) ) &
              - ( (a(2,1)*a(1,2)-a(1,1)*a(2,2))   &
                 *(a(3,1)*a(2,3)-a(2,1)*a(3,3)) )
case default
  print *,'determinant: n>3 not implemented.'
  stop
end select
return
end function determinant

!-----------------------------------------------------------------------------

function deflate_poly(a, n, c) result(b)
! Divide a polynomial by a monomial factor (x-c)
! a(1:n) are the coefficients of the polynomial in ascending power order
! b(1:n-1) are the coefficients of the divided polynomial
! b(n) is the remainder
integer, intent(in) :: n
double precision,intent(in) :: a(n), c
double precision :: b(n)

! Local variables
integer  :: i
double precision :: r

r = a(n)
do i=n-1,1,-1
  b(i) = r
  r = a(i)+r*c
end do
b(n) = r
return
end function deflate_poly

!-----------------------------------------------------------------------------

double precision function midangle(th0,th1, reflex)
! Bisector of the angles th0 and th1 in the range 0 -> twopi
! If reflex is present and .true., return the bisector the forms
! an reflex angle with the two input angles.
implicit none
double precision, intent(in) :: th0, th1
logical, intent(in), optional :: reflex
! Local variables
logical reflex_l

if (present(reflex)) then
  reflex_l = reflex
else
  reflex_l = .false.
endif

if (reflex_l) then

  if (abs(th1-th0).gt.pi) then
    midangle = 0.5d0*(th0+th1)
  else 
    midangle = mod(pi+0.5d0*(th0+th1), twopi)
  endif

else

  if (abs(th1-th0).gt.pi) then
    midangle = 0.5d0*(th0+th1) - pi
    if (midangle.lt.0.d0) midangle = midangle + twopi
  else 
    midangle = 0.5d0*(th0+th1)
  endif

endif

end function midangle
                                                                        
!-------------------------------------------------------------------------------

double precision function triangle_area(p_1,p_2,p_3,signed)
! Area of a triangle defined by three points
! p_1 = (x_1, y_1), p_2 = (x_2, y_2) and p_3 = (x_3, y_3)
! If signed=.true. then the result is positive/negative according to whether
! the points are clockwise/counterclockwise
implicit none
double precision, intent(in) :: p_1(2), p_2(2),p_3(2)
logical, optional,intent(in) :: signed

! Local variables
logical :: signed_l

if (present(signed)) then
  signed_l = signed
else
  signed_l = .false.
endif

if (signed_l) then
 triangle_area = 0.5d0*(p_1(1)*(p_2(2)-p_3(2)) &
               - p_2(1)*(p_1(2)-p_3(2)) + p_3(1)*(p_1(2)-p_2(2)))
else
 triangle_area = abs(0.5d0*(p_1(1)*(p_2(2)-p_3(2)) &
               - p_2(1)*(p_1(2)-p_3(2)) + p_3(1)*(p_1(2)-p_2(2))))
endif

return
end function triangle_area

!-----------------------------------------------------------------------------

function ucircle_line_intersect(qline) result(t)
!
! Intersection of a line with unit circle centred on the origin
!
! The line is defined by the parametric equation
!  (x, y) = (x_0, y_0) + t*(dx, dy) 
!
! qline = [x_0, y_0, dx, dy]
!
! Values of t for which the line intersects the unit circle are returned in 
! t(1:2)
!     
! If there is only one intersection then t(1) = t(2)
!     
! If there are no intersections then t(1) = -huge(0.d0) and t(2) = huge(0.d0)
!
implicit none
double precision, intent(in) :: qline(4) ! Parametric equation of the line
double precision :: t(2)                 ! Line parameter values at intercept

! Local variables
double precision :: a,b,c,d,r,w

a = qline(3)**2 + qline(4)**2
b = 2.0d0*(qline(1)*qline(3) + qline(2)*qline(4))
c = qline(1)**2 + qline(2)**2 - 1.0d0
d = b**2 - 4.0d0*a*c
if (d < 0.0) then
  t(1) = -huge(0.d0)
  t(2) =  huge(0.d0)
else
  r = sqrt(d)
  w = -(b + sign(r,b))
  t(1) = 2.0d0*c/w
  t(2) = 0.5d0*w/a
endif
return 
end function ucircle_line_intersect
                                                                        
!-------------------------------------------------------------------------------

double precision function arc_area(p_1,p_2)
! Area of the arc defined by a unit circle centred on the origin and the lines
! between two points on its circumference, p_1 = (x_1, y_1) and p_2 = (x_2, y_2)
! N.B. requires hypot(x_1, y_2) = hypot(x_2, y_2) = 1 - not checked!
implicit none
double precision, intent(in) :: p_1(2), p_2(2)

! Local variables
double precision :: t_1, t_2, dt
double precision, parameter ::p_0(2) = (/0.d0,0.d0/)

t_1 = atan2(p_1(2),p_1(1))
t_2 = atan2(p_2(2),p_2(1))
dt = mod(abs(t_1 - t_2), twopi)
if (dt > pi) dt = twopi - dt
arc_area = 0.5d0*dt - triangle_area(p_0,p_1,p_2)
return
end function arc_area
                                                                        
!-------------------------------------------------------------------------------

subroutine circles_intersect(circle_1, circle_2, xy, circles_flags)
! Number and positions of intersections points of two circles on a plane
implicit none
double precision, intent(in)  :: circle_1(3) ! (x_c, y_c) and radius of circle 1
double precision, intent(in)  :: circle_2(3) ! (x_c, y_c) and radius of circle 2
double precision, intent(out) :: xy(2,2)     ! (x,y) for intersections 
integer, intent(out)  :: circles_flags ! Bit mask for flags
! circles_flags is a bit mask with the following flags
! b_circles_no_overlap  (but may have a tangent point)
! b_circles_1_inside_2  (and may have a tangent point)
! b_circles_2_inside_1  (and may have a tangent point)
! b_circles_tangent     (there is a tangent point, xy(1:2,1) = xy(1:2,2)
! b_circles_intersect   (no tangent points, two distinct intersection points)
! b_circles_identical   (xy not defined in this case)
! b_circles_invalid     (radius <=  0 for either circle)
!
! To test these bit flags, use btest(circles_flags, <bit flag name>)

! Local variables
double precision :: a, dx, dy, d, h, rx, ry, x2, y2, r_1, r_2

r_1 = circle_1(3)
r_2 = circle_2(3)
circles_flags = 0 
if ((r_1 <= 0.0d0).or.(r_2 <= 0.0d0)) then
  xy(:,:) = not_set_dble
  circles_flags = ibset(circles_flags, b_circles_invalid)
  return
endif

dx = circle_2(1) - circle_1(1)
dy = circle_2(2) - circle_1(2)
d = hypot(dx,dy)

if ((d == 0.0d0).and.(r_1 == r_2)) then
  xy(:,:) = not_set_dble
  circles_flags = ibset(circles_flags, b_circles_identical)
  return
endif

if (d > (r_1 + r_2)) then
  xy(:,:) = not_set_dble
  circles_flags = ibset(circles_flags, b_circles_no_overlap)
  return
endif

if (d < abs(r_2 - r_1)) then
  xy(:,:) = not_set_dble
  if (r_1 < r_2) then
    circles_flags = ibset(circles_flags, b_circles_1_inside_2)
  else
    circles_flags = ibset(circles_flags, b_circles_2_inside_1)
  endif
  return
endif

a = (r_1**2 - r_2**2 + d**2) / (2.0d0 * d)
x2 = circle_1(1) + dx * a/d
y2 = circle_1(2) + dy * a/d

if (d == (r_1 + r_2)) then
  xy(1:2,1) = (/ x2, y2 /)
  xy(1:2,2) = (/ x2, y2 /)
  circles_flags = ibset(circles_flags, b_circles_tangent)
  circles_flags = ibset(circles_flags, b_circles_no_overlap)
  return
endif

if (d == abs(r_1 - r_2)) then
  xy(1:2,1) = (/ x2, y2 /)
  xy(1:2,2) = (/ x2, y2 /)
  circles_flags = ibset(circles_flags, b_circles_tangent)
  if (r_1 < r_2) then
    circles_flags = ibset(circles_flags, b_circles_1_inside_2)
  else
    circles_flags = ibset(circles_flags, b_circles_2_inside_1)
  endif
  return
endif

h = sqrt(r_1**2 - a**2)
rx = -dy * h/d
ry =  dx * h/d
xy(1:2,1) = (/ x2 + rx, y2 + ry /)
xy(1:2,2) = (/ x2 - rx, y2 - ry /)
circles_flags = ibset(circles_flags, b_circles_intersect)

return

end subroutine circles_intersect
                                                         
                                                                        
!-------------------------------------------------------------------------------

subroutine sph_circ_intersect(circle_1, circle_2, xpts, circles_flags)
! Number and positions of intersections points of two circles on a sphere
! See https://gis.stackexchange.com/questions/48937/calculating-intersection-of-two-circles
! All input/output angles are in radians, including the radii of the circles.
implicit none
double precision, intent(in)  :: circle_1(3) ! (lon, lat) and radius of circle 1
double precision, intent(in)  :: circle_2(3) ! (lon, lat) and radius of circle 2
double precision, intent(out) :: xpts(2,2)   ! (lon,lat) for intersections 
integer, intent(out)  :: circles_flags ! Bit mask for flags
! circles_flags is a bit mask with the following flags
! b_circles_no_overlap  (but may have a tangent point)
! b_circles_1_inside_2  (and may have a tangent point)
! b_circles_2_inside_1  (and may have a tangent point)
! b_circles_tangent     (there is a tangent point, xy(1:2,1) = xy(1:2,2)
! b_circles_intersect   (no tangent points, two distinct intersection points)
! b_circles_identical   (xy not defined in this case)
! b_circles_invalid     (parameters invalid for either circle)
!
! To test these bit flags, use btest(circles_flags, <bit flag name>)

! Local variables
double precision :: x_1(3), x_2(3), n(3), x_0(3)
double precision :: r_1, r_2, q, d, a, b, t

circles_flags = 0 
r_1 = circle_1(3)
r_2 = circle_2(3)

if ( (r_1 <= 0.0d0).or.(r_1 >= halfpi).or. &
     (r_2 <= 0.0d0).or.(r_2 >= halfpi) ) then
  xpts(:,:) = not_set_dble
  circles_flags = ibset(circles_flags, b_circles_invalid)
  return
endif

d = angular_distance(circle_1(1),circle_1(2),circle_2(1),circle_2(2))

if ((d <  epsilon(0.0d0)).and.(abs(r_1 - r_2)< epsilon(0.0d0))) then
    xpts(:,:) = not_set_dble
    circles_flags = ibset(circles_flags, b_circles_identical)
    return
endif

if (d > (r_1+r_2)) then
  xpts(:,:) = not_set_dble
  circles_flags = ibset(circles_flags, b_circles_no_overlap)
  return
endif

if (d < abs(r_2 - r_1)) then
  xpts(:,:) = not_set_dble
  if (r_1 < r_2) then
    circles_flags = ibset(circles_flags, b_circles_1_inside_2)
  else
    circles_flags = ibset(circles_flags, b_circles_2_inside_1)
  endif
  return
endif


x_1 = (/cos(circle_1(1))*cos(circle_1(2)), &
        sin(circle_1(1))*cos(circle_1(2)), &
        sin(circle_1(2)) /)
x_2 = (/cos(circle_2(1))*cos(circle_2(2)), &
        sin(circle_2(1))*cos(circle_2(2)), &
        sin(circle_2(2)) /)

q = dot_product(x_1, x_2)
a = (cos(r_1) - cos(r_2)*q) / (1.0d0 - q**2)
b = (cos(r_2) - cos(r_1)*q) / (1.0d0 - q**2)

n = cross3(x_1, x_2)

x_0 = a*x_1 + b*x_2

t = sqrt((1.0d0 - min(1.0d0,dot_product(x_0,x_0)))/dot_product(n,n))

if (t == 0.0d0) then
    circles_flags = ibset(circles_flags, b_circles_tangent)
endif

x_1 = x_0 + t*n
x_2 = x_0 - t*n

xpts(1,1) =  atan2(x_1(2), x_1(1))
xpts(2,1) =  atan2(x_1(3), hypot(x_1(1), x_1(2))) 
xpts(1,2) =  atan2(x_2(2), x_2(1))
xpts(2,2) =  atan2(x_2(3), hypot(x_2(1), x_2(2))) 
      
return

end subroutine sph_circ_intersect
                                                         
!-------------------------------------------------------------------------------

double precision function distance(xy_1, xy_2) 
implicit none
double precision, intent(in) :: xy_1(2), xy_2(2)
distance = hypot(xy_1(1)-xy_2(1), xy_1(2)-xy_2(2))
return
end function distance

!-------------------------------------------------------------------------------

double precision function angular_distance(lon_1, lat_1, lon_2, lat_2)
! Great circle distance in radians between points (lon1, lat1) and 
! (lon2, lat2), also in radians. 
! See https://www.movable-type.co.uk/scripts/gis-faq-5.1.html
implicit none
double precision, intent(in)  :: lon_1, lat_1, lon_2, lat_2

double precision :: dlon, dlat, a
dlon = lon_2-lon_1
dlat = lat_2-lat_1
a = sin(0.5d0*dlat)**2 + cos(lat_1)*cos(lat_2)*sin(0.5d0*dlon)**2
angular_distance = 2.0d0 * asin(min(1.0d0,sqrt(a)))
return
end function angular_distance

!----------------------------------------------------------------------

function parfit(dx,dy) result (a)

double precision, intent(in) :: dx(3),dy(3)
double precision :: a(3)
! Returns the coefficients of the unique parabola intersecting three
! points. Calculation done in double precision.
! If no such parabola exists, all three coefficients are set to bad_dble
!
!-----------------------------------------------------------------------

double precision :: dv1,dv2,d
      
dv1 = dx(3) * (dx(3) - dx(1))
dv2 = dx(2) * (dx(2) - dx(1))
d =  dv1*(dx(1)-dx(2)) - dv2*(dx(1)-dx(3)) 
if (abs(d) < epsilon(0.0d0)) then
  a(1:3) = bad_dble
  return
endif
a(1) = ( dv1*(dx(1)*dy(2)-dx(2)*dy(1)) -  &
         dv2*(dx(1)*dy(3)-dx(3)*dy(1)) )/d

dv1 = dx(3)*dx(3)-dx(1)*dx(1)
dv2 = dx(2)*dx(2)-dx(1)*dx(1)
d =  dv1*(dx(2)-dx(1)) - dv2*(dx(3)-dx(1))  
if (abs(d) < epsilon(0.0d0)) then
  a(1:3) = bad_dble
  return
endif
a(2) = ( dv1*(dy(2)-dy(1)) - dv2*(dy(3)-dy(1)) )/d

dv1 = dx(3) - dx(1)
dv2 = dx(2) - dx(1) 
d =  dv1*(dx(2)*dx(2)-dx(1)*dx(1)) - dv2*(dx(3)*dx(3)-dx(1)*dx(1))
if (abs(d) < epsilon(0.0d0)) then
  a(1:3) = bad_dble
  return
endif
a(3) = ( dv1*(dy(2)-dy(1)) - dv2*(dy(3)-dy(1)) )/d
      
return

end function parfit

!-----------------------------------------------------------------------------      

double precision function unitfunc(u,v,npar,par,verbose) 
implicit none
integer, intent(in)  :: npar, verbose
double precision, intent(in) :: u,v,par(npar)
unitfunc = 1.0d0
return
end function unitfunc

!----------------------------------------------------------------------

integer function verbose_for_calls(verbose)
! When debugging, reduce verbose by one for all function/subroutine calls.
integer, intent(in) :: verbose
if (verbose >= v_debug) then
  verbose_for_calls = verbose -1
else
  verbose_for_calls = verbose
endif
end function verbose_for_calls

!----------------------------------------------------------------------

function eval_poly_deriv(x,np,p,verbose) result(px_dpdx)
! Evaluate a polynomial and its derivative
! Argument list suitable for calling by rtsafe or zbrent for root finding.
implicit none
double precision, intent(in) :: x     ! Value at which to evaluate polynomial
integer, intent(in)          :: np    ! No. of coefficients
double precision, intent(in) :: p(np) ! Array of coefficients, ascending powers
integer, intent(in)          :: verbose 
double precision             :: px_dpdx(2) ! p(x) and dp/dx

! Local variables
integer :: i
double precision :: px,dpdx

if (verbose >= v_debug) print *,'eval_poly_deriv: np,p =',real(np),real(p(1:np)) 
select case (np)
case(:0) 
  px = bad_dble
  dpdx = bad_dble
case(1) 
  px = p(1)
  dpdx = 0.d0
case(2) 
  px = x*p(2)+p(1)
  dpdx = p(2)
case(3) 
  px = x*(x*p(3)+p(2))+p(1)
  dpdx = 2.d0*x*p(3)+p(2)
case(4) 
  px = x*(x*(x*p(4)+p(3))+p(2))+p(1)
  dpdx = x*(x*p(4)*3.d0+p(3)*2.d0)+p(2)
case(5) 
  px = x*(x*(x*(x*p(5)+p(4))+p(3))+p(2))+p(1)
  dpdx = x*(x*(x*p(5)*4.d0+p(4)*3.d0)+p(3)*2.d0)+p(2)
case  default
  px = p(np)*x + p(np-1)
  dpdx = p(np)
  do i=np-2,1,-1
    dpdx = dpdx*x*px
    px = px*x + p(i)
  end do
end select
px_dpdx(1) = px
px_dpdx(2) = dpdx
return
end function eval_poly_deriv
                                                                        
!----------------------------------------------------------------------

double precision function brent(ax, bx, cx, f, npar, par, tol, xmin, verbose)

! Given a function F, and given a bracketing triplet of abscissas AX,BX,CX
! such that BX is between AX and CX, and F(BX) is less than both F(AX)
! and F(CX), this routine isolates the minimum to a fractional precision
! of about TOL using Brent's method. The abscissa of the minimum is returned
! as XMIN, and the minimum function value is returned as BRENT, the returned
! function value

implicit  none
integer, intent(in) :: npar, verbose
double precision,intent(in) :: ax, bx, cx, tol, par(:)
double precision,intent(out) :: xmin

! Maximum allowed number of iterations, golden ratio and a small
! number to prevent trying for a fractional accuracy in a minimum
! which is zero
double precision, parameter :: cgold = 0.3819660d0
double precision, parameter :: zeps = 1.0d-10

integer, parameter  ::  itmax=100 
integer :: iter, verbose1
double precision :: w, x, e, fx, fv, fw, fu, u, v
double precision :: xm, tol1, tol2, r, q, p, a, b, d, etemp

interface
  function f(fx,fnpar, fpar, fverbose)
  double precision             :: f
  double precision, intent(in) :: fx
  integer, intent(in)  :: fnpar
  double precision, intent(in) :: fpar(fnpar)
  integer, intent(in)  :: fverbose
  end function f
end interface

verbose1 = verbose_for_calls(verbose)
brent = not_set_dble
d = not_set_dble
a = min(ax, cx)
b = max(ax, cx)
v = bx
w = v
x = v
e = 0.0d0
fx = f(x,npar,par,verbose1)
fv = fx
fw = fx
do iter = 1, itmax
  xm = 0.5d0*(a+b)
  tol1 = tol*abs(x)+zeps
  tol2 = 2.0d0*tol1
  if(abs(x-xm) <= (tol2-0.5d0*(b-a))) then
    xmin = x
    brent = fx
    return
  endif
  if(abs(e) > tol1) then
    r = (x-w)*(fx-fv)
    q = (x-v)*(fx-fw)
    p = (x-v)*q - (x-w)*r
    q = 2.0d0*(q-r)
    if(q > 0.0d0) p = - p
    q = abs(q)
    etemp = e
    e = d
    if(abs(p) >= abs(.5*q*etemp).or. &
       p <= q*(a-x) .or. p >= q*(b-x)) then
       e = merge(a-x, b-x, p >= q*(b-x))
       d = cgold*e
     else
       d = p/q
       u=x+d
       if ( u-a < tol2 .or. b-u < tol2) d = sign(tol1, xm-x)
     endif
   else
     e = merge(a-x, b-x, x >= xm)
     d = cgold*e
   endif
   u = merge(x+d, x+sign(tol1,d), abs(d) >= tol1)
   fu = f(u,npar,par,verbose1)
   if (fu <= fx) then
     if ( u >= x) then
       a = x
     else
       b = x
     endif
     call shft(v,w,x,u)
     call shft(fv,fw,fx,fu)
   else
     if (u < x) then
       a = u
     else
       b = u
     endif
     if (fu <= fw .or. w == x) then
       v=w
       fv=fw
       w=u
       fw=fu
     else if (fu <= fv .or. v == x .or. v == w) then
       v=u
       fv=fu
     endif
   endif
 end do
 if (verbose >= v_warn) print *,'brent: failed to converge'
 contains

 subroutine shft(a,b,c,d)
 double precision, intent(out) :: a
 double precision, intent(inout) :: b, c
 double precision, intent(in) :: d
 a=b
 b=c
 c=d
 end subroutine shft
 end function brent

!----------------------------------------------------------------------

double precision function acot(x)
implicit none
double precision, intent(in) :: x  ! Input value

if (x < 0.d0 ) then
  acot = -halfpi - atan(x)
else
  acot =  halfpi - atan(x)
endif 
end function acot

!----------------------------------------------------------------------

function cross3(a, b) result(x)
implicit none
double precision, intent(in) :: a(3), b(3)
double precision :: x(3)

x = (/ a(2) * b(3) - a(3) * b(2), &
       a(3) * b(1) - a(1) * b(3), &
       a(1) * b(2) - a(2) * b(1) /)

end function cross3

!-----------------------------------------------------------------------------
end module utils
