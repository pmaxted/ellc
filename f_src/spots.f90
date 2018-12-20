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

module spots
! Module for dealing with star spots in ellc
!
! HISTORY
! -------
! Jan 2016
! First version 
! p.maxted@keele.ac.uk
!

use constants
use utils     
use ellipse

implicit none

contains

subroutine eker(l,b,i,r,a,u1,u2,phi,nphi,df,ii,ifail)
implicit none
integer,intent(in)  :: nphi
double precision, intent(in) :: l,b,i,phi(nphi),r,a,u1,u2
double precision, intent(out) :: df(nphi)
integer,intent(out) :: ii(nphi), ifail
!
! Analytic light curve for a rotating star with a single circular spot 
!
! Implementation of formulae from Eker, Z., 1994ApJ...420..373E with 
! correction for equation for K_3 from erratum 1994ApJ...430..438E.
!
! Input values
!  l        = longitude of spot centre (radians)
!  b        = latitude of spot centre (radians)
!  i        = inclination of rotation axis (radians)
!  r        = angular radius of spot (radians)
!  a        = Spot contrast ratio (a=Is/Ip).
!  u1       = linear limb darkening coefficient.
!  u2       = quadratic limb darkening coefficient.
!  phi(nphi) = array of rotation phase values in radians, i.e. twopi*(t-t_0)/P
!  nphi      = number of phase values
!
! N.B. Eker uses a limb darkening law of the form I_0[1 - u_1.mu + u_2.mu^2], 
! i.e., u_2 is the negative of the normal quadratic limb darkening coefficient.
!
! Output values
!  df(nphi)  = light curve
!  ii(nphi)  = spot position flag
!  ifail    = status flag
!
! Return values in the array df are (Fp + Fs)/F where
!  Fp is the flux from the unspotted photosphere
!  Fs is the flux from the spot
!  F is the flux from the star without a spot.
! 
! Spot position flag ii() is as follows
!   0 = spot not visible
!   1 = Spot is on the limb and less than half the spot is visible
!   2 = Spot is on the limb and more than half the spot is visible.
!   3 = Spot is completely visible
!
! Return value of ifail is the sum of the following flags.
!   0 => all ok
!   1 => r >= pi/2
!   2 => a < 0.0
!   4 => nphi < 1
!

! Local variables
integer :: iphi
double precision :: cosisinb,sinicosb
double precision :: sinr2,sinr3,sinr4,sin2r
double precision :: cosr,cos2r,cosr3,cosr4
! qn = I_n/pi, ql = I_l/pi, qq = I_q/pi
double precision :: qn,ql,qq
double precision :: th0,Cl,Cq,K1,K2,K3,K4
double precision :: t13b1,t13b2,t13c1,t13c2
double precision :: costh0,sinth0,sinth02,tanth0
double precision :: phi0,r0,sinphi0,cosphi0,sinr0,cosr0,sin2r0,cos2r0
double precision :: f0,fn,fl,fq

ifail = 0 
if (r >= halfpi) ifail = ifail + 1
if (a < 0.0d0)  ifail = ifail + 2
if (nphi < 1) ifail = ifail + 4
if (ifail /= 0) return

! Take constant terms in the calculation out of the do loop

cosisinb = cos(i)*sin(b)
sinicosb = sin(i)*cos(b)
cosr = cos(r)
cosr3 = cosr**3
cosr4 = cosr**4
sin2r = sin(2.0d0*r)
cos2r = cos(2.0d0*r)
sinr2 = 1.0d0 - cosr**2
sinr3  = sin(r)**3
sinr4  = sinr2**2
f0 = (a-1.0d0)/(1.0d0-u1/3.0d0+u2/6.0d0)
fn = (1.0d0-u1+u2)*f0
fl = (u1-2.0d0*u2)*f0
fq = u2*f0
t13b1 = 2.0d0*(1.0d0-cosr3)/3.0d0
t13b2 = cosr*sinr2
t13c1 = 0.5d0*(1.0d0-cosr4)
t13c2 = 0.75d0*sinr4

do iphi=1,nphi

! (20)
  costh0 = cosisinb+sinicosb*cos(phi(iphi)-l)
  th0 = acos(costh0)
  sinth02 = 1.0d0 - costh0**2
  sinth0 = sqrt(sinth02)
  tanth0 = sinth0/costh0

  if (th0 >= (halfpi+r)) then ! Spot is completely around the back side

    df(iphi) = 1.0d0
    ii(iphi) = 0 

  else ! Spot is visible or partially visible

    ii(iphi) = 3
    if ((u1 == 0.0d0).and.(u2 == 0.0d0)) then
      qn = sinr2*costh0
      ql = 0.0d0
      qq = 0.0d0  
    else if (u2 == 0.0d0) then
      qn = sinr2*costh0
      ql = t13b1 - t13b2*sinth02
      qq = 0.0d0  
    else
      qn = sinr2*costh0
      ql = t13b1 - t13b2*sinth02
      qq = t13c1*costh0**3 + t13c2*costh0*sinth02
    endif

    if (th0 > halfpi) then
    ! Spot is on the limb and less than half the spot is visible...

      ii(iphi) = 1
      cosphi0 = -1.0d0/(tanth0*tan(r))
      phi0 = acos(cosphi0)
      sinphi0 = sin(phi0)
      qn = (phi0*qn - asin(cosr/sinth0) - 0.5d0*sinth0*sinphi0*sin2r)/pi + 0.5d0
      if ((u1 /= 0.0d0).or.(u2 /= 0.0d0)) then
        r0 = abs(th0-halfpi)
        sinr0 = sin(r0)
        cosr0 = cos(r0)
        sin2r0 = sin(2.0d0*r0)
        cos2r0 = cos(2.0d0*r0)
        ! (19a)
        ql = (phi0/3.0d0*(cosr3-cosr0**3)* (1.0d0-3.0d0*costh0**2) -  &
             (phi0 + sinphi0*cosphi0)*(cosr - cosr0)*sinth02 +        &
             4.0d0/3.0d0*sinphi0*(sinr3-sinr0**3)*sinth0*costh0 +     &
             sinphi0*cosphi0/3.0d0*(cosr3 - cosr0**3)*sinth02)/pi
        if (u2 /= 0.0d0) then
          K1 = 0.25d0*phi0*(cosr0**4-cosr4)
          K2 = -0.125d0*sinphi0*(r0 - r + 0.5d0*(sin2r*cos2r-sin2r0*cos2r0))
          K3 = 0.125d0*(phi0+sinphi0*cosphi0)*(sinr4-sinr0**4)
          K4 = (sinphi0-sinphi0**3/3.0d0)*(0.375d0*(r-r0) +   &
               0.0625d0*(sin2r*(cos2r-4.0d0)-sin2r0*(cos2r0-4.0d0)))
          qq = (2.0d0*costh0**3*K1 +     &
                6.0d0*costh0*sinth0* (costh0*K2 + sinth0*K3) +   &
                2.0d0*sinth0**3*K4 )/pi
        endif
      endif

      ! Spot is on the limb and more than half the spot is visible...
    else if (th0 > (halfpi-r)) then

      ii(iphi) = 2
      cosphi0 = -1.0d0/(tanth0*tan(r))
      phi0 = acos(cosphi0)
      sinphi0 = sin(phi0)
      qn = (phi0*qn - asin(cosr/sinth0) - 0.5d0*sinth0*sinphi0*sin2r)/pi + 0.5d0
      if ((u1 /= 0.0d0).or.(u2 /= 0.0d0)) then
        r0 = abs(th0-halfpi)
        sinr0 = sin(r0)
        cosr0 = cos(r0)
        sin2r0 = sin(2.0d0*r0)
        cos2r0 = cos(2.0d0*r0)
        ! (18a)
        Cl = (pi-phi0)/3.0d0*(cosr3-cosr0**3)*(1.0d0-3.0d0*costh0**2) -  &
             (pi-phi0-sinphi0*cosphi0)*(cosr-cosr0)*sinth02 -   &
             4.0d0/3.0d0*sinphi0*(sinr3-sinr0**3)*sinth0*costh0 -   &
             sinphi0*cosphi0*(cosr3-cosr0**3)*sinth02/3.0d0
        ql = ql - Cl/pi
        if (u2 /= 0.0d0) then
          K1 = 0.25d0*(pi-phi0)*(cosr0**4-cosr4)
          K2 = 0.125d0*sinphi0*(r0 - r + 0.5d0*(sin2r*cos2r-sin2r0*cos2r0))
          K3 = 0.125d0*(pi-phi0+sinphi0*cosphi0)*(sinr4-sinr0**4)
          K4 = -(sinphi0-sinphi0**3/3.0d0)*(0.375d0*(r-r0) + &
                 0.0625d0*(sin2r*(cos2r-4.0d0)-sin2r0*(cos2r0-4.0d0)))
          ! (18b)
          Cq = 2.0d0*costh0**3*K1 + 6.0d0*costh0*sinth0* &
               (costh0*K2 + sinth0*K3) + 2.0d0*sinth0**3*K4
          qq = qq - Cq/pi
        endif
      endif
    endif
    ! (12c) (which reduces to (12b) and (12a) if u1=0 and/or u2=0)
    df(iphi) = 1.0d0 + fn*qn + fl*ql + fq*qq
  endif

enddo

end subroutine eker

!------------------------------------------------------------------------------

double precision function spot_limb_eclipse(ellipse, ellipse_s, xy_tng, &
                                            large, verbose)
! Fraction of a spot on the limb of a star that is eclipsed by an ellipse.
!
! Returns bad_dble is algorithm fails
! All input angles are in radians
!
implicit none
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(in) :: ellipse(n_ell_par) ! Occulting ellipse 
double precision, intent(in) :: ellipse_s(n_ell_par) ! Ellipse for spot
! (x,y) positions of tangent points between projected spot and circle
double precision, intent(in) :: xy_tng(2,2) 
logical, intent(in) :: large                       ! See below ...
integer, intent(in) :: verbose

! A circular spot on the limb in projection is seen as an ellipse that is
! interior to the circle and tangent to it at two points.
!
! large  = .true. to take the larger lens-shaped area between the tangent
!          ellipse and the circle.
!        = .false. to take the smaller crescent-shaped area between the
!          tangent ellipse and circle.
!
! Local variables
double precision :: area_flags_c(2), area_flags_s(2) 
double precision :: area_c, area_s, area_e
double precision :: xy_on_c_mid_tng(2),test_c_mid_tng,xyc_e(2),xyc_s(2)
double precision :: ell_tng_arc_area, t_tng(2)
double precision :: crescent_area,lens_area
double precision :: test_e,test_c,test_1,test_2,test_s,test_m,t_s_e(4),t_s_s(4)
double precision :: xy_int_c(2,4), xy_int_s(2,4), t_c_e(4), t_c_c(4)
double precision :: ell_arc_area, circ_arc_area, ell_arc_area_s, area_t
double precision :: xy_off(2,4), xy_on(2,4), t_on(4), t_off(4),xy_off_c(2,4)
double precision :: t_on_e(4), t_off_e(4),t_on_s(4), t_off_s(4), d_k(3)
double precision :: t_on_c_e(4), t_off_c_e(4),t_on_c_c(4), t_off_c_c(4)
double precision :: t_m, xy_m(2), xy_m_c(2), xy_on_c(2,4)
logical :: lcswitch=.false.
integer :: flag_c, flag_s, k(3)
integer :: n_onside,n_offside,i,i_on_1,i_on_2,j_on,k_in,j_off,i_on,i_off
integer :: n_on_s,n_off_s,n_on_c,n_off_c,i_c_1,i_c_2,i_s_1,i_s_2
integer :: verbose1

! Parameters
double precision, parameter :: t_eps = 1.0d-6 ! Limit for test values

verbose1 = verbose_for_calls(verbose)
spot_limb_eclipse = bad_dble

xyc_s(1:2) = (/ellipse_s(i_ell_x_c), ellipse_s(i_ell_y_c)/)

if (verbose >= v_debug) then
  print *,'spot_limb_eclipse: xy_tng,1 = ',xy_tng(1:2,1)
  print *,'                   xy_tng,2 = ',xy_tng(1:2,2)
endif

area_flags_c = ell_ucircle_overlap(ellipse,verbose1,xy_int_c,t_c_e,t_c_c)
area_c = area_flags_c(1)
flag_c = int(area_flags_c(2))
area_flags_s = ell_ell_overlap(ellipse,ellipse_s,verbose1,xy_int_s,t_s_e,t_s_s)
area_s = area_flags_s(1)
flag_s = int(area_flags_s(2))
if (verbose >= v_debug)  then
  print *,'spot_limb_eclipse: flag_c, flag_s = ',flag_c, flag_s
endif

! Ellipse does not overlap circle - 
if (btest(flag_c, b_ell_no_overlap)) then
  spot_limb_eclipse = 0.0d0
  if (verbose >= v_debug) print *,'spot_limb_eclipse: 0, ',spot_limb_eclipse
  return
endif

! Ellipse completely covers the circle
if (btest(flag_c, b_ell_circle_inside_ellipse)) then
  spot_limb_eclipse = 1.0d0
  if (verbose >= v_debug) print *,'spot_limb_eclipse: 1, ',spot_limb_eclipse
  return
endif


! Value of test_c_mid_tng is used to test on which side of the line-of-tangents
! various points lie.
t_m = atan2(xy_tng(2,1)+xy_tng(2,2),xy_tng(1,1)+xy_tng(1,2))
xy_on_c_mid_tng(1:2) = (/cos(t_m), sin(t_m)/)
test_c_mid_tng = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                               xy_on_c_mid_tng,signed=.true.)

! If ellipse is completely inside the circle and does not overlap the 
! tangent ellipse, we can do a quick check here whether the result is 0.
! This depends on which side of the line between the tangent points ellipse 1
! lies. To test this, calculate the (signed) areas of the triangles formed
! by the tangent points and i. the mid-point between the tangents on the circle
! ii. the centre of the occulting ellipse. If they have different signs then 
! the occulting ellipse cannot cover any part of the required region, so the
! result is 0.
xyc_e = (/ ellipse(i_ell_x_c), ellipse(i_ell_y_c) /)
if (btest(flag_c, b_ell_ellipse_inside_circle).and. &
    btest(flag_s, b_ell_no_overlap) ) then
  test_c = triangle_area(xy_tng(1:2,1), xy_tng(1:2,2), xyc_e,signed=.true.)
  if ((test_c > 0.0d0).neqv.(test_c_mid_tng > 0.0d0)) then
    spot_limb_eclipse = 0.0d0
    if (verbose >= v_debug) print *,'spot_limb_eclipse: 2, ',spot_limb_eclipse
    return
  endif
endif

! If ellipse is completely inside the tangent ellipse and large is false,
! the results is also 0.
if ((.not.large).and.(btest(flag_s,b_ell_1_inside_2))) then
  spot_limb_eclipse = 0.0d0
  if (verbose >= v_debug) print *,'spot_limb_eclipse: 3, ',spot_limb_eclipse
  return
endif

! Area of smaller elliptical arc for the tangent ellipse defined by the tangent 
! points
t_tng(1) = ell_xy_to_t(xy_tng(1:2,1), ellipse_s)
t_tng(2) = ell_xy_to_t(xy_tng(1:2,2), ellipse_s)
ell_tng_arc_area = ell_sector_area(t_tng(1), t_tng(2), ellipse_s) & 
                 - triangle_area(xy_tng(1:2,1), xy_tng(1:2,2), xyc_s)

! Area of the crescent between tangent ellipse and the circle
crescent_area = arc_area(xy_tng(1:2,1),xy_tng(1:2,2)) - ell_tng_arc_area

! Area of the lens-shaped area between tangent ellipse and the circle
lens_area = crescent_area + ellipse_s(i_ell_area)

! If the occulting ellipse is completely inside the circle and does not overlap
! the  tangent ellipse, we already checked for the case where the result is 0,
! so the occulting ellipse must lie in the crescent shaped region between the 
! tangent ellipse and the circle. 
if (btest(flag_c, b_ell_ellipse_inside_circle).and. &
    btest(flag_s, b_ell_no_overlap) ) then
  if (large) then
    spot_limb_eclipse = ellipse(i_ell_area)/lens_area
    if (verbose >= v_debug) print *,'spot_limb_eclipse: 4, ',spot_limb_eclipse
    return
  else
    spot_limb_eclipse = ellipse(i_ell_area)/crescent_area
    if (verbose >= v_debug) print *,'spot_limb_eclipse: 5, ',spot_limb_eclipse
    return
  endif
endif

! If the occulting ellipse is completely inside the tangent ellipse. We
! already dealt with the case large is .false., here is the result if it is 
! .true.. 
if (btest(flag_s,b_ell_1_inside_2)) then
  spot_limb_eclipse = ellipse(i_ell_area)/lens_area
  if (verbose >= v_debug) print *,'spot_limb_eclipse: 6, ',spot_limb_eclipse
  return
endif

! Tangent ellipse is completely inside ellipse 1 
if (btest(flag_s, b_ell_2_inside_1)) then

  ! If there are two intersection points between the circle and the occulting
  ! ellipse then result depends on which side of the line connecting the 
  ! tangent points they lie. Only need to check one of the intersection points 
  ! because we know here that tangent ellipse is completely inside the occulting
  ! ellipse so they must both be on the same side of this line.
  if (btest(flag_c, b_ell_two_intersects)) then
    test_1 = triangle_area(xy_tng(1:2,1), xy_tng(1:2,2), xy_int_c(1:2,1), &
                           signed=.true.)
    if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
      ! Area of the smaller circular arc defined by the intersection points of 
      ! the occulting ellipse and the circle
      circ_arc_area =  arc_area(xy_int_c(1:2,1),xy_int_c(1:2,2))
      ! Area of the smaller elliptical arc on the occulting ellipse defined by 
      ! the intersection points 
      ell_arc_area = ell_sector_area(t_c_e(1),t_c_e(2),ellipse) - &
                     triangle_area(xyc_e, xy_int_c(1:2,1),xy_int_c(1:2,2))
      if (large) then
        spot_limb_eclipse = 1.0d0 - (circ_arc_area-ell_arc_area)/lens_area
        if (verbose >= v_debug) print *,'spot_limb_eclipse: 7, ', &
          spot_limb_eclipse
        return
      else
        spot_limb_eclipse = 1.0d0 - (circ_arc_area-ell_arc_area)/crescent_area
        if (verbose >= v_debug) print *,'spot_limb_eclipse: 8, ', &
          spot_limb_eclipse
        return
      endif

    else

      ! Region of interest is completely covered by the occulting ellipse
      spot_limb_eclipse = 1.0d0
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 9',spot_limb_eclipse
      return

    endif
  endif
  ! Still dealing here with the case of the tangent ellipse completely inside
  ! the occulting ellipse.
  !
  ! We have dealt with the case of two intersections between the occulting
  ! ellipse and the circle, and the case of the circle inside the occulting
  ! ellipse, so the remaining possibility is that there are four 
  ! intersections. 
  !
  ! Count the number of circle/ellipse intersection points that are on the
  ! same  side of the line between the tangent points as the mid-point on the 
  ! circle between the tangent points. Store the first/only two for later use.
  n_onside = 0
  do i=1,4
    test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_int_c(1:2,i),signed=.true.)
    if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
      select case (n_onside)
      case (0)
        n_onside = 1
        i_on_1 = i 
      case (1)
        n_onside = 2
        i_on_2 = i
      case default
        n_onside = n_onside + 1
      end select
    endif
  enddo 

  select case (n_onside)

  case(0)
    spot_limb_eclipse = 1.0d0
    if (verbose >= v_debug) print *,'spot_limb_eclipse: 10',spot_limb_eclipse
    return

  case(2)
    ! Area of elliptical arc on occulting ellipse between the two
    ! intersection points on the same side of line-of-tangents as mid-point
    ! on circle between tangent points.
    ell_arc_area = ell_sector_area(t_c_e(i_on_1),t_c_e(i_on_2),ellipse) - &
              triangle_area(xyc_e, xy_int_c(1:2,i_on_1),xy_int_c(1:2,i_on_2))
    ! Area of circular arc between the same two intersection points.
    circ_arc_area =  arc_area(xy_int_c(1:2,i_on_1),xy_int_c(1:2,i_on_2))
    if (large) then
      spot_limb_eclipse = 1.0d0- (circ_arc_area-ell_arc_area)/lens_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 11', spot_limb_eclipse
      return
    else
      spot_limb_eclipse = 1.0d0- (circ_arc_area-ell_arc_area)/crescent_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 12', spot_limb_eclipse
      return
    endif

  case(4)
    if (large) then
      spot_limb_eclipse = 1.0d0 - (pi-area_c)/lens_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 13', spot_limb_eclipse
      return
    else
      spot_limb_eclipse = 1.0d0 - (pi-area_c)/crescent_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 14', spot_limb_eclipse
      return
    endif

  case default
    if (verbose >= v_debug) then
      print *,'spot_limb_eclipse: n_onside = ',n_onside
      print *,ellipse
      print *,ellipse_s
    endif
    spot_limb_eclipse = bad_dble
    return

  end select

endif

! Occulting ellipse has four intersections with the circle and two with the
! tangent ellipse
if (btest(flag_c, b_ell_four_intersects).and. &
    btest(flag_s, b_ell_two_intersects) ) then

  ! Determine on which side of the line-of-tangents the ellipse-circle
  ! intersections occur. 
  ! n_onside = no. of points on SAME side as mid-point on circle between
  ! ellipse-circle tangent points, positions xy_on(1:2,1:n_onside)
  ! n_offside = no. of points on OPPOSITE side as mid-point on circle between
  ! ellipse-circle tangent points, positions xy_off(1:2,1:n_offside)
  n_onside = 0
  n_offside = 0
  do i=1,4
    test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_int_c(1:2,i),signed=.true.)
    if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
      n_onside =  n_onside + 1
      xy_on(1:2,n_onside) = xy_int_c(1:2,i)
      t_on(n_onside) = t_c_e(i)
    else
      n_offside = n_offside + 1
      xy_off(1:2,n_offside) = xy_int_c(1:2,i)
      t_off(n_offside) = t_c_e(i)
    endif
  enddo 

  select case (n_onside)

  case(0)
    if (large) then
      spot_limb_eclipse = area_s/lens_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 17', spot_limb_eclipse
      return
    else
      spot_limb_eclipse = 0.0d0
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 18', spot_limb_eclipse
      return
    endif

  case(1)
    ! Ellipse-ellipse intersection points lie on different sides of the 
    ! line-of-tangents. Find which one is on the same side as the mid-point on
    ! the circle between the tangent points ("on-side").
    test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_int_s(1:2,1),signed=.true.)
    if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
      j_on =  1
    else
      j_on =  2
    endif
    ! Which tangent point is inside the occulting ellipse?
    if (ell_point_is_inside(xy_tng(1:2,1),ellipse)) then
      k_in = 1
    else
      k_in = 2
    endif
    ! Calculate the area of the three-sided shape defined by the two elliptical
    ! arcs and the circular arc between by the tangent point, the ellipse-
    ! ellipse intersection point and the ellipse-circle intersection point.
    circ_arc_area =  arc_area(xy_on(1:2,1),xy_tng(1:2,k_in))
    ell_arc_area = ell_sector_area(t_on(1),t_s_e(j_on),ellipse) - &
                   triangle_area(xy_on(1:2,1), xy_int_s(1:2,j_on), xyc_e)
    ell_arc_area_s = ell_sector_area(t_tng(k_in),t_s_s(j_on),ellipse_s) - &
                     triangle_area(xy_tng(1:2,k_in), xy_int_s(1:2,j_on), xyc_s)
    area_t = ell_arc_area + circ_arc_area - ell_arc_area_s &
           + triangle_area(xy_on(1:2,1),xy_tng(1:2,k_in),xy_int_s(1:2,j_on))
    if (large) then
      spot_limb_eclipse = (area_s + area_t)/lens_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 19', spot_limb_eclipse
      return
    else
      spot_limb_eclipse = area_t/crescent_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 20', spot_limb_eclipse
      return
    endif

  case(2)
    ! Area of the circular arc defined by the two intersection points of
    ! the occulting ellipse and the circle that are on the same side of the
    ! line-of-tangents as the mid-point on the circle between the tangent points
    !
    circ_arc_area =  arc_area(xy_on(1:2,1),xy_on(1:2,2))
    ! Area of elliptical arc for occulting ellipse defined by the same 
    ! intersection points
    ell_arc_area = ell_sector_area(t_on(1),t_on(2),ellipse) - &
                   triangle_area(xy_on(1:2,1), xy_on(1:2,2), xyc_e)

    ! Two possibilities here, either the ellipse-ellipse intersection points are
    ! on the same side of the line-of-tangents as the mid-point on the circle
    ! between the tangent points, or both are on the opposite side. Only need to
    ! test one of the intersection points.
    test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_int_s(1:2,1),signed=.true.)
    if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then

      if (large) then
        spot_limb_eclipse = 1.0d0 - (circ_arc_area - ell_arc_area)/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 21', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse = 1.0d0 -  &
        (circ_arc_area-ell_arc_area - (ellipse_s(i_ell_area)-area_s)) &
        /crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 22', spot_limb_eclipse
        endif
        return
      endif

    else

      if (large) then
        spot_limb_eclipse = 1.0d0 - &
        (circ_arc_area - ell_arc_area+ellipse_s(i_ell_area)-area_s)/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 23', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse = 1.0d0 -  (circ_arc_area-ell_arc_area) &
        /crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 24', spot_limb_eclipse
        endif
        return
      endif

    endif

  case (3)
    ! Find which ellipse-ellipse intersection point is on the opposite side of 
    ! the line-of-tangents to the mid-point on the circle between the tangent
    ! points.
    test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_int_s(1:2,1),signed=.true.)
    if ((test_c_mid_tng > 0.0d0).neqv.(test_1 > 0.0d0)) then
      j_off =  1
    else
      j_off =  2
    endif
    ! Which tangent point is inside the occulting ellipse?
    if (ell_point_is_inside(xy_tng(1:2,1),ellipse)) then
      k_in = 1
    else
      k_in = 2
    endif
    ! Calculate the area of the three-sided shape defined by the two elliptical
    ! arcs and the circular arc between by the tangent point inside the
    ! occulting ellipse, the ellipse-ellipse intersection point on the
    ! opposite side to the line-of-tangents to the mid-point on the circle 
    ! between the tangent points, and the ellipse-circle intersection point on 
    ! the same side of the line-of-tangents.

   ! Find the mid-point on the circle between the tangent points inside the
   ! occulting ellipse and the ellipse-circle intersection point on the opposite
   ! side to the line-of-tangents to the mid-point on the circle between the
   ! tangent points.
   t_m = midangle(atan2(xy_tng(2,k_in),xy_tng(1,k_in)), &
                  atan2(xy_off(2,1),xy_off(1,1)))
   xy_m_c = (/cos(t_m), sin(t_m)/)
   ! Test on which side of the line-of-tangents this point lies and use this
   ! to determine whether to calculate the larger or smaller of the two possible
   ! circular arc areas defined by the intersection/tangent points.
    test_c = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_m_c,signed=.true.)
    if ((test_c_mid_tng > 0.0d0).eqv.(test_c > 0.0d0)) then
      circ_arc_area =  pi - arc_area(xy_off(1:2,1),xy_tng(1:2,k_in))
      lcswitch = .true.
    else
      circ_arc_area =  arc_area(xy_off(1:2,1),xy_tng(1:2,k_in))
      lcswitch = .false.
    endif

    ! Area of elliptical arc for occulting ellipse. Here again we need to test
    ! the mid-point between the intersection points to determine whether to use
    ! the larger or the smaller area.
    xy_m = ell_t_to_xy(midangle(t_off(1),t_s_e(j_off)),ellipse)
    test_e = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_m,signed=.true.)
    if (((test_e > 0.0d0).eqv.(test_c > 0.0d0)).neqv.lcswitch) then
      ell_arc_area = ell_sector_area(t_off(1),t_s_e(j_off),ellipse) - &
                     triangle_area(xy_off(1:2,1), xy_int_s(1:2,j_off), xyc_e)
    else
      ell_arc_area = ellipse(i_ell_area) &
                   - ell_sector_area(t_off(1),t_s_e(j_off),ellipse) &
                   + triangle_area(xy_off(1:2,1), xy_int_s(1:2,j_off), xyc_e)
    endif

    ! Area of elliptical arc for tangent ellipse - same test with mid-point
    xy_m = ell_t_to_xy(midangle(t_tng(k_in),t_s_s(j_off)),ellipse_s)
    test_s = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_m,signed=.true.)
    if (((test_s > 0.0d0).eqv.(test_c > 0.0d0)).neqv.lcswitch) then
      ell_arc_area_s = ell_sector_area(t_tng(k_in),t_s_s(j_off),ellipse_s)  &
                     - triangle_area(xy_tng(1:2,k_in),xy_int_s(1:2,j_off),xyc_s)
    else
      ell_arc_area_s = ellipse_s(i_ell_area)  &
                     - ell_sector_area(t_tng(k_in),t_s_s(j_off),ellipse_s)  &
                     + triangle_area(xy_tng(1:2,k_in),xy_int_s(1:2,j_off),xyc_s)
    endif
    area_t = ell_arc_area + circ_arc_area - ell_arc_area_s &
           + triangle_area(xy_off(1:2,1),xy_tng(1:2,k_in),xy_int_s(1:2,j_off))
    if (large) then
      spot_limb_eclipse = (area_c - area_t)/lens_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 25', spot_limb_eclipse
      endif
      return
    else
      spot_limb_eclipse = (area_c-area_s-area_t)/crescent_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 26', spot_limb_eclipse
      endif
      return
    endif

  case (4)
  ! Two possible orientations of the occulting ellipse and the tangent ellipse 
  ! here that can be distinguished by whether the point on the circle opposite
  ! the mid-point between the tangents (anti-mid-point?) is inside or outside
  ! the occulting ellipse.
    
  if (ell_point_is_inside(-xy_m_c, ellipse)) then
    if (large) then
      spot_limb_eclipse = 1.0d0  &
      - (pi - area_c - ellipse_s(i_ell_area)- area_s)/lens_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 27', spot_limb_eclipse
      endif
      return
    else
      spot_limb_eclipse = 1.0d0 - (pi - area_c - ellipse_s(i_ell_area)- area_s)&
      / crescent_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 28', spot_limb_eclipse
      endif
      return
    endif
  else
    if (large) then
      spot_limb_eclipse = area_c/lens_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 29', spot_limb_eclipse
      endif
      return
    else
      spot_limb_eclipse = (area_c-area_s)/crescent_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 30', spot_limb_eclipse
      endif
      return
    endif
  endif

  end select
endif

! Occulting ellipse intersects circle at two points but has no overlap with
! tangent  ellipse.
if (btest(flag_c, b_ell_two_intersects).and. &
    btest(flag_s, b_ell_no_overlap) ) then

  ! Test if ellipse-circle intersection points are on the same side of the 
  ! line of tangents as the mid-point on the circle between the tangent points.
  ! (Only need to test one of the points)
  test_1 = triangle_area(xy_tng(1:2,1), xy_tng(1:2,2), xy_int_c(1:2,1), &
                         signed=.true.)
  if ((test_c_mid_tng > 0.0d0).neqv.(test_1 > 0.0d0)) then
    spot_limb_eclipse = 0.0d0
    if (verbose >= v_debug) print *,'spot_limb_eclipse: 31', spot_limb_eclipse
    return
  endif

  if (large) then
    spot_limb_eclipse = area_c/lens_area
    if (verbose >= v_debug) then
      print *,'spot_limb_eclipse: 32', spot_limb_eclipse
    endif
    return
  else
    spot_limb_eclipse = area_c/crescent_area
    if (verbose >= v_debug) then
      print *,'spot_limb_eclipse: 33', spot_limb_eclipse
    endif
    return
  endif
endif


! Occulting ellipse intersects circle at two points and tangent ellipse at two
! points.

if (btest(flag_c, b_ell_two_intersects).and. &
    btest(flag_s, b_ell_two_intersects) ) then

  ! Values to determine on which side of the line-of-tangents each
  ! ellipse-ellipse intersection occurs. 
  test_1 = triangle_area(xy_tng(1:2,1), xy_tng(1:2,2), xy_int_s(1:2,1), &
                         signed=.true.)
  test_2 = triangle_area(xy_tng(1:2,1), xy_tng(1:2,2), xy_int_s(1:2,2), &
                         signed=.true.)
                         
  ! Both ellipse-ellipse intersections on same side of line-of-tangents as 
  ! mid-point on circle between tangent points 
  if ( ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) .and. &
       ((test_c_mid_tng > 0.0d0).eqv.(test_2 > 0.0d0)) ) then

    xy_m = ell_t_to_xy(midangle(t_c_e(1),t_c_e(2)),ellipse)
    if  (hypot(xy_m(1),xy_m(2)) < 1.0d0) then

      if (large) then
        spot_limb_eclipse = area_c/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 34', spot_limb_eclipse
         !!do i_on=1,2
         !!  write (34,*) xy_int_s(1:2,i_on)
         !!end do
         !!write (35,*) xy_tng(1:2,1)
         !!write (35,*) xy_tng(1:2,2)
         !!do i_on=1,10001
         !!  t_m = 2*pi*i_on/10001.0d0
         !!  write (36,*) cos(t_m),sin(t_m),&
         !!    ell_t_to_xy(t_m,ellipse_s), &
         !!    ell_t_to_xy(t_m,ellipse)
         !!end do
        endif
        return
      else
        spot_limb_eclipse = (area_c-area_s)/crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 35', spot_limb_eclipse
        endif
        return
      endif
    else
      if (large) then
        spot_limb_eclipse = 1.0D0 - (ellipse_s(i_ell_area)-area_s) &
                          / lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 36', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse = 1.0d0 - &
           (area_c - area_s) / crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 37', spot_limb_eclipse
        endif
        return
      endif
    endif

  endif

  ! Both ellipse-ellipse intersection points on the opposite side to mid-point
  ! on circle between ellipse-circle tangent point
  if ( ((test_c_mid_tng > 0.0d0).neqv.(test_1 > 0.0d0)) .and. &
       ((test_c_mid_tng > 0.0d0).neqv.(test_2 > 0.0d0)) ) then

    if (ell_point_is_inside(xy_on_c_mid_tng, ellipse)) then
      if (large) then
        spot_limb_eclipse = 1.0D0 - (ellipse_s(i_ell_area)-area_s) &
                          / lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 38', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse = 1.0d0
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 39', spot_limb_eclipse
        endif
        return
      endif
    else
      if (large) then
        spot_limb_eclipse = area_s/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 40', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse = 0.0d0
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 41', spot_limb_eclipse
        endif
        return
      endif
    endif

  endif

  ! Remaining case for two intersections of occulting ellipse with circle and
  ! tangent ellipse is ellipse-ellipse intersections on opposite side of the
  ! line-of-tangents.
  !
    
  ! Which ellipse-circle intersection point is "on-side"?
  test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                         xy_int_c(1:2,1),signed=.true.)
  if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
    i_on =  1
  else
    i_on =  2
  endif

  ! Which ellipse-ellipse intersection point is "on-side"?
  test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                         xy_int_s(1:2,1),signed=.true.)
  if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
    j_on =  1
  else
    j_on =  2
  endif

  ! Which tangent point is inside the occulting ellipse?
  if (ell_point_is_inside(xy_tng(1:2,1),ellipse)) then
    k_in = 1
  else
    k_in = 2
  endif

  circ_arc_area =  arc_area(xy_int_c(1:2,i_on),xy_tng(1:2,k_in))
  ell_arc_area = ell_sector_area(t_c_e(i_on),t_s_e(j_on),ellipse) - &
                 triangle_area(xy_int_c(1:2,i_on), xy_int_s(1:2,j_on), xyc_e)
  ell_arc_area_s = ell_sector_area(t_tng(k_in),t_s_s(j_on),ellipse_s) - &
                   triangle_area(xy_tng(1:2,k_in), xy_int_s(1:2,j_on), xyc_s)
  area_t = ell_arc_area + circ_arc_area - ell_arc_area_s &
    + triangle_area(xy_int_c(1:2,i_on),xy_tng(1:2,k_in),xy_int_s(1:2,j_on))

  if (large) then
    spot_limb_eclipse = (area_s + area_t)/lens_area
    if (verbose >= v_debug) print *,'spot_limb_eclipse: 42', spot_limb_eclipse
    return
  else
    spot_limb_eclipse = area_t/crescent_area
    if (verbose >= v_debug) print *,'spot_limb_eclipse: 43', spot_limb_eclipse
    return
  endif

endif

! Occulting ellipse intersects tangent ellipse at two points and is completely
! inside the circle
if (btest(flag_c, b_ell_ellipse_inside_circle).and. &
    btest(flag_s, b_ell_two_intersects) ) then

  ! Determine to which side of line-of-tangents occulting ellipse lies
  test_e = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                         xyc_e,signed=.true.)
  if ((test_c_mid_tng > 0.0d0).eqv.(test_e > 0.0d0)) then
    if (large) then
      spot_limb_eclipse = ellipse(i_ell_area)/lens_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 44', spot_limb_eclipse
      return
    else
      spot_limb_eclipse = (ellipse(i_ell_area)-area_s)/crescent_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 45', spot_limb_eclipse
      return
    endif
  else
    if (large) then
      spot_limb_eclipse = area_s/lens_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 46', spot_limb_eclipse
      return
    else
      spot_limb_eclipse = 0.0D0
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 47', spot_limb_eclipse
      return
    endif
  endif
endif

! Occulting ellipse intersects the tangent ellipse at four points and is
! completely inside circle
if (btest(flag_c, b_ell_ellipse_inside_circle).and. &
    btest(flag_s, b_ell_four_intersects) ) then
  ! Count the number of circle/ellipse intersection points that are on the
  ! same  side of the line between the tangent points as the mid-point on the 
  ! circle between the tangent points. Store the first/only two for later use.
  n_onside = 0
  do i=1,4
    test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_int_s(1:2,i),signed=.true.)
    if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
      select case (n_onside)
      case (0)
        n_onside = 1
        i_on_1 = i 
      case (1)
        n_onside = 2
        i_on_2 = i
      case default
        n_onside = n_onside + 1
      end select
    endif
  enddo 
  select case (n_onside)

  case (0)
    if (large) then
      spot_limb_eclipse = ellipse(i_ell_area)/lens_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 48', spot_limb_eclipse
      return
    else
      spot_limb_eclipse = (ellipse(i_ell_area) - area_s )/crescent_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 49', spot_limb_eclipse
      return
    endif

  case (2)
    ! Calculate area of the region between the two ellipses on the same side of
    ! the line-of-tangents as the mid-point on the circle between the tangent
    ! points.
    ell_arc_area = ell_sector_area(t_s_e(i_on_1),t_s_e(i_on_2),ellipse) - &
      triangle_area(xy_int_s(1:2,i_on_1), xy_int_s(1:2,i_on_2), xyc_e)
    ell_arc_area_s = ell_sector_area(t_s_s(i_on_1),t_s_s(i_on_2),ellipse_s) - &
      triangle_area(xy_int_s(1:2,i_on_1), xy_int_s(1:2,i_on_2), xyc_s)
    area_t = ell_arc_area - ell_arc_area_s
    if (large) then
      spot_limb_eclipse = (area_s + area_t)/lens_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 50', spot_limb_eclipse
      return
    else
      spot_limb_eclipse =  area_t/crescent_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 51', spot_limb_eclipse
      return
    endif

  case (4)
    if (large) then
      spot_limb_eclipse = area_s/lens_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 52', spot_limb_eclipse
      return
    else
      spot_limb_eclipse =  0.0d0
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 53', spot_limb_eclipse
      return
    endif

  end select

endif

! Occulting ellipse intersects the tangent ellipse at four points and the
! circle at two points.
if (btest(flag_c, b_ell_two_intersects).and. &
    btest(flag_s, b_ell_four_intersects) ) then

  ! Determine on which side of the line-of-tangents the ellipse-ellipse
  ! intersections occur. 
  ! n_onside = no. of points on SAME side as mid-point on circle between
  ! ellipse-circle tangent points, positions xy_on(1:2,1:n_onside)
  ! n_offside = no. of points on OPPOSITE side as mid-point on circle between
  ! ellipse-circle tangent points, positions xy_off(1:2,1:n_offside)
  n_onside = 0
  n_offside = 0
  do i=1,4
    test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_int_s(1:2,i),signed=.true.)
    if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
      n_onside =  n_onside + 1
      xy_on(1:2,n_onside) = xy_int_s(1:2,i)
      t_on_e(n_onside) = t_s_e(i)
      t_on_s(n_onside) = t_s_s(i)
    else
      n_offside = n_offside + 1
      xy_off(1:2,n_offside) = xy_int_s(1:2,i)
      t_off_e(n_offside) = t_s_e(i)
      t_off_s(n_offside) = t_s_s(i)
    endif
  enddo 

  select case (n_offside)

  case (0)
    ! Two possible orientations of the occulting ellipse and the tangent 
    ! ellipse possible here that can be distinguished by whether the point on
    ! the circle opposite the mid-point between the tangents (anti-mid-point?) 
    ! is inside or outside ellipse 1.
    if (ell_point_is_inside(-xy_m_c, ellipse)) then

      if (large) then
        spot_limb_eclipse = 1.0d0 - (pi - area_c)/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 55', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse =  1.0d0 - &
          (pi-area_c-ellipse_s(i_ell_area)+area_s)/crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 56', spot_limb_eclipse
        endif
        return
      endif

    else

      if (large) then
        spot_limb_eclipse = area_c/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 57', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse =  (area_c-area_s)/crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 58', spot_limb_eclipse
        endif
        return
      endif

    endif

  case (1)
    ! Which ellipse-circle intersection point is "off-side"?
    test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_int_c(1:2,1),signed=.true.)
    if ((test_c_mid_tng > 0.0d0).neqv.(test_1 > 0.0d0)) then
      i_off =  1
    else
      i_off =  2
    endif

    ! Which tangent point is inside the occulting ellipse?
    if (ell_point_is_inside(xy_tng(1:2,1),ellipse)) then
      k_in = 1
    else
      k_in = 2
    endif

    ! Find the mid-point on the circle between xy_int_c(1:2,i_off) and 
    ! xy_tng(1:2,k_in) then check whether this is inside or outside the
    ! occulting ellipse.
    t_m = midangle(atan2(xy_tng(2,k_in),xy_tng(1,k_in)), &
                   atan2(xy_int_c(2,i_off),xy_int_c(1,i_off)))
    xy_m_c = (/cos(t_m), sin(t_m)/)
    if (ell_point_is_inside(xy_m_c, ellipse)) then
      circ_arc_area =  arc_area(xy_int_c(1:2,i_off),xy_tng(1:2,k_in))
      lcswitch = .false.
    else
      circ_arc_area =  pi - arc_area(xy_int_c(1:2,i_off),xy_tng(1:2,k_in))
      lcswitch = .true.
    endif
    ! Area of elliptical arc for occulting ellipse depends on whether mid-point
    ! between xy_int_c(1:2,i_off) and xy_off(1:2,1) is on the same side of the
    ! line between these points as the tangent point inside the occulting
    ! ellipse.
    xy_m = ell_t_to_xy(midangle(t_off_e(1), t_c_e(i_off)) ,ellipse)
    test_m = triangle_area(xy_int_c(1:2,i_off), xy_off(1:2,1), &
      xy_m, signed=.true.)
    test_e = triangle_area(xy_int_c(1:2,i_off), xy_off(1:2,1), &
      xy_tng(1:2,k_in), signed=.true.)
    ! These test values may be unreliable if they are very small so..
    if (min(abs(test_m),abs(test_e)) < t_eps) then
      ell_arc_area = 0.0d0
    else
      if ((test_m > 0.0d0).neqv.(test_e > 0.0d0)) then
        ell_arc_area = ell_sector_area(t_off_e(1), t_c_e(i_off) ,ellipse) &
          - triangle_area(xy_off(1:2,1),xy_int_c(1:2,i_off),xyc_e)
      else
        ell_arc_area = ellipse(i_ell_area)   &
          -  ell_sector_area(t_off_e(1), t_c_e(i_off) ,ellipse) &
          + triangle_area(xy_off(1:2,1),xy_int_c(1:2,i_off),xyc_e)
      endif
    endif
    ! Area of elliptical arc for tangent ellipse depends on whether the
    ! mid-point between xy_tng(1:2,k_in) and xy_off(1:2,1) is on the same side
    ! or opposite side of the line between these points as xy_m_c
    xy_m = ell_t_to_xy(midangle(t_tng(k_in), t_off_s(1)), ellipse_s)
    test_m = triangle_area(xy_tng(1:2,k_in), xy_off(1:2,1), xy_m, &
      signed=.true.)
    test_c = triangle_area(xy_tng(1:2,k_in), xy_off(1:2,1), xy_m_c, &
      signed=.true.)
    if (((test_m > 0.0D0).eqv.(test_c > 0.0D0)).neqv.lcswitch) then
      ell_arc_area_s = ell_sector_area(t_tng(k_in), t_off_s(1), ellipse_s)  & 
        - triangle_area(xy_tng(1:2,k_in), xy_off(1:2,1), xyc_s)
    else
      ell_arc_area_s = ellipse_s(i_ell_area) &
        - ell_sector_area(t_tng(k_in), t_off_s(1), ellipse_s)  & 
        + triangle_area(xy_tng(1:2,k_in), xy_off(1:2,1), xyc_s)
    endif
    area_t = ell_arc_area + circ_arc_area - ell_arc_area_s & 
      + triangle_area(xy_tng(1:2,k_in),xy_off(1:2,1),xy_int_c(1:2,i_off))
    if (large) then
      spot_limb_eclipse = (area_c-area_t)/lens_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 59', spot_limb_eclipse
      endif
      return
    else
      spot_limb_eclipse =  (area_c-area_t-area_s)/crescent_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 60', spot_limb_eclipse
      endif
      return
    endif
    
  case (2)
    ! Determine whether ellipse-circle intersection points are on the same side
    ! of the line-of-tangents as the mid-point on the circle between the two
    ! tangent points. (Only need to test one of the points).
    test_1 = triangle_area(xy_tng(1:2,1), xy_tng(1:2,2), xy_int_c(1:2,1), &
                           signed=.true.)
    if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
      ! The ellipse-circle intersection points are on the same side ...
      !  Area of the ellipse arc defined by the two ellipse-ellipse intersection
      ! points on the opposite side of the line of tangents to the midpoint on
      ! the circle between the two tangent points ("offside").
      !  This arc may be either the larger or the smaller of the two arcs
      ! defined by the intersection points - determine which one to use by
      ! testing whether the mid-point between the intersection points is on
      ! the same or the opposite side of the line-of-tangents at the midpoint on
      ! circle between the two tangent points.
      xy_m = ell_t_to_xy(midangle(t_off_e(1), t_off_e(2)), ellipse)
      test_m = triangle_area(xy_tng(1:2,1), xy_tng(1:2,2), xy_m, signed=.true.)
      if ((test_m > 0.0D0).eqv.(test_c_mid_tng > 0.0D0)) then
        ell_arc_area = ellipse(i_ell_area) &
        - ell_sector_area(t_off_e(1), t_off_e(2), ellipse)  & 
        + triangle_area(xy_off(1:2,1), xy_off(1:2,2), xyc_e)
      else
        ell_arc_area = ell_sector_area(t_off_e(1), t_off_e(2), ellipse)  & 
        - triangle_area(xy_off(1:2,1), xy_off(1:2,2), xyc_e)
      endif
      ! Area of the smaller elliptical arc of the tangent ellipse between the
      ! same two ellipse-ellipse intersection points.
      ell_arc_area_s = ell_sector_area(t_off_s(1), t_off_s(2), ellipse_s)  &
        - triangle_area(xy_off(1:2,1), xy_off(1:2,2), xyc_s)
      if (large) then
        spot_limb_eclipse = (area_c - ell_arc_area + ell_arc_area_s)/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 61', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse = (area_c - ell_arc_area + ell_arc_area_s - area_s) &
                            / crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 62', spot_limb_eclipse
        endif
        return
      endif
    else
      ! Area of the ellipse arc defined by the two ellipse-ellipse intersection
      ! points on the same side of the line of tangents to the midpoint on the
      ! circle between the two tangent points. Of the two possible ellipse arcs,
      ! use the one that also lies on the same side of the line of tangents as
      ! the midpoint on the circle between the two tangent points.
      
      xy_m = ell_t_to_xy(midangle(t_on_e(1), t_on_e(2)), ellipse)
      test_m = triangle_area(xy_tng(1:2,1), xy_tng(1:2,2), xy_m, signed=.true.)
      if ((test_m > 0.0D0).eqv.(test_c_mid_tng > 0.0D0)) then
        ell_arc_area = ell_sector_area(t_on_e(1), t_on_e(2), ellipse)  & 
        - triangle_area(xy_on(1:2,1), xy_on(1:2,2), xyc_e)
      else
        ell_arc_area = ellipse(i_ell_area) &
        - ell_sector_area(t_on_e(1), t_on_e(2), ellipse)  & 
        + triangle_area(xy_on(1:2,1), xy_on(1:2,2), xyc_e)
      endif
      ! Area of the smaller elliptical arc of the tangent ellipse between the
      ! same two ellipse-ellipse intersection points.
      ell_arc_area_s = ell_sector_area(t_on_s(1), t_on_s(2), ellipse_s)  &
        - triangle_area(xy_on(1:2,1), xy_on(1:2,2), xyc_s)
      if (large) then
        spot_limb_eclipse = (area_s + ell_arc_area - ell_arc_area_s)/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 63', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse = (ell_arc_area- ell_arc_area_s) / crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 64', spot_limb_eclipse
        endif
        return
      endif
       
    endif

  case (3)
    ! Which ellipse-circle intersection point is "on-side"?
    test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_int_c(1:2,1),signed=.true.)
    if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
      i_on =  1
    else
      i_on =  2
    endif

    ! Which tangent point is inside the occulting ellipse?
    if (ell_point_is_inside(xy_tng(1:2,1),ellipse)) then
      k_in = 1
    else
      k_in = 2
    endif

    circ_arc_area =  arc_area(xy_int_c(1:2,i_on),xy_tng(1:2,k_in))
    ell_arc_area = ell_sector_area(t_c_e(i_on),t_on_e(1),ellipse) - &
                   triangle_area(xy_int_c(1:2,i_on), xy_on(1:2,1), xyc_e)
    ell_arc_area_s = ell_sector_area(t_tng(k_in),t_on_s(1),ellipse_s) - &
                     triangle_area(xy_tng(1:2,k_in), xy_on(1:2,1), xyc_s)
    area_t = ell_arc_area + circ_arc_area - ell_arc_area_s &
      + triangle_area(xy_int_c(1:2,i_on),xy_tng(1:2,k_in),xy_on(1:2,1))

    if (large) then
      spot_limb_eclipse = (area_s + area_t)/lens_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 65', spot_limb_eclipse
      return
    else
      spot_limb_eclipse = area_t/crescent_area
      if (verbose >= v_debug) print *,'spot_limb_eclipse: 66', spot_limb_eclipse
      return
    endif

  case (4)
    if (large) then
      spot_limb_eclipse = area_s/lens_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 67', spot_limb_eclipse
      endif
      return
    else
      spot_limb_eclipse =  0.0d0
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 68 ', spot_limb_eclipse
      endif
      return
    endif

  end select

endif 

! Occulting ellipse intersects the circle at four points but does not overlap
! the tangent ellipse.
if (btest(flag_c, b_ell_four_intersects).and. &
    btest(flag_s, b_ell_no_overlap) ) then
  ! Test on which side of the line-of-tangents the ellipse-circle intersects
  ! lie (only need to test one)
  test_1 = triangle_area(xy_tng(1:2,1), xy_tng(1:2,2), xy_int_c(1:2,1), &
                         signed=.true.)
  if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
    if (large) then
      spot_limb_eclipse = area_c/lens_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 69', spot_limb_eclipse
      endif
      return
    else
      spot_limb_eclipse =  area_c/crescent_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 70 ', spot_limb_eclipse
      endif
      return
    endif
  else 
    spot_limb_eclipse =  0.0d0
    if (verbose >= v_debug) print *,'spot_limb_eclipse: 71 ', spot_limb_eclipse
  endif

endif

if (btest(flag_c, b_ell_four_intersects).and. &
    btest(flag_s, b_ell_four_intersects) ) then
    
  ! Determine on which side of the line-of-tangents the ellipse-ellipse
  ! intersections occur. 
  ! n_onside = no. of points on SAME side as mid-point on circle between
  ! ellipse-circle tangent points, positions xy_on(1:2,1:n_onside)
  ! n_offside = no. of points on OPPOSITE side as mid-point on circle between
  ! ellipse-circle tangent points, positions xy_off(1:2,1:n_offside)
  n_on_s = 0
  n_off_s = 0
  do i=1,4
    test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), &
                           xy_int_s(1:2,i),signed=.true.)
    if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
      n_on_s =  n_on_s + 1
      xy_on(1:2,n_on_s) = xy_int_s(1:2,i)
      t_on_e(n_on_s) = t_s_e(i)
      t_on_s(n_on_s) = t_s_s(i)
    else
      n_off_s = n_off_s + 1
      xy_off(1:2,n_off_s) = xy_int_s(1:2,i)
      t_off_e(n_off_s) = t_s_e(i)
      t_off_s(n_off_s) = t_s_s(i)
    endif
  enddo 
  
  select case(n_off_s)

  case(0)
    ! Two cases here depending on whether there are 2 or 4 ellipse-circle
    ! intersection points on the same side of the line-of-tangents as the
    ! mid-point on the circle between the tangent points.
    n_on_c = 0
    n_off_c = 0
    do i=1,4
      test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), xy_int_c(1:2,i),&
                             signed=.true.)
      if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
        n_on_c = n_on_c + 1
      else
        n_off_c = n_off_c + 1
        xy_off_c(1:2,n_off_c) = xy_int_c(1:2,i)
        t_off_c_e(n_off_c) = t_c_e(i)
        t_off_c_c(n_off_c) = t_c_c(i)
      endif
    end do

    select case (n_on_c)
    case (2)
      if (large) then
        ! Area between circle and ellipse opposite the mid-point on the circle
        ! between the line of tangents depends on whether the mid-point on the
        ! circle between the "off-side" ellipse-circle intersection points is
        ! on the same side of the line between these points as the mid-point
        ! on the circle between the tangent points.
        t_m = midangle(t_off_c_c(1), t_off_c_c(2))
        xy_m = (/cos(t_m), sin(t_m)/)
        test_c = triangle_area(xy_off_c(1:2,1), xy_off_c(1:2,2), xy_m, &
                               signed=.true.)
        test_m = triangle_area(xy_off_c(1:2,1), xy_off_c(1:2,2), xy_m_c, &
                               signed=.true.)
        if ((test_c > 0.0d0).eqv.(test_m > 0.0d0)) then
          circ_arc_area = pi - arc_area(xy_off_c(1:2,1), xy_off_c(1:2,2))
        else 
          circ_arc_area = arc_area(xy_off_c(1:2,1), xy_off_c(1:2,2))
        endif

        ell_arc_area = ell_sector_area(t_off_c_e(1),t_off_c_e(2),ellipse) - &
                 triangle_area(xy_off_c(1:2,1), xy_off_c(1:2,2), xyc_e)
        spot_limb_eclipse =  1.0d0 - &
          (pi-area_c-circ_arc_area+ell_arc_area-ellipse_s(i_ell_area)+area_s) &
          / lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 72 ', spot_limb_eclipse
        endif
        return
        
      else
        spot_limb_eclipse =  (ellipse_s(i_ell_area) - area_s)/crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 73 ', spot_limb_eclipse
        endif
      endif

      case (4)
      if (large) then
        spot_limb_eclipse = area_s/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 74', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse =  (area_c-area_s)/crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 75 ', spot_limb_eclipse
        endif
        return
      endif

    end select

  case (1)
    ! Two possibilities here depending on whether there are 1 or 3
    ! ellipse-circle intersection points on the same side of the
    ! line-of-tangents as the mid-point on circle between tangent points.
    n_on_c = 0
    n_off_c = 0
    do i=1,4
      test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), xy_int_c(1:2,i),&
                             signed=.true.)
      if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
        n_on_c = n_on_c + 1
        xy_on_c(1:2,n_on_c) = xy_int_c(1:2,i)
        t_on_c_e(n_on_c) = t_c_e(i)
        t_on_c_c(n_on_c) = t_c_c(i)
      else
        n_off_c = n_off_c + 1
        xy_off_c(1:2,n_off_c) = xy_int_c(1:2,i)
        t_off_c_e(n_off_c) = t_c_e(i)
        t_off_c_c(n_off_c) = t_c_c(i)
      endif
    end do

    ! Which tangent point is inside the occulting ellipse?
    if (ell_point_is_inside(xy_tng(1:2,1),ellipse)) then
      k_in = 1
    else
      k_in = 2
    endif

    select case (n_on_c)

    case (1) 
      ! Find which of the three ellipse-ellipse intersection points on the same
      ! side of the line-of-tangents as the mid-point on circle between the
      ! tangent point is closest to the ellipse-circle intersection point on the
      ! same side of the line-of-tangents ("on-side").
      do i=1,3
        d_k(i) = distance(xy_on(1:2,i), xy_on_c(1:2,1))
      end do
      call heapsort(d_k,k)
      ! Calculate the area of the three-sided shape defined by the two
      ! elliptical arcs and the circular arc between by the tangent point inside
      ! the occulting ellipse, the on-side ellipse-circle intersection point and
      ! the nearest on-side ellipse-ellipse intersection point to this
      ! ellipse-circle intersection point.
      
      ! Area of the circular arc
      circ_arc_area = arc_area(xy_on_c(1:2,1), xy_tng(1:2,k_in))

      ! Area of elliptical arc for occulting ellipse
      ell_arc_area = ell_sector_area(t_on_c_e(1),t_on_e(k(1)),ellipse) - &
        triangle_area(xy_on_c(1:2,1), xy_on(1:2,k(1)), xyc_e)
      
      ! Area of elliptical arc for tangent ellipse
      ell_arc_area_s = ell_sector_area(t_tng(k_in),t_on_s(k(1)),ellipse_s) - &
        triangle_area(xy_tng(1:2, k_in), xy_on(1:2,k(1)), xyc_s)

      area_t = ell_arc_area + circ_arc_area - ell_arc_area_s  &
      + triangle_area(xy_on_c(1:2,1), xy_tng(1:2,k_in), xy_on(1:2,k(1)))

      ! Area defined by the two elliptical arcs between their two on-side 
      ! intersections that are not closest to the interior tangent point.
      ell_arc_area = ell_sector_area(t_on_e(k(2)),t_on_e(k(3)),ellipse) - &
        triangle_area(xy_on(1:2,k(2)), xy_on(1:2,k(3)), xyc_e)
      ell_arc_area_s = ell_sector_area(t_on_s(k(2)),t_on_s(k(3)),ellipse_s) - &
        triangle_area(xy_on(1:2,k(2)), xy_on(1:2,k(3)), xyc_s)
      area_e = ell_arc_area - ell_arc_area_s

      if (large) then
        spot_limb_eclipse = (ellipse_s(i_ell_area)+area_t+area_e)/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 76', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse =  (area_t+area_e)/crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 77 ', spot_limb_eclipse
        endif
        return
      endif
      
    case (3)
      ! Calculate the area of the three-sided shape defined by the two
      ! elliptical arcs and the circular arc between the tangent point, the
      ! ellipse-ellipse intersection point and the ellipse-circle intersection
      ! point.
      
      ! Find the mid-point on the circle between offside ellipse-circle
      ! intersection and the interior tangent point and check on which side of
      ! the line-of-tangents it lies.
      t_m = midangle(atan2(xy_tng(2,k_in),xy_tng(1,k_in)), t_off_c_c(1))
      xy_m = (/cos(t_m), sin(t_m)/)
      test_c = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), xy_m, &
                             signed=.true.)
      if ((test_c > 0.0d0).eqv.(test_c_mid_tng > 0.0d0)) then
        circ_arc_area = pi - arc_area(xy_off_c(1:2,1), xy_tng(1:2,k_in))
      else 
        circ_arc_area = arc_area(xy_off_c(1:2,1), xy_tng(1:2,k_in))
      endif

      ! Area of elliptical arc for occulting ellipse
      ell_arc_area = ell_sector_area(t_off_c_e(1),t_off_e(1),ellipse) - &
        triangle_area(xy_off_c(1:2,1), xy_off(1:2,1), xyc_e)

      ! Area of elliptical arc for tangent ellipse depends on whether the
      ! mid-point between xy_tng(1:2,k_in) and xy_off(1:2,1) is on the same side
      ! or opposite side of the line between these points as xy_m_c
      xy_m = ell_t_to_xy(midangle(t_tng(k_in), t_off_s(1)), ellipse_s)
      test_m = triangle_area(xy_tng(1:2,k_in), xy_off(1:2,1), xy_m, &
        signed=.true.)
      test_c = triangle_area(xy_tng(1:2,k_in), xy_off(1:2,1), xy_m_c, &
        signed=.true.)
      if (((test_m > 0.0D0).eqv.(test_c > 0.0D0)).neqv.lcswitch) then
        ell_arc_area_s = ell_sector_area(t_tng(k_in), t_off_s(1), ellipse_s)  & 
          - triangle_area(xy_tng(1:2,k_in), xy_off(1:2,1), xyc_s)
      else
        ell_arc_area_s = ellipse_s(i_ell_area) &
          - ell_sector_area(t_tng(k_in), t_off_s(1), ellipse_s)  & 
          + triangle_area(xy_tng(1:2,k_in), xy_off(1:2,1), xyc_s)
      endif

      area_t = ell_arc_area + circ_arc_area - ell_arc_area_s + &
               triangle_area(xy_off_c(1:2,1), xy_off(1:2,1),xy_tng(1:2,k_in))

      if (large) then
        spot_limb_eclipse = (area_c - area_t)/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 78', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse =  (area_c-area_t - area_s)/crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 79 ', spot_limb_eclipse
        endif
        return
      endif

    end select
    
  case (2) !  = n_off_s, number of off-side ellipse-ellipse intersections
    ! Still in 4 ellipse-ellipse and 4 ellipse-circle intersections case
    ! Four possible configurations here distinguished by whether the tangent
    ! points are inside the occulting ellipse or not - only need to check one of
    ! the tangent points - and the number of circle-ellipse intersection points
    ! on the same side of the line of tangents as the mid-point on the circle
    ! between the tangent points ("on-side"). 
    n_on_c = 0
    do i=1,4
      test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), xy_int_c(1:2,i),&
                             signed=.true.)
      if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
        n_on_c = n_on_c + 1
        xy_on_c(1:2,n_on_c) = xy_int_c(1:2,i)
        t_on_c_e(n_on_c) = t_c_e(i)
        t_on_c_c(n_on_c) = t_c_c(i)
        if (n_on_c == 2) exit
      endif
    end do
    if (ell_point_is_inside(xy_tng(1:2,1),ellipse)) then
        ! Calculate the sum of the areas of the two three-sided shapes defined
        ! by the two elliptical arcs and a circular arc between each tangent
        ! point and the nearest on-side ellipse-ellipse intersection point and
        ! the nearest on-side ellipse-circle intersection point.
        
        if (distance(xy_tng(1:2,1), xy_on_c(1:2,1)) < & 
            distance(xy_tng(1:2,1), xy_on_c(1:2,2)) ) then
          i_c_1 = 1
          i_c_2 = 2
        else
          i_c_1 = 2
          i_c_2 = 1
        endif
        
        if (distance(xy_tng(1:2,1), xy_on(1:2,1)) < & 
            distance(xy_tng(1:2,1), xy_on(1:2,2)) ) then
          i_s_1 = 1
          i_s_2 = 2
        else
          i_s_1 = 2
          i_s_2 = 1
        endif
        ! First tangent point
        circ_arc_area = arc_area(xy_on_c(1:2,i_c_1), xy_tng(1:2,1))
        ell_arc_area = &
          ell_sector_area(t_on_c_e(i_c_1),t_on_e(i_s_1),ellipse) - &
          triangle_area(xy_on_c(1:2,i_c_1), xy_on(1:2,i_s_1), xyc_e)
        ell_arc_area_s = ell_sector_area(t_on_s(i_s_1),t_tng(1),ellipse_s) - &
          triangle_area(xy_on(1:2,i_s_1), xy_tng(1:2,1), xyc_s)
        area_t =  circ_arc_area + ell_arc_area - ell_arc_area_s + &
          triangle_area(xy_on_c(1:2,i_c_1), xy_tng(1:2,1),xy_on(1:2,i_s_1))

        ! Second tangent point
        circ_arc_area = arc_area(xy_on_c(1:2,i_c_2), xy_tng(1:2,2))
        ell_arc_area = &
          ell_sector_area(t_on_c_e(i_c_2),t_on_e(i_s_2),ellipse) - &
          triangle_area(xy_on_c(1:2,i_c_2), xy_on(1:2,i_s_2), xyc_e)
        ell_arc_area_s = ell_sector_area(t_on_s(i_s_2),t_tng(2),ellipse_s) - &
          triangle_area(xy_on(1:2,i_s_2), xy_tng(1:2,2), xyc_s)
        area_t =  area_t + circ_arc_area + ell_arc_area - ell_arc_area_s + &
          triangle_area(xy_on_c(1:2,i_c_2), xy_tng(1:2,2),xy_on(1:2,i_s_2))

        if (large) then
          spot_limb_eclipse = (area_s +  area_t)/lens_area
          if (verbose >= v_debug) then
            print *,'spot_limb_eclipse: 80', spot_limb_eclipse
          endif
          return
        else
          spot_limb_eclipse =  area_t/crescent_area
          if (verbose >= v_debug) then
            print *,'spot_limb_eclipse: 81', spot_limb_eclipse
          endif
          return
        endif

      return

    else ! tangent points are not inside the occulting ellipse

      n_on_c = 0
      n_off_c = 0
      do i=1,4
        test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), xy_int_c(1:2,i),&
                               signed=.true.)
        if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
          n_on_c = n_on_c + 1
          xy_on_c(1:2,n_on_c) = xy_int_c(1:2,i)
          t_on_c_e(n_on_c) = t_c_e(i)
          t_on_c_c(n_on_c) = t_c_c(i)
        else
          n_off_c = n_off_c + 1
          xy_off_c(1:2,n_off_c) = xy_int_c(1:2,i)
          t_off_c_e(n_off_c) = t_c_e(i)
          t_off_c_c(n_off_c) = t_c_c(i)
        endif
      end do
      select case (n_on_c)
      case (0) 
        ! Areas of elliptical arcs defined by the ellipse-ellipse intersection
        ! points on the same side of the line of tangents as the mid-point on
        ! the circle between the tangent points ("on-side").

        ! For the occulting ellipse, need to check whether the mid-point on the
        ! ellipse between the on-side intersection points is also on-side
        xy_m = ell_t_to_xy(midangle(t_on_e(1), t_on_e(2)), ellipse)
        test_m = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), xy_m,&
                               signed=.true.)
        if ((test_m > 0.0D0).eqv.(test_c_mid_tng > 0.0D0)) then
          ell_arc_area = ell_sector_area(t_on_e(1), t_on_e(2), ellipse)  & 
          - triangle_area(xy_on(1:2,1), xy_on(1:2,2), xyc_e)
        else
          ell_arc_area = ellipse(i_ell_area) &
          - ell_sector_area(t_on_e(1), t_on_e(2), ellipse)  & 
          + triangle_area(xy_on(1:2,1), xy_on(1:2,2), xyc_e)
        endif
        ! Area of the smaller elliptical arc of the tangent ellipse between the
        ! same two ellipse-ellipse intersection points.
        ell_arc_area_s = ell_sector_area(t_on_s(1), t_on_s(2), ellipse_s)  & 
          - triangle_area(xy_on(1:2,1), xy_on(1:2,2), xyc_s)
        if (large) then
          spot_limb_eclipse = (area_s - ell_arc_area - ell_arc_area_s)/lens_area
          if (verbose >= v_debug) then
            print *,'spot_limb_eclipse: 82', spot_limb_eclipse
          endif
          return
        else
          spot_limb_eclipse =  (ell_arc_area - ell_arc_area_s)/crescent_area
          if (verbose >= v_debug) then
            print *,'spot_limb_eclipse: 83 ', spot_limb_eclipse
          endif
          return
        endif

      case (2) 
        ! Calculate the sum of the areas of the two three-sided shapes defined
        ! by the two elliptical arcs and a circular arc between each tangent
        ! point and the nearest on-side ellipse-ellipse intersection point and
        ! the nearest on-side ellipse-circle intersection point.
        
        if (distance(xy_tng(1:2,1), xy_on_c(1:2,1)) < & 
            distance(xy_tng(1:2,1), xy_on_c(1:2,2)) ) then
          i_c_1 = 1
          i_c_2 = 2
        else
          i_c_1 = 2
          i_c_2 = 1
        endif
        
        if (distance(xy_tng(1:2,1), xy_on(1:2,1)) < & 
            distance(xy_tng(1:2,1), xy_on(1:2,2)) ) then
          i_s_1 = 1
          i_s_2 = 2
        else
          i_s_1 = 2
          i_s_2 = 1
        endif
        ! First tangent point
        circ_arc_area = arc_area(xy_on_c(1:2,i_c_1), xy_tng(1:2,1))
        ell_arc_area = &
          ell_sector_area(t_on_c_e(i_c_1),t_on_e(i_s_1),ellipse) - &
          triangle_area(xy_on_c(1:2,i_c_1), xy_on(1:2,i_s_1), xyc_e)
        ell_arc_area_s = ell_sector_area(t_on_s(i_s_1),t_tng(1),ellipse_s) - &
          triangle_area(xy_on(1:2,i_s_1), xy_tng(1:2,1), xyc_s)
        area_t =  circ_arc_area - ell_arc_area - ell_arc_area_s + &
          triangle_area(xy_on_c(1:2,i_c_1), xy_tng(1:2,1),xy_on(1:2,i_s_1))

        ! Second tangent point
        circ_arc_area = arc_area(xy_on_c(1:2,i_c_2), xy_tng(1:2,2))
        ell_arc_area = &
          ell_sector_area(t_on_c_e(i_c_2),t_on_e(i_s_2),ellipse) - &
          triangle_area(xy_on_c(1:2,i_c_2), xy_on(1:2,i_s_2), xyc_e)
        ell_arc_area_s = ell_sector_area(t_on_s(i_s_2),t_tng(2),ellipse_s) - &
          triangle_area(xy_on(1:2,i_s_2), xy_tng(1:2,2), xyc_s)
        area_t =  area_t + circ_arc_area - ell_arc_area - ell_arc_area_s + &
          triangle_area(xy_on_c(1:2,i_c_2), xy_tng(1:2,2),xy_on(1:2,i_s_2))

        if (large) then
          spot_limb_eclipse = (area_s + crescent_area - area_t)/lens_area
          if (verbose >= v_debug) then
            print *,'spot_limb_eclipse: 84', spot_limb_eclipse
          endif
          return
        else
          spot_limb_eclipse =  1.0d0 - area_t/crescent_area
          if (verbose >= v_debug) then
            print *,'spot_limb_eclipse: 85 ', spot_limb_eclipse
          endif
          return
        endif

      case (4) 
        ! Areas of elliptical arcs defined by the ellipse-ellipse intersection
        ! points on the same side of the line of tangents as the mid-point on
        ! the circle between the tangent points ("off-side").

        ! For the occulting ellipse, need to check whether the mid-point on the
        ! ellipse between the off-side intersection points is also off-side
        xy_m = ell_t_to_xy(midangle(t_off_e(1), t_off_e(2)), ellipse)
        test_m = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), xy_m,&
                               signed=.true.)
        if ((test_m > 0.0D0).eqv.(test_c_mid_tng > 0.0D0)) then
          ell_arc_area = ellipse(i_ell_area) &
          - ell_sector_area(t_off_e(1), t_off_e(2), ellipse)  & 
          + triangle_area(xy_off(1:2,1), xy_off(1:2,2), xyc_e)
        else
          ell_arc_area = ell_sector_area(t_off_e(1), t_off_e(2), ellipse)  & 
          - triangle_area(xy_off(1:2,1), xy_off(1:2,2), xyc_e)
        endif
        ! Area of the smaller elliptical arc of the tangent ellipse between the
        ! same two ellipse-ellipse intersection points.
        ell_arc_area_s = ell_sector_area(t_off_s(1), t_off_s(2), ellipse_s)  & 
          - triangle_area(xy_off(1:2,1), xy_off(1:2,2), xyc_s)
        if (large) then
          spot_limb_eclipse = (area_c - ell_arc_area + ell_arc_area_s)/lens_area
          if (verbose >= v_debug) then
            print *,'spot_limb_eclipse: 86', spot_limb_eclipse
          endif
          return
        else
          spot_limb_eclipse = (area_c - ell_arc_area + ell_arc_area_s - area_s)&
            /crescent_area
          if (verbose >= v_debug) then
            print *,'spot_limb_eclipse: 87 ', spot_limb_eclipse
          endif
          return
        endif

      end select ! n_on_c

    endif

  case (3) !  = n_off_s, number of off-side ellipse-ellipse intersections
    ! Still in 4 ellipse-ellipse and 4 ellipse-circle intersections case
    !
    ! Two possibilities here depending on whether there are 1 or 3
    ! ellipse-circle intersection points on the same side of the
    ! line-of-tangents as the mid-point on circle between tangent points
    ! ("on-side")
    n_on_c = 0
    n_off_c = 0
    do i=1,4
      test_1 = triangle_area(xy_tng(1:2,1),xy_tng(1:2,2), xy_int_c(1:2,i),&
                             signed=.true.)
      if ((test_c_mid_tng > 0.0d0).eqv.(test_1 > 0.0d0)) then
        n_on_c = n_on_c + 1
        xy_on_c(1:2,n_on_c) = xy_int_c(1:2,i)
        t_on_c_e(n_on_c) = t_c_e(i)
        t_on_c_c(n_on_c) = t_c_c(i)
      else
        n_off_c = n_off_c + 1
        xy_off_c(1:2,n_off_c) = xy_int_c(1:2,i)
        t_off_c_e(n_off_c) = t_c_e(i)
        t_off_c_c(n_off_c) = t_c_c(i)
      endif
    end do

    ! Which tangent point is inside the occulting ellipse?
    if (ell_point_is_inside(xy_tng(1:2,1),ellipse)) then
      k_in = 1
    else
      k_in = 2
    endif

    select case (n_off_c)

    case (1) 
      ! Find which of the three ellipse-ellipse intersection points on the
      ! off-side is closest to the off-side ellipse-circle intersection point
      do i=1,3
        d_k(i) = distance(xy_off(1:2,i), xy_off_c(1:2,1))
      end do
      call heapsort(d_k,k)
      ! Calculate the area of the three-sided shape defined by the two
      ! elliptical arcs and the circular arc between by the tangent point inside
      ! the occulting ellipse, the off-side ellipse-circle intersection point
      ! and the nearest off-side ellipse-ellipse intersection point to this
      ! ellipse-circle intersection point.
      
      ! Area of the circular arc
      circ_arc_area = arc_area(xy_off_c(1:2,1), xy_tng(1:2,k_in))

      ! Area of elliptical arc for occulting ellipse
      ell_arc_area = ell_sector_area(t_off_c_e(1),t_off_e(k(1)),ellipse) - &
        triangle_area(xy_off_c(1:2,1), xy_off(1:2,k(1)), xyc_e)
      
      ! Area of elliptical arc for tangent ellipse
      ell_arc_area_s = ell_sector_area(t_tng(k_in),t_off_s(k(1)),ellipse_s) - &
        triangle_area(xy_tng(1:2, k_in), xy_off(1:2,k(1)), xyc_s)

      area_t = ell_arc_area + circ_arc_area - ell_arc_area_s  &
      + triangle_area(xy_off_c(1:2,1), xy_tng(1:2,k_in), xy_off(1:2,k(1)))

      ! Area defined by the two elliptical arcs between their two off-side 
      ! intersections that are not closest to the interior tangent point.
      ell_arc_area = ell_sector_area(t_off_e(k(2)),t_off_e(k(3)),ellipse) - &
        triangle_area(xy_off(1:2,k(2)), xy_off(1:2,k(3)), xyc_e)
      ell_arc_area_s = ell_sector_area(t_off_s(k(2)),t_off_s(k(3)),ellipse_s)- &
        triangle_area(xy_off(1:2,k(2)), xy_off(1:2,k(3)), xyc_s)
      area_e = ell_arc_area - ell_arc_area_s

      if (large) then
        spot_limb_eclipse = (area_c - area_t -area_e)/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 88', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse =  (area_c - area_t -area_e-area_s)/crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 89 ', spot_limb_eclipse
        endif
        return
      endif
      
    case (3) ! n_off_c = no. of offside ellipse, circle intersects
      ! Calculate the area of the three-sided shape defined by the two
      ! elliptical arcs and the circular arc between the tangent point, the
      ! ellipse-ellipse intersection point and the ellipse-circle intersection
      ! point.
      
      circ_arc_area = arc_area(xy_on_c(1:2,1), xy_tng(1:2,k_in))
      ell_arc_area = ell_sector_area(t_on_c_e(1),t_on_e(1),ellipse) - &
        triangle_area(xy_on_c(1:2,1), xy_on(1:2,1), xyc_e)
      ell_arc_area_s = ell_sector_area(t_tng(k_in), t_on_s(1), ellipse_s)  & 
        - triangle_area(xy_tng(1:2,k_in), xy_on(1:2,1), xyc_s)

      area_t = ell_arc_area + circ_arc_area - ell_arc_area_s + &
               triangle_area(xy_on_c(1:2,1), xy_on(1:2,1),xy_tng(1:2,k_in))

      if (large) then
        spot_limb_eclipse = (area_s + area_t)/lens_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 90', spot_limb_eclipse
        endif
        return
      else
        spot_limb_eclipse =  area_t/crescent_area
        if (verbose >= v_debug) then
          print *,'spot_limb_eclipse: 91 ', spot_limb_eclipse
        endif
        return
      endif

    end select
    
  case (4) !  = n_off_s, number of off-side ellipse-ellipse intersections
    ! Still in 4 ellipse-ellipse and 4 ellipse-circle intersections case

    if (large) then
      spot_limb_eclipse = area_s/lens_area
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 92', spot_limb_eclipse
      endif
      return
    else
      spot_limb_eclipse = 0.0d0
      if (verbose >= v_debug) then
        print *,'spot_limb_eclipse: 93 ', spot_limb_eclipse
      endif
      return
    endif
  end select

  endif

end function spot_limb_eclipse

!------------------------------------------------------------------------------

subroutine project_spot(alpha_s, beta_s, gamma_s, ellipse_s, uv_tng)
!
! Ellipse that is the projection of a circular spot on a spherical star onto
! the plane of the sky, and the coordinates of the tangent points of this
! ellipse with the projected limb of the star for spots on the limb.
!
! The spot is assumed to be on a spherical star of unit radius centred on the
! origin. The spot edge is defined by the intersection of a cone with its vertex
! at the origin and with opening angle gamma_s.
!
! The spot longitude alpha_s is measured from the u-axis towards the v-axis in
! the u-v plane.  The star is viewed from the positive w direction. The spot
! latitude beta_s is measured from the u-v plane towards the w-axis
!
implicit none
double precision, intent(in) :: alpha_s, beta_s, gamma_s ! Spot coords/size
!f2py integer, parameter :: n_ell_par = 14
double precision, intent(out) :: ellipse_s(n_ell_par) ! Projected ellipse
double precision, intent(out) :: uv_tng(2,2)          ! Tangent points

! Local variables
double precision :: sina, cosa, sinb, cosb, sing, cosg
double precision :: a_p, b_p, u_c, v_c, sint_t, cost_t, d
double precision :: ellpar_s(5) ! Ellipse parameters for projected spot

! Projection of the spot onto the u-v plane is an ellipse
sina = sin(alpha_s)
sinb = sin(beta_s)
sing = sin(gamma_s)
cosa = cos(alpha_s)
cosb = cos(beta_s)
cosg = cos(gamma_s)
! ellpar_s = [a_p, b_p, u_c, v_c, phi]
a_p = sing
b_p = abs(sinb)*sing
u_c = cosa*cosb*cosg
v_c = sina*cosb*cosg
ellpar_s = (/a_p, b_p, u_c, v_c, halfpi + alpha_s/)
ellipse_s = ell_init_from_par(ellpar_s)

! For a spot on the limb there are two tangent points between the projected
! ellipse and the unit circle corresponding to z=0. Calculate their positions.
d = cosb*sing
if ( d == 0.0d0) then
  uv_tng(:,:) = bad_dble
else
  sint_t = -abs(sinb*cosg/d)
  if (abs(sint_t) <= 1.0d0) then
    cost_t = sqrt(1.0d0 - sint_t**2)
    uv_tng(1,1) = -a_p*sina*cost_t  - b_p*cosa*sint_t + u_c
    uv_tng(2,1) =  a_p*cosa*cost_t  - b_p*sina*sint_t + v_c
    cost_t = -cost_t
    uv_tng(1,2) = -a_p*sina*cost_t  - b_p*cosa*sint_t + u_c
    uv_tng(2,2) =  a_p*cosa*cost_t  - b_p*sina*sint_t + v_c
  else
    uv_tng(:,:) = not_set_dble
  endif
endif


end subroutine project_spot

end module spots
