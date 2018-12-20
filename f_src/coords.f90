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

module coords
implicit none

!
! Coordinate transformation for ellc 
!
! Assumes that the stellar rotation axis and orbital angular momentum axes are
! parallel.
! 
! The star is defined by an ellipsoid which is offset by a distance D from the
! centre-of-mass of the star towards the companion.
!
! Star coordinate system is (x, y, z) with origin at the star's centre-of-mass
! - x towards the companion
! - y in the orbital plane with star motion towards negative y
! - z parallel to rotation/orbital angular momentum vector
!
! Ellipsoid coordinate system is (x', y', z') with origin at the centre of the
! ellipsoid, i.e., x' = x - D, y' = y, z' = z.
!
! Observer coordinate system (u, v, w) with origin at the centre-of-mass of the
! binary
! - u in the plane of the sky perpendicular to the angular momentum vector
! - v in the plane of the sky parallel to the angular momentum vector
! - w pointing towards the observer
!
! The ellipsoid defining the star projected onto the plane of the sky is an
! ellipse. The projected coordinate system is (s, t, p) with the same origin as
! (x',y',z') and the same orientation as (u, v, w).
! 
! The inclination, incl, is the angle between the orbital angular momentum
! vector and the line of sight.
!
! The angle phi is measured from the x-axis to the projection of the w-axis
! in the x-y plane.
!
!  The distance from the centre-of-mass of the binary to the centre-of-mass for 
! star 1 is r_1, and similarly for star_2, r_2.
!
! HISTORY
! -------
!
! Sep 2016
! First version 
! p.maxted@keele.ac.uk
!
!

! *** All angles in radians ***
public xyz2stp  ! Star to projected
public stp2xyz  ! Projected to star
public uvw2stp  ! Observer to projected

contains

function xyz2stp(xyz,D,incl,phi) result (stp)
implicit none

double precision,intent(in) :: xyz(3),D,incl,phi
double precision            :: stp(3)

! Local variables
double precision :: xp, yp, zp, sini, cosi, sinphi, cosphi

sini = sin(incl)
cosi = cos(incl)
sinphi = sin(phi)
cosphi = cos(phi)
xp = xyz(1) - D
yp = xyz(2)
zp = xyz(3)
stp(1) =  xp*sinphi      + yp*cosphi
stp(2) = -xp*cosphi*cosi + yp*sinphi*cosi + zp*sini
stp(3) =  xp*cosphi*sini - yp*sinphi*sini + zp*cosi

end function xyz2stp

!--------------------------------------------------------------------------

function stp2xyz(stp,D,incl,phi) result (xyz)
implicit none

double precision,intent(in) :: stp(3),D,incl,phi
double precision            :: xyz(3)

! Local variables
double precision :: s, t, p, sini, cosi, sinphi, cosphi

sini = sin(incl)
cosi = cos(incl)
sinphi = sin(phi)
cosphi = cos(phi)
s = stp(1)
t = stp(2)
p = stp(3)
xyz(1) =  s*sinphi - t*cosphi*cosi + p*cosphi*sini + D
xyz(2) =  s*cosphi + t*sinphi*cosi - p*sinphi*sini
xyz(3) =             t*sini        + p*cosi

end function stp2xyz

!--------------------------------------------------------------------------

function uvw2stp(uvw,D,incl,phi,r_i) result (stp)
implicit none

double precision,intent(in) :: uvw(3),D,incl,phi, r_i
! r_i is distance from binary c-of-m to c-of-m for star i 
double precision            :: stp(3)

! Local variables
double precision :: u, v, w, sini, cosi, sinphi, cosphi, s_0, t_0, p_0

sini = sin(incl)
cosi = cos(incl)
sinphi = sin(phi)
cosphi = cos(phi)

! Origin of (u, v, w) coordinate system in (s, t, p) coordinate system
s_0 =  (r_i-D)*sinphi
t_0 = -(r_i-D)*cosphi*cosi
p_0 =  (r_i-D)*cosphi*sini
u = uvw(1)
v = uvw(2)
w = uvw(3)
stp(1) =  s_0 + u
stp(2) =  t_0 + v
stp(3) =  p_0 + w

end function uvw2stp

end module coords
