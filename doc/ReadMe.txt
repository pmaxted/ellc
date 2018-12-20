
 The coordiate system used in the program changed in version 1.5.0 compared to
the description in Maxted A&A 591, A111, 2016. This was required to fix a bug
with the flux-weighted radial velocity in cases where detailed
heating/reflection was used. 

 See below for a description of the coordinate systems displayed in coords.png

-Pierre Maxted
 September 2016


 Coordinate transformation for ellc 

 Assumes that the stellar rotation axis and orbital angular momentum axes are
 parallel.
 
 The star is defined by an ellipsoid which is offset by a distance D from the
 centre-of-mass of the star towards the companion.

 Star coordinate system is (x, y, z) with origin at the star's centre-of-mass
 - x towards the companion
 - y in the orbital plane with star motion towards negative y
 - z parallel to rotation/orbital angular momentum vector

 Ellipsoid coordinate system is (x', y', z') with origin at the centre of the
 ellipsoid, i.e., x' = x - D, y' = y, z' = z.

 Observer coordinate system (u, v, w) with origin at the centre-of-mass of the
 binary
 - u in the plane of the sky perpendicular to the angular momentum vector
 - v in the plane of the sky parallel to the angular momentum vector
 - w pointing towards the observer

 The ellipsoid defining the star projected onto the plane of the sky is an
 ellipse. The projected coordinate system is (s, t, p) with the same origin as
 (x',y',z') and the same orientation as (u, v, w).
 
 The inclination, incl, is the angle between the orbital angular momentum
 vector and the line of sight.

 The angle phi is measured from the x-axis to the projection of the w-axis
 in the x-y plane.

