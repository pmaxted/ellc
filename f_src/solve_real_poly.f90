MODULE Solve_Real_Poly
! CACM Algorithm 493 by Jenkins & Traub

! Compliments of netlib   Sat Jul 26 11:57:43 EDT 1986
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-06-02  Time: 10:42:22
! 
! Code converted to a form suitable for calling from python using f2py by
! restricting polynomial degree to 4 (enough for this application).
! Also converted all REALs to DOUBLE PRECISION
! p.maxted@keele.ac.uk
! 29 Nov 2015

IMPLICIT NONE
! INTEGER, PARAMETER  :: dp = SELECTED_REAL_KIND(14, 60)

! COMMON /global/ p, qp, k, qk, svk, sr, si, u, v, a, b, c, d, a1,  &
!    a2, a3, a6, a7, e, f, g, h, szr, szi, lzr, lzi, eta, are, mre, n, nn

DOUBLE PRECISION, SAVE  :: p(5), qp(5), k(5), qk(5), svk(5)
DOUBLE PRECISION, SAVE :: sr, si, u, v, a, b, c, d, a1, a2, a3, a6,  &
                          a7, e, f, g, h, szr, szi, lzr, lzi
DOUBLE PRECISION, SAVE             :: eta, are, mre
INTEGER, SAVE          :: n, nn

!PRIVATE
!PUBLIC  :: dp, rpoly


CONTAINS


SUBROUTINE rpoly(op, degree, zeror, zeroi, fail)

! Finds the zeros of a real polynomial
! op  - double precision vector of coefficients in order of
!       decreasing powers.
! degree   - integer degree of polynomial.
! zeror, zeroi - output double precision vectors of real and imaginary parts
!                of the zeros.
! fail  - output logical parameter, true only if leading coefficient is zero
!         or if rpoly has found fewer than degree zeros.
!         In the latter case degree is reset to the number of zeros found.

! To change the size of polynomials which can be solved, reset the dimensions
! of the arrays in the common area and in the following declarations.
! The subroutine uses single precision calculations for scaling, bounds and
! error calculations.  All calculations for the iterations are done in
! double precision.


DOUBLE PRECISION, INTENT(IN)    :: op(5)
INTEGER, INTENT(IN OUT)  :: degree
DOUBLE PRECISION, INTENT(OUT)   :: zeror(4), zeroi(4)
LOGICAL, INTENT(OUT)     :: fail

DOUBLE PRECISION :: temp(5)
DOUBLE PRECISION :: pt(5)

DOUBLE PRECISION :: t, aa, bb, cc, factor
DOUBLE PRECISION      :: lo, MAX, MIN, xx, yy, cosr, sinr, xxx, x, sc, bnd,  &
             xm, ff, df, dx, infin, smalno, base
INTEGER   :: cnt, nz, i, j, jj, l, nm1
LOGICAL   :: zerok

! The following statements set machine constants used in various parts of the
! program.  The meaning of the four constants are...
! eta     the maximum relative representation error which can be described
!         as the smallest positive floating point number such that
!         1.d0+eta is greater than 1.
! infiny  the largest floating-point number.
! smalno  the smallest positive floating-point number if the exponent range
!         differs in single and double precision then smalno and infin should
!         indicate the smaller range.
! base    the base of the floating-point number system used.

base = RADIX(0.0d0)
eta = EPSILON(1.0d0)
infin = HUGE(0.0d0)
smalno = TINY(0.0d0)

temp = 0

! are and mre refer to the unit error in + and * respectively.
! They are assumed to be the same as eta.
are = eta
mre = eta
lo = smalno / eta

! Initialization of constants for shift rotation
xx = SQRT(0.5)
yy = -xx
cosr = -0.069756474d0
sinr = 0.99756405d0
fail = .false.
n = degree
nn = n + 1

! Algorithm fails if the leading coefficient is zero.
IF (op(1) == 0.d0) THEN
  fail = .true.
  degree = 0
  RETURN
endif

! Remove the zeros at the origin if any
10 IF (op(nn) == 0.0D0) THEN
  j = degree - n + 1
  zeror(j) = 0.d0
  zeroi(j) = 0.d0
  nn = nn - 1
  n = n - 1
  GO TO 10
endif

! Make a copy of the coefficients
p(1:nn) = op(1:nn)

! Start the algorithm for one zero
30 IF (n <= 2) THEN
  IF (n < 1) RETURN

! calculate the final zero or pair of zeros
  IF (n /= 2) THEN
    zeror(degree) = -p(2) / p(1)
    zeroi(degree) = 0.0D0
    RETURN
  endif
  CALL quad(p(1), p(2), p(3), zeror(degree-1), zeroi(degree-1),  &
            zeror(degree), zeroi(degree))
  RETURN
endif

! Find largest and smallest moduli of coefficients.
MAX = 0
MIN = infin
DO  i = 1, nn
  x = ABS( DBLE(p(i)) )
  IF (x > MAX) MAX = x
  IF (x /= 0. .AND. x < MIN) MIN = x
END DO

! Scale if there are large or very small coefficients computes a scale
! factor to multiply the coefficients of the polynomial.
! The scaling is done to avoid overflow and to avoid undetected underflow
! interfering with the convergence criterion.
! The factor is a power of the base
sc = lo / MIN
IF (sc <= 1.0) THEN
  IF (MAX < 10.) GO TO 60
  IF (sc == 0.) sc = smalno
ELSE
  IF (infin/sc < MAX) GO TO 60
endif
l = INT(LOG(sc) / LOG(base) + .5)
factor = (base*1.0D0) ** l
IF (factor /= 1.d0) THEN
  p(1:nn) = factor * p(1:nn)
endif

! compute lower bound on moduli of zeros.
60 pt(1:nn) = ABS(p(1:nn))
pt(nn) = -pt(nn)

! compute upper estimate of bound
x = EXP((LOG(-pt(nn)) - LOG(pt(1))) / n)
IF (pt(n) /= 0.) THEN
! if newton step at the origin is better, use it.
  xm = -pt(nn) / pt(n)
  IF (xm < x) x = xm
endif

! chop the interval (0,x) until ff .le. 0
80 xm = x * .1
ff = pt(1)
DO  i = 2, nn
  ff = ff * xm + pt(i)
END DO
IF (ff > 0.) THEN
  x = xm
  GO TO 80
endif
dx = x

! do newton iteration until x converges to two decimal places
100 IF (ABS(dx/x) > .005) THEN
  ff = pt(1)
  df = ff
  DO  i = 2, n
    ff = ff * x + pt(i)
    df = df * x + ff
  END DO
  ff = ff * x + pt(nn)
  dx = ff / df
  x = x - dx
  GO TO 100
endif
bnd = x

! compute the derivative as the intial k polynomial
! and do 5 steps with no shift
nm1 = n - 1
DO  i = 2, n
  k(i) = (nn-i) * p(i) / n
END DO
k(1) = p(1)
aa = p(nn)
bb = p(n)
zerok = k(n) == 0.d0
DO  jj = 1, 5
  cc = k(n)
  IF (.NOT.zerok) THEN
! use scaled form of recurrence if value of k at 0 is nonzero
    t = -aa / cc
    DO  i = 1, nm1
      j = nn - i
      k(j) = t * k(j-1) + p(j)
    END DO
    k(1) = p(1)
    zerok = ABS(k(n)) <= ABS(bb) * eta * 10.
  ELSE
! use unscaled form of recurrence
    DO  i = 1, nm1
      j = nn - i
      k(j) = k(j-1)
    END DO
    k(1) = 0.d0
    zerok = k(n) == 0.d0
  endif
END DO

! save k for restarts with new shifts
temp(1:n) = k(1:n)

! loop to select the quadratic  corresponding to each
! new shift
DO  cnt = 1, 20
! Quadratic corresponds to a double shift to a non-real point and its complex
! conjugate.  The point has modulus bnd and amplitude rotated by 94 degrees
! from the previous shift
  xxx = cosr * xx - sinr * yy
  yy = sinr * xx + cosr * yy
  xx = xxx
  sr = bnd * xx
  si = bnd * yy
  u = -2.0D0 * sr
  v = bnd

! second stage calculation, fixed quadratic
  CALL fxshfr(20*cnt,nz)
  IF (nz /= 0) THEN

! The second stage jumps directly to one of the third stage iterations and
! returns here if successful.
! Deflate the polynomial, store the zero or zeros and return to the main
! algorithm.
    j = degree - n + 1
    zeror(j) = szr
    zeroi(j) = szi
    nn = nn - nz
    n = nn - 1
    p(1:nn) = qp(1:nn)
    IF (nz == 1) GO TO 30
    zeror(j+1) = lzr
    zeroi(j+1) = lzi
    GO TO 30
  endif

! If the iteration is unsuccessful another quadratic
! is chosen after restoring k
  k(1:nn) = temp(1:nn)
END DO

! Return with failure if no convergence with 20 shifts
fail = .true.
degree = degree - n
RETURN
END SUBROUTINE rpoly




SUBROUTINE fxshfr(l2, nz)

! Computes up to  l2  fixed shift k-polynomials, testing for convergence in
! the linear or quadratic case.  Initiates one of the variable shift
! iterations and returns with the number of zeros found.
! l2 - limit of fixed shift steps
! nz - number of zeros found

INTEGER, INTENT(IN)   :: l2
INTEGER, INTENT(OUT)  :: nz

DOUBLE PRECISION :: svu, svv, ui, vi, s
DOUBLE PRECISION :: betas, betav, oss, ovv, ss, vv, ts, tv, ots, otv, tvv, tss
INTEGER   :: TTYPE, j, iflag
LOGICAL   :: vpass, spass, vtry, stry

nz = 0
betav = .25
betas = .25
oss = sr
ovv = v

! Evaluate polynomial by synthetic division
CALL quadsd(nn, u, v, p, qp, a, b)
CALL calcsc(TTYPE)
DO  j = 1, l2
! calculate next k polynomial and estimate v
  CALL nextk(TTYPE)
  CALL calcsc(TTYPE)
  CALL newest(TTYPE, ui, vi)
  vv = vi

! Estimate s
  ss = 0
  IF (k(n) /= 0.d0) ss = -p(nn) / k(n)
  tv = 1
  ts = 1
  IF (j /= 1 .AND. TTYPE /= 3) THEN
! Compute relative measures of convergence of s and v sequences
    IF (vv /= 0.) tv = ABS((vv-ovv)/vv)
    IF (ss /= 0.) ts = ABS((ss-oss)/ss)

! If decreasing, multiply two most recent convergence measures
    tvv = 1
    IF (tv < otv) tvv = tv * otv
    tss = 1
    IF (ts < ots) tss = ts * ots

! Compare with convergence criteria
    vpass = tvv < betav
    spass = tss < betas
    IF (spass .OR. vpass) THEN

! At least one sequence has passed the convergence test.
! Store variables before iterating
      svu = u
      svv = v
      svk(1:n) = k(1:n)
      s = ss

! Choose iteration according to the fastest converging sequence
      vtry = .false.
      stry = .false.
      IF (spass .AND. ((.NOT.vpass) .OR. tss < tvv)) GO TO 40
      20 CALL quadit(ui, vi, nz)
      IF (nz > 0) RETURN

! Quadratic iteration has failed. flag that it has
! been tried and decrease the convergence criterion.
      vtry = .true.
      betav = betav * .25

! Try linear iteration if it has not been tried and
! the s sequence is converging
      IF (stry.OR.(.NOT.spass)) GO TO 50
      k(1:n) = svk(1:n)
      40 CALL realit(s, nz, iflag)
      IF (nz > 0) RETURN

! Linear iteration has failed.  Flag that it has been
! tried and decrease the convergence criterion
      stry = .true.
      betas = betas * .25
      IF (iflag /= 0) THEN

! If linear iteration signals an almost double real
! zero attempt quadratic interation
        ui = -(s+s)
        vi = s * s
        GO TO 20
      endif

! Restore variables
      50 u = svu
      v = svv
      k(1:n) = svk(1:n)

! Try quadratic iteration if it has not been tried
! and the v sequence is converging
      IF (vpass .AND. (.NOT.vtry)) GO TO 20

! Recompute qp and scalar values to continue the second stage
      CALL quadsd(nn, u, v, p, qp, a, b)
      CALL calcsc(TTYPE)
    endif
  endif
  ovv = vv
  oss = ss
  otv = tv
  ots = ts
END DO
RETURN
END SUBROUTINE fxshfr




SUBROUTINE quadit(uu, vv, nz)

! Variable-shift k-polynomial iteration for a quadratic factor, converges
! only if the zeros are equimodular or nearly so.
! uu,vv - coefficients of starting quadratic
! nz - number of zero found

DOUBLE PRECISION, INTENT(IN)  :: uu
DOUBLE PRECISION, INTENT(IN)  :: vv
INTEGER, INTENT(OUT)   :: nz

DOUBLE PRECISION :: ui, vi
DOUBLE PRECISION :: mp, omp, ee, relstp, t, zm
INTEGER   :: TTYPE, i, j
LOGICAL   :: tried

nz = 0
tried = .false.
u = uu
v = vv
j = 0
relstp = 0
omp = 0

! Main loop
10 CALL quad(1.d0, u, v, szr, szi, lzr, lzi)

! Return if roots of the quadratic are real and not
! close to multiple or nearly equal and  of opposite sign.
IF (ABS(ABS(szr)-ABS(lzr)) > .01D0*ABS(lzr)) RETURN

! Evaluate polynomial by quadratic synthetic division
CALL quadsd(nn, u, v, p, qp, a, b)
mp = ABS(a-szr*b) + ABS(szi*b)

! Compute a rigorous  bound on the rounding error in evaluting p
zm = SQRT(ABS( DBLE(v)))
ee = 2. * ABS( DBLE(qp(1)))
t = -szr * b
DO  i = 2, n
  ee = ee * zm + ABS( DBLE(qp(i)) )
END DO
ee = ee * zm + ABS( DBLE(a) + t)
ee = (5.*mre+4.*are) * ee - (5.*mre+2.*are) * (ABS( DBLE(a)  + t) +  &
     ABS( DBLE(b))*zm) + 2. * are * ABS(t)

! Iteration has converged sufficiently if the
! polynomial value is less than 20 times this bound
IF (mp <= 20.*ee) THEN
  nz = 2
  RETURN
endif
j = j + 1

! Stop iteration after 20 steps
IF (j > 20) RETURN
IF (j >= 2) THEN
  IF (.NOT.(relstp > .01 .OR. mp < omp .OR. tried)) THEN

! A cluster appears to be stalling the convergence.
! five fixed shift steps are taken with a u,v close to the cluster
    IF (relstp < eta) relstp = eta
    relstp = SQRT(relstp)
    u = u - u * relstp
    v = v + v * relstp
    CALL quadsd(nn, u, v, p, qp, a, b)
    DO  i = 1, 5
      CALL calcsc(TTYPE)
      CALL nextk(TTYPE)
    END DO
    tried = .true.
    j = 0
  endif
endif
omp = mp

! Calculate next k polynomial and new u and v
CALL calcsc(TTYPE)
CALL nextk(TTYPE)
CALL calcsc(TTYPE)
CALL newest(TTYPE,ui,vi)

! If vi is zero the iteration is not converging
IF (vi == 0.d0) RETURN
relstp = ABS((vi-v)/vi)
u = ui
v = vi
GO TO 10
END SUBROUTINE quadit




SUBROUTINE realit(sss,nz,iflag)

! Variable-shift h polynomial iteration for a real
! zero.
! sss   - starting iterate
! nz    - number of zero found
! iflag - flag to indicate a pair of zeros near real axis.

DOUBLE PRECISION, INTENT(IN OUT)  :: sss
INTEGER, INTENT(OUT)       :: nz, iflag

DOUBLE PRECISION :: pv, kv, t, s
DOUBLE PRECISION :: ms, mp, omp, ee
INTEGER   :: i, j

nz = 0
s = sss
iflag = 0
j = 0
t = 0
omp = 0

! Main loop
10 pv = p(1)

! Evaluate p at s
qp(1) = pv
DO  i = 2, nn
  pv = pv * s + p(i)
  qp(i) = pv
END DO
mp = ABS(pv)

! Compute a rigorous bound on the error in evaluating p
ms = ABS(s)
ee = (mre/(are+mre)) * ABS( DBLE(qp(1)))
DO  i = 2, nn
  ee = ee * ms + ABS( DBLE(qp(i)))
END DO

! Iteration has converged sufficiently if the
! polynomial value is less than 20 times this bound
IF (mp <= 20.*((are+mre)*ee - mre*mp)) THEN
  nz = 1
  szr = s
  szi = 0.d0
  RETURN
endif
j = j + 1

! Stop iteration after 10 steps
IF (j > 10) RETURN
IF (j >= 2) THEN
  IF (ABS(t) <= .001*ABS(s-t) .AND. mp > omp) THEN
! A cluster of zeros near the real axis has been encountered,
! return with iflag set to initiate a quadratic iteration
    iflag = 1
    sss = s
    RETURN
  endif
endif

! Return if the polynomial value has increased significantly
omp = mp

! Compute t, the next polynomial, and the new iterate
kv = k(1)
qk(1) = kv
DO  i = 2, n
  kv = kv * s + k(i)
  qk(i) = kv
END DO
IF (ABS(kv) > ABS(k(n))*10.*eta) THEN
! Use the scaled form of the recurrence if the value of k at s is nonzero
  t = -pv / kv
  k(1) = qp(1)
  DO  i = 2, n
    k(i) = t * qk(i-1) + qp(i)
  END DO
ELSE
! Use unscaled form
  k(1) = 0.0D0
  DO  i = 2, n
    k(i) = qk(i-1)
  END DO
endif
kv = k(1)
DO  i = 2, n
  kv = kv * s + k(i)
END DO
t = 0.d0
IF (ABS(kv) > ABS(k(n))*10.*eta) t = -pv / kv
s = s + t
GO TO 10
END SUBROUTINE realit




SUBROUTINE calcsc(TTYPE)

! This routine calculates scalar quantities used to
! compute the next k polynomial and new estimates of
! the quadratic coefficients.
! type - integer variable set here indicating how the
! calculations are normalized to avoid overflow

INTEGER, INTENT(OUT)  :: TTYPE

! Synthetic division of k by the quadratic 1,u,v
CALL quadsd(n, u, v, k, qk, c, d)
IF (ABS(c) <= ABS(k(n))*100.*eta) THEN
  IF (ABS(d) <= ABS(k(n-1))*100.*eta) THEN
    TTYPE = 3
! type=3 indicates the quadratic is almost a factor of k
    RETURN
  endif
endif

IF (ABS(d) >= ABS(c)) THEN
  TTYPE = 2
! type=2 indicates that all formulas are divided by d
  e = a / d
  f = c / d
  g = u * b
  h = v * b
  a3 = (a+g) * e + h * (b/d)
  a1 = b * f - a
  a7 = (f+u) * a + h
  RETURN
endif
TTYPE = 1
! type=1 indicates that all formulas are divided by c
e = a / c
f = d / c
g = u * e
h = v * b
a3 = a * e + (h/c+g) * b
a1 = b - a * (d/c)
a7 = a + g * d + h * f
RETURN
END SUBROUTINE calcsc




SUBROUTINE nextk(TTYPE)

! Computes the next k polynomials using scalars computed in calcsc.

INTEGER, INTENT(IN)  :: TTYPE

DOUBLE PRECISION :: temp
INTEGER   :: i

IF (TTYPE /= 3) THEN
  temp = a
  IF (TTYPE == 1) temp = b
  IF (ABS(a1) <= ABS(temp)*eta*10.) THEN
! If a1 is nearly zero then use a special form of the recurrence
    k(1) = 0.d0
    k(2) = -a7 * qp(1)
    DO  i = 3, n
      k(i) = a3 * qk(i-2) - a7 * qp(i-1)
    END DO
    RETURN
  endif

! Use scaled form of the recurrence
  a7 = a7 / a1
  a3 = a3 / a1
  k(1) = qp(1)
  k(2) = qp(2) - a7 * qp(1)
  DO  i = 3, n
    k(i) = a3 * qk(i-2) - a7 * qp(i-1) + qp(i)
  END DO
  RETURN
endif

! Use unscaled form of the recurrence if type is 3
k(1) = 0.d0
k(2) = 0.d0
DO  i = 3, n
  k(i) = qk(i-2)
END DO
RETURN
END SUBROUTINE nextk




SUBROUTINE newest(TTYPE,uu,vv)

! Compute new estimates of the quadratic coefficients
! using the scalars computed in calcsc.

INTEGER, INTENT(IN)     :: TTYPE
DOUBLE PRECISION, INTENT(OUT)  :: uu
DOUBLE PRECISION, INTENT(OUT)  :: vv

DOUBLE PRECISION :: a4, a5, b1, b2, c1, c2, c3, c4, temp

! Use formulas appropriate to setting of type.
IF (TTYPE /= 3) THEN
  IF (TTYPE /= 2) THEN
    a4 = a + u * b + h * f
    a5 = c + (u+v*f) * d
  ELSE
    a4 = (a+g) * f + h
    a5 = (f+u) * c + v * d
  endif

! Evaluate new quadratic coefficients.
  b1 = -k(n) / p(nn)
  b2 = -(k(n-1)+b1*p(n)) / p(nn)
  c1 = v * b2 * a1
  c2 = b1 * a7
  c3 = b1 * b1 * a3
  c4 = c1 - c2 - c3
  temp = a5 + b1 * a4 - c4
  IF (temp /= 0.d0) THEN
    uu = u - (u*(c3+c2)+v*(b1*a1+b2*a7)) / temp
    vv = v * (1.+c4/temp)
    RETURN
  endif
endif

! If type=3 the quadratic is zeroed
uu = 0.d0
vv = 0.d0
RETURN
END SUBROUTINE newest




SUBROUTINE quadsd(nn, u, v, p, q, a, b)

! Divides p by the quadratic  1,u,v  placing the
! quotient in q and the remainder in a,b.

INTEGER, INTENT(IN)     :: nn
DOUBLE PRECISION, INTENT(IN)   :: u, v, p(nn)
DOUBLE PRECISION, INTENT(OUT)  :: q(nn), a, b

DOUBLE PRECISION  :: c
INTEGER    :: i

b = p(1)
q(1) = b
a = p(2) - u * b
q(2) = a
DO  i = 3, nn
  c = p(i) - u * a - v * b
  q(i) = c
  b = a
  a = c
END DO
RETURN
END SUBROUTINE quadsd




SUBROUTINE quad(a, b1, c, sr, si, lr, li)

! Calculate the zeros of the quadratic a*z**2+b1*z+c.
! The quadratic formula, modified to avoid overflow, is used to find the
! larger zero if the zeros are real and both zeros are complex.
! The smaller real zero is found directly from the product of the zeros c/a.

DOUBLE PRECISION, INTENT(IN)             :: a, b1, c
DOUBLE PRECISION, INTENT(OUT)            :: sr, si, lr, li

DOUBLE PRECISION :: b, d, e

IF (a /= 0.d0) GO TO 20
sr = 0.d0
IF (b1 /= 0.d0) sr = -c / b1
lr = 0.d0
10 si = 0.d0
li = 0.d0
RETURN

20 IF (c == 0.d0) THEN
  sr = 0.d0
  lr = -b1 / a
  GO TO 10
endif

! Compute discriminant avoiding overflow
b = b1 / 2.d0
IF (ABS(b) >= ABS(c)) THEN
  e = 1.d0 - (a/b) * (c/b)
  d = SQRT(ABS(e)) * ABS(b)
ELSE
  e = a
  IF (c < 0.d0) e = -a
  e = b * (b/ABS(c)) - e
  d = SQRT(ABS(e)) * SQRT(ABS(c))
endif
IF (e >= 0.d0) THEN

! Real zeros
  IF (b >= 0.d0) d = -d
  lr = (-b+d) / a
  sr = 0.d0
  IF (lr /= 0.d0) sr = (c/lr) / a
  GO TO 10
endif
! complex conjugate zeros
sr = -b / a
lr = sr
si = ABS(d/a)
li = -si
RETURN
END SUBROUTINE quad

END MODULE Solve_Real_Poly
