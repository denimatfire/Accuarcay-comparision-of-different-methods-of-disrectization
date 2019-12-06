# Accuracy-comparision-of-different-methods-of-disrectization
A sharp 2-D Gaussian cone described by the scalar φ is rotated about the origin (x = 0, y = 0) in a unit domain −0.5 ≤ x, y ≤ 0.5. The motion of the scalar is governed by the convection diffusion equation
The velocity is chosen as u = −y, v = x, which makes the scalar rotate about the origin. The
diffusion term is switched off by taking Γ = 0. Thus an initial φ-profile remains undistorted as
it rotates. Time-integrate the convection diffusion equation for three complete rotations of the
cone, i.e., up to t = 6π considering the following

Initial profile for φ is φ(x, y, t = 0) = a exp[−ω{(x − xc)2 + (y − yc)2}] with a =
5, ω = 1500, xc = −0.25, yc = 0. Thus, peak value of the cone should remain φmax = 5
irrespective of number of rotations it completes

Use 2nd-order Adams Bashforth method and 4th-order Runge-Kutta scheme for the time
integration with ∆t = 10−3

Use the following for the spatial discretization of the convective terms: (a) 1st-order
upwind scheme, (b) 2nd-order central difference scheme, (c) 3rd-order upwind scheme,
(d) 4th-order central difference scheme, (e) the 4th-order high resolution central explicit
scheme of Tam and Webb (J. Comput. Phys. 1993), (f) 5th-order symmetric upwind
scheme and (g) the 6th-order tridiagonal central compact scheme

Use the following for the spatial discretization of the convective terms: (a) 1st-order
upwind scheme, (b) 2nd-order central difference scheme, (c) 3rd-order upwind scheme,
(d) 4th-order central difference scheme, (e) the 4th-order high resolution central explicit
scheme of Tam and Webb (J. Comput. Phys. 1993), (f) 5th-order symmetric upwind
scheme and (g) the 6th-order tridiagonal central compact scheme

Use TDMA for the compact scheme
