function U = genBessel(X,Y,z,wavlen,n,theta)
% =========================================================================
% Generate Bessel beams.
% -------------------------------------------------------------------------
% Input:    - X / Y  : 2D coordinates (mm).
%           - z      : Axial location (mm).
%           - wavlen : Wavelength (mm).
%           - n      : Topological charge.
%           - theta  : Axicon angle (rad).
% Output:   - U      : Complex amplitude of the wavefront.
% =========================================================================

k = 2*pi/wavlen;
[phi,rho] = cart2pol(X,Y);

kz = k*cos(theta);   % longitudinal wavenumber
kr = k*sin(theta);   % transverse wavenumber

U = exp(1i*kz*z) .* besselj(n,kr*rho) .* exp(1i*n*phi);

end