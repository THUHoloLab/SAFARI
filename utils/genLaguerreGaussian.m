function U = genLaguerreGaussian(X,Y,z,wavlen,w0,l,p)
% =========================================================================
% Generate Laguerre-Gaussian beams.
% -------------------------------------------------------------------------
% Input:    - X / Y  : 2D coordinates (mm).
%           - z      : Axial location (mm).
%           - wavlen : Wavelength (mm).
%           - w0     : Waist radius (mm).
%           - l      : Azimuthal index.
%           - p      : Radial index.
% Output:   - U      : Complex amplitude of the wavefront.
% =========================================================================

zR = pi*w0^2/wavlen;
N  = abs(l) + 2*p;
k  = 2*pi/wavlen;

w   = @(z) w0*sqrt(1+(z/zR).^2);
psi = @(z) (N+1)*atan(z/zR);

% convert to polar coordinates
[phi,rho] = cart2pol(X,Y);

% define the normalizing constant
C = sqrt(2*factorial(p)/(pi*factorial(p+abs(l))));

% calculate the wavefront
U = C .* 1./w(z) * (rho*sqrt(2)./w(z)).^abs(l) .*exp(-rho.^2./w(z).^2) ...
    .* myLaguerre(p,abs(l),2*rho.^2./w(z).^2)...
    .*exp(-1i*k*rho.^2./2./R(z)).*exp(-1i*l*phi).*exp(1i*psi(z));


function v = R(z)
if z == 0
    v = inf;
else
    v = z*(1+(zR./z).^2);
end
end

end