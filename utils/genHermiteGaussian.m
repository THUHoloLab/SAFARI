function U = genHermiteGaussian(X,Y,z,wavlen,w0,m,n)
% =========================================================================
% Generate Hermite-Gaussian beams.
% -------------------------------------------------------------------------
% Input:    - X / Y  : 2D coordinates (mm).
%           - z      : Axial location (mm).
%           - wavlen : Wavelength (mm).
%           - w0     : Waist radius (mm).
%           - m / n  : Mode indices.
% Output:   - U      : Complex amplitude of the wavefront.
% =========================================================================

k = 2*pi/wavlen;
rho = sqrt(X.^2+Y.^2);
zr = pi*w0^2/wavlen;
w = w0*sqrt(1+(z/zr)^2);
Hx = polyval(myHermite(m),sqrt(2)*X/w);
Hy = polyval(myHermite(n),sqrt(2)*Y/w);
rc = sqrt(2^(1-n-m)/(pi*factorial(n)*factorial(m)))/w;

% calculate the wavefront
U = rc * Hx.*Hy.*exp(1i*(n+m+1)*atan(z/zr)).*exp(-rho.^2/w^2).*exp(-1i*k*rho.^2/(2*R(z)))*exp(1i*k*z);

% normalization
U = U/sqrt(sum(sum(abs(U).^2)));

function v = R(z)
if z == 0
    v = inf;
else
    v = z*(1 + (zr./z).^2);
end
end

end