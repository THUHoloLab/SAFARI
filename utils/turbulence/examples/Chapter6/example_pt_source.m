% example_pt_source.m

D = 8e-3;   % diameter of the observation aperture [m]
wvl = 1e-6;  % optical wavelength [m]
k = 2*pi / wvl; % optical wavenumber [rad/m]
Dz = 1;     % propagation distance [m]
arg = D/(wvl*Dz);
delta1 = 1/(10*arg); % source-plane grid spacing [m]
delta2 = D/100; % observation-plane grid spacing [m]
N = 1024;        % number of grid points
% source-plane coordinates
[x1 y1] = meshgrid((-N/2 : N/2-1) * delta1);
[theta1 r1] = cart2pol(x1, y1);
pt = exp(-i*k/(2*Dz) * r1.^2) * arg^2 ...
    .* sinc(arg*x1) .* sinc(arg*y1);
[x2 y2 Uout] = ang_spec_prop(pt, wvl, delta1, delta2, Dz);