% example_square_one_step_prop_samp.m

D1 = 2e-3;  % diam of the source aperture [m]
D2 = 3e-3;  % diam of the obs-plane region of interest [m]
delta1 = D1 / 50;   % want at least 50 grid pts across ap
wvl = 1e-6;  % optical wavelength [m]
k = 2*pi / wvl;
Dz = 0.5;        % propagation distance [m]
% minimum number of grid points
Nmin = D1 * wvl*Dz / (delta1 * (wvl*Dz - D2*delta1));
N = 2^ceil(log2(Nmin));    % number of grid pts per side
% source plane
[x1 y1] = meshgrid((-N/2 : N/2-1) * delta1);
ap = rect(x1/D1) .* rect(y1/D1);
% simulate the propagation
[x2 y2 Uout] = one_step_prop(ap, wvl, delta1, Dz);

% analytic result for y2=0 slice
Uout_an ...
    = fresnel_prop_square_ap(x2(N/2+1,:), 0, D1, wvl, Dz);