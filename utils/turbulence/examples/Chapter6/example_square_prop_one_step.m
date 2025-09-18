% example_square_prop_one_step.m

N = 1024;    % number of grid points per side
L = 1e-2;   % total size of the grid [m]
delta1 = L / N;  % grid spacing [m]
D = 2e-3;   % diameter of the aperture [m]
wvl = 1e-6;  % optical wavelength [m]
k = 2*pi / wvl;
Dz = 1;     % propagation distance [m]

[x1 y1] = meshgrid((-N/2 : N/2-1) * delta1);
ap = rect(x1/D) .* rect(y1/D);
[x2 y2 Uout] = one_step_prop(ap, wvl, delta1, Dz);

% analytic result for y2=0 slice
Uout_an ...
    = fresnel_prop_square_ap(x2(N/2+1,:), 0, D, wvl, Dz);