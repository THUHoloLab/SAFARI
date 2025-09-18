% example_square_prop_ang_spec.m

D1 = 2e-3;   % diameter of the source aperture [m]
D2 = 4e-3;   % diameter of the observation aperture [m]
wvl = 1e-6;  % optical wavelength [m]
k = 2*pi / wvl;
Dz = 0.1;     % propagation distance [m]
delta1 = 9.4848e-6;
delta2 = 28.1212e-6;
Nmin = D1/(2*delta1) + D2/(2*delta2) ...
    + (wvl*Dz)/(2*delta1*delta2);
% bump N up to the next power of 2 for efficient FFT
N = 2^ceil(log2(Nmin));

[x1 y1] = meshgrid((-N/2 : N/2-1) * delta1);
ap =rect(x1/D1) .* rect(y1/D1);
[x2 y2 Uout] = ang_spec_prop(ap, wvl, delta1, delta2, Dz);

% analytic result for y2=0 slice
Uout_an ...
    = fresnel_prop_square_ap(x2(N/2+1,:), 0, D1, wvl, Dz);