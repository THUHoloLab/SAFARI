% example_fraunhofer_circ.m

N = 512;    % number of grid points per side
L = 7.5e-3;   % total size of the grid [m]
d1 = L / N;  % source-plane grid spacing [m]
D = 1e-3;   % diameter of the aperture [m]
wvl = 1e-6;  % optical wavelength [m]
k = 2*pi / wvl;
Dz = 20;     % propagation distance [m]

[x1 y1] = meshgrid((-N/2 : N/2-1) * d1);
Uin = circ(x1, y1, D);
[Uout x2 y2] = fraunhofer_prop(Uin, wvl, d1, Dz);

% analytic result
Uout_th = exp(i*k/(2*Dz)*(x2.^2+y2.^2)) ...
    / (i*wvl*Dz) * D^2*pi/4 ...
    .* jinc(D*sqrt(x2.^2+y2.^2)/(wvl*Dz));