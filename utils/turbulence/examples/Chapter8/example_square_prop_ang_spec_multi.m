% example_square_prop_ang_spec_multi.m

D1 = 2e-3;   % diameter of the source aperture [m]
D2 = 6e-3;   % diameter of the observation aperture [m]
wvl = 1e-6;  % optical wavelength [m]
k = 2*pi / wvl; % optical wavenumber [rad/m]
z = 2;     % propagation distance [m]
delta1 = D1/30; % source-plane grid spacing [m]
deltan = D2/30; % observation-plane grid spacing [m]
N = 128;        % number of grid points
n = 5;          % number of partial propagations
% switch from total distance to individual distances
z = (1:n) * z / n;
% source-plane coordinates
[x1 y1] = meshgrid((-N/2 : N/2-1) * delta1);
ap = rect(x1/D1) .* rect(y1/D1);    % source aperture
[x2 y2 Uout] = ...
    ang_spec_multi_prop_vac(ap, wvl, delta1, deltan, z);

% analytic result for y2=0 slice
Dz = z(end); % switch back to total distance
Uout_an ...
    = fresnel_prop_square_ap(x2(N/2+1,:), 0, D1, wvl, Dz);