% example_pt_source_vac_prop.m

delta1 = d1;    % source-plane grid spacing [m]
deltan = d2;    % observation-plane grid spacing [m]
n = nscr;         % number of planes

% coordinates
[x1 y1] = meshgrid((-N/2 : N/2-1) * delta1);
[theta1 r1] = cart2pol(x1, y1);

% point source
pt = exp(-i*k/(2*R) * r1.^2) / D1^2 ...
    .* sinc(x1/D1) .* sinc(y1/D1) ...
    .* exp(-(r1/(4*D1)).^2);
% partial prop planes
z = (1 : n-1) * Dz / (n-1);

% simulate vacuum propagation
sg = exp(-(x1/(0.47*N*d1)).^16) ...
    .* exp(-(y1/(0.47*N*d1)).^16);
t = repmat(sg, [1 1 n]);
[xn yn Uvac] = ang_spec_multi_prop(pt, wvl, ...
    delta1, deltan, z, t);
% collimate the beam
Uvac = Uvac .* exp(-i*pi/(wvl*R)*(xn.^2+yn.^2));