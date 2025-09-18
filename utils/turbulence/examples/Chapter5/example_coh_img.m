% example_coh_img.m

N = 256;    % number of grid points per side
L = 0.1;      % total size of the grid [m]
D = 0.07;   % diameter of pupil [m]
delta = L / N;  % grid spacing [m]
wvl = 1e-6; % optical wavelength [m]
z = 0.25;   % image distance [m]
% pupil-plane coordinates
[x y] = meshgrid((-N/2 : N/2-1) * delta);
[theta r] = cart2pol(x, y);
% wavefront aberration
W = 0.05 * zernike(4, 2*r/D, theta);
% complex pupil function
P = circ(x, y, D) .* exp(i * 2*pi * W);
% amplitude spread function
h = ft2(P, delta);
delta_u = wvl * z / (N*delta);
% image-plane coordinates
[u v] = meshgrid((-N/2 : N/2-1) * delta_u);
% object (same coordinates as h)
obj = (rect((u-1.4e-4)/5e-5) + rect(u/5e-5) ...
    + rect((u+1.4e-4)/5e-5)) .* rect(v/2e-4);
% convolve the object with the ASF to simulate imaging
img = myconv2(obj, h, 1);