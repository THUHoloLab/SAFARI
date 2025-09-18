% example_zernike_projection.m

N = 32;    % number of grid points per side
L = 2;      % total size of the grid [m]
delta = L / N;  % grid spacing [m]
% cartesian & polar coordinates
[x y] = meshgrid((-N/2 : N/2-1) * delta);
[theta r] = cart2pol(x, y);
% unit circle aperture
ap = circ(x, y, 2);
% 3 Zernike modes
z2 = zernike(2, r, theta) .* ap;
z4 = zernike(4, r, theta) .* ap;
z21 = zernike(21, r, theta) .* ap;
% create the aberration
W = 0.5 *  z2 + 0.25 * z4 - 0.6 * z21;
% find only grid points within the aperture
idx = logical(ap);
% perform linear indexing in column-major order
W = W(idx);
Z = [z2(idx) z4(idx) z21(idx)];
% solve the system of equations to compute coefficients
A = Z \ W