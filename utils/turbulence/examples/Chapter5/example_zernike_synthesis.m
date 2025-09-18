% example_zernike_synthesis.m

N = 40;     % number of grid points per side
L = 2;      % total size of the grid [m]
delta = L / N;  % grid spacing [m]
% cartesian & polar coordinates
[x y] = meshgrid((-N/2 : N/2-1) * delta);
[theta r] = cart2pol(x, y);
% unit circle aperture
ap = circ(x, y, 2);
% indices of grid points in aperture
idxAp = logical(ap);
% create atmospheric phase screen
r0 = L / 20;
screen = ft_phase_screen(r0, N, delta, inf, 0) ...
    / (2*pi) .* ap;
W = screen(idxAp);   % perform linear indexing

%%% analyze screen
nModes = 100;   % number of Zernike modes
% create matrix of Zernike polynomial values
Z = zeros(numel(W), nModes);
for idx = 1 : nModes
    temp = zernike(idx, r, theta);
    Z(:,idx) = temp(idxAp);
end
% compute mode coefficients
A = Z \ W;
% synthesize mode-limited screen
W_prime = Z*A;
% reshape mode-limited screen into 2-D for display
scr = zeros(N);
scr(idxAp) = W_prime;