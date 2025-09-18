% example_corr2_rect_rect.m

N = 256;     % number of samples
L = 16;      % grid size [m]
delta = L / N;  % sample spacing [m]
F = 1/L;    % frequency-domain grid spacing [1/m]
x = (-N/2 : N/2-1) * delta;
[x y] = meshgrid(x);
w = 2;      % width of rectangle
A = rect(x/w) .* rect(y/w);  % signal
mask = ones(N);
% perform discrete correlation
C = corr2_ft(A, A, mask, delta);
% analytic correlation
C_cont = w^2*tri(x/w) .* tri(y/w);