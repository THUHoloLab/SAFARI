% example_strfcn2_rect.m

N = 256;     % number of samples
L = 16;      % grid size [m]
delta = L / N;  % sample spacing [m]
F = 1/L;    % frequency-domain grid spacing [1/m]
x = (-N/2 : N/2-1) * delta;
[x y] = meshgrid(x);
w = 2;      % width of rectangle
A = rect(x/w) .* rect(y/w);  % signal
mask = ones(N);
% perform discrete structure function
C = str_fcn2_ft(A, mask, delta) / delta^2;
% continuous structure function
C_cont = 2 * w^2 * (1 - tri(x/w) .* tri(y/w));