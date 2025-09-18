% example_gradient_ft.m
N = 64;     % number of samples
L = 6;      % grid size [m]
delta = L/N;    % grid spacing [m]
x = (-N/2 : N/2-1) * delta;
[x y] = meshgrid(x);
g = exp(-(x.^2 + y.^2));
% discrete derivatives
[gx_samp gy_samp] = gradient_ft(g, delta);
gx_samp = real(gx_samp);
gy_samp = real(gy_samp);
% analytic derivatives
gx = -2*x.*exp(-(x.^2+y.^2));
gy = -2*y.*exp(-(x.^2+y.^2));