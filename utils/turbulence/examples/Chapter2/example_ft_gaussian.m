% example_ft_gaussian.m

% function values to be used in DFT
L = 5;     % spatial extent of the grid
N = 32;     % number of samples
delta = L / N; % sample spacing
x = (-N/2 : N/2-1) * delta;
f = (-N/2 : N/2-1) / (N*delta);
a = 1;
% sampled function & its DFT
g_samp = exp(-pi*a*x.^2); % function samples
g_dft = ft(g_samp, delta); % DFT
% analytic function & its continuous FT
M = 1024;
x_cont = linspace(x(1), x(end), M);
f_cont = linspace(f(1), f(end), M);
g_cont = exp(-pi*a*x_cont.^2);
g_ft_cont = exp(-pi*f_cont.^2/a)/a;