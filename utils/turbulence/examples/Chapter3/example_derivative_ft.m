% example_derivative_ft.m

N = 64;     % number of samples
L = 6;      % grid size [m]
delta = L/N;    % grid spacing [m]
x = (-N/2 : N/2-1) * delta;
w = 3;      % size of window (or region of interest) [m]
window = rect(x/w); % window function [m]
g = x.^5 .* window; % function
% discrete derivatives
gp_samp = real(derivative_ft(g, delta, 1)) .* window;
gpp_samp = real(derivative_ft(g, delta, 2)) .* window;
% analytic derivatives
gp = 5*x.^4 .* window;
gpp = 20*x.^3 .* window;