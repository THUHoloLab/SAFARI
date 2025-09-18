function [Uout x2 y2] = ...
    lens_against_ft(Uin, wvl, d1, f)
% function [Uout x2 y2] = ...
%     lens_against_ft(Uin, wvl, d1, f)

    N = size(Uin, 1);   % assume square grid
    k = 2*pi/wvl;    % optical wavevector
    fX = (-N/2 : 1 : N/2 - 1) / (N * d1);
    % observation plane coordinates
    [x2 y2] = meshgrid(wvl * f * fX);
    clear('fX');
    
    % evaluate the Fresnel-Kirchhoff integral but with
    % the quadratic phase factor inside cancelled by the
    % phase of the lens
    Uout = exp(i*k/(2*f)*(x2.^2 + y2.^2)) ...
        / (i*wvl*f) .* ft2(Uin, d1);