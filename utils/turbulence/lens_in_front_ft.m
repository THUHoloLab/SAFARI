function [x2 y2 Uout] ...
    = lens_in_front_ft(Uin, wvl, d1, f, d)
% function [x2 y2 U_out] ...
%     = lens_in_front_ft(Uin, wvl, d1, f, d)

    N = size(Uin, 1);   % assume square grid
    k = 2*pi/wvl;    % optical wavevector
    fX = (-N/2 : 1 : N/2 - 1) / (N * d1);
    % observation plane coordinates
    [x2 y2] = meshgrid(wvl * f * fX);
    clear('fX');
    
    % evaluate the Fresnel-Kirchhoff integral but with
    % the quadratic phase factor inside cancelled by the
    % phase of the lens
    Uout = 1 / (i*wvl*f)...
        .* exp(i*k/(2*f) * (1-d/f) * (x2.^2 + y2.^2)) ...
        .* ft2(Uin, d1);