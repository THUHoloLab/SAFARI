function [x2 y2 Uout] ...
    = lens_behind_ft(Uin, wvl, d1, f)
% function [x2 y2 Uout] ...
%     = lens_behind_ft(Uin, wvl, d1, d, f)

    N = size(Uin, 1);   % assume square grid
    k = 2*pi/wvl;    % optical wavevector
    fX = (-N/2 : 1 : N/2 - 1) / (N * d1);
    % observation plane coordinates
    [x2 y2] = meshgrid(wvl * d * fX);
    clear('fX');
    
    % evaluate the Fresnel-Kirchhoff integral but with
    % the quadratic phase factor inside cancelled by the
    % phase of the lens
    Uout = f/d * 1  / (i*wvl*d)...
        .* exp(i*k/(2*d)*(x2.^2 + y2.^2)) .* ft2(Uin, d1);