function [Uout x2 y2] = ...
    fraunhofer_prop(Uin, wvl, d1, Dz)
% function [Uout x2 y2] = ...
%     fraunhofer_prop(Uin, wvl, d1, Dz)

    N = size(Uin, 1);   % assume square grid
    k = 2*pi / wvl;  % optical wavevector
    fX = (-N/2 : N/2-1) / (N*d1);
    % observation-plane coordinates
    [x2 y2] = meshgrid(wvl * Dz * fX);
    clear('fX');
    Uout = exp(i*k/(2*Dz)*(x2.^2+y2.^2)) ...
        / (i*wvl*Dz) .* ft2(Uin, d1);