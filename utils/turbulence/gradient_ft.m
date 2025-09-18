function [gx gy] = gradient_ft(g, delta)
% function [gx gy] = gradient_ft(g, delta)

    N = size(g, 1);   % number of samples per side in g
    % grid spacing in the frequency domain
    F = 1/(N*delta);
    fX = (-N/2 : N/2-1) * F;   % frequency values
    [fX fY] = meshgrid(fX);
    
    gx = ift2(i*2*pi*fX .* ft2(g, delta), F);
    gy = ift2(i*2*pi*fY .* ft2(g, delta), F);