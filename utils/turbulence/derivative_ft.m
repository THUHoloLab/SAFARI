function der = derivative_ft(g, delta, n)
% function der = derivative_ft(g, delta, n)

    N = length(g);   % number of samples in g
    % grid spacing in the frequency domain
    F = 1/(N*delta);
    f_X = (-N/2 : N/2-1) * F;   % frequency values
    
    der = ift((i*2*pi*f_X).^n .* ft(g, delta), F);