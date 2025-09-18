function Z = GenZernike(n,m,N,dx,R_max)
    
    % Define Coordinate Grid %
    x         = (-N/2 : N/2-1).*dx ; % Define grid spacing
    x         = x + dx/2           ; % Centre grid
    [X,Y]     = meshgrid(x)        ; % Define Cartesian grid coordinates
    [Phi,Rho] = cart2pol(X,Y)      ; % Define polar grid coordinates 
    Rho = Rho./max(max(R_max))     ; % Rescale radial coordinate so Zernike lies on a unit circle
    
    % Define Noll Indexing Number %    
    j1 = n*(n+1)/2;    
    j2 = abs(m);    
    if mod(n,4) == 2 || mod(n,4) == 3    
        if m < 0
            j3 = 0;
        elseif m >= 0 
            j3 = 1;
        end
    elseif mod(n,4) == 0 || mod(n,4) == 1
        if m > 0
            j3 = 0;
        elseif m <= 0 
            j3 = 1;
        end
    end
    j = j1 + j2 + j3 ;
    
    % Define the Radial component %    
    R = zeros(N,N);    
    for s = 0:(n-abs(m))/2
        numerator = (-1).^s .* factorial(n-s) .* Rho.^(n-(2*s));    
        denominator = factorial(s) .*factorial( (n+m)./2 - s) .*factorial( (n-m)./2 - s);
        R = R + numerator/denominator;
    end
    
    % Define the Azimuthal Component and Final Zernike Polynomial %
    Z = zeros(N,N);
    if m == 0
        Z = sqrt(n+1) .* R;
    elseif m ~= 0 && mod(j,2) == 0 % j is even and greater than 0
        Z = sqrt(n+1) .* R .*sqrt(2) .* cos(m.*Phi) ;
    elseif m ~= 0 && mod(j,2) == 1 % j is even and greater than 0
        Z = sqrt(n+1) .* R .*sqrt(2) .* sin(m.*Phi) ;
    end
    
    % Eliminate areas outside of unit circle
    Z(Rho>1) = 0;
end
