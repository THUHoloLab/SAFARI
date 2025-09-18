function [u2] = AngularSpecProp(u1,Lx,Ly,lambda,z)

    % u1 is the original field
    % L is length of original plane
    % lambda is the wavelength
    % z is the propagation distance
    
    % Create coordinate grids
    
    Mx = size(u1,2); % get size of original field
    My= size(u1,1); % get size of original field
    
    dx = Lx/Mx; % spacing of the starting spatial grid
    dy = Ly/My; % spacing of the starting spatial grid
    fx = -1/(2*dx):1/Lx:1/(2*dx)-1/Lx;
    fy = -1/(2*dy):1/Ly:1/(2*dy)-1/Ly;% frequency coordinates
    [Fx,Fy] = meshgrid(fx,fy);
    
    % Create propogation phase in freqeucnyy domain
    
    H = exp(-1i*pi*lambda*z*(Fx.^2+Fy.^2)); % define porpagation pahse
    H = fftshift(H); % shift the propagation phase for Matlab's fft function
    
    % Propagate field
    
    U1 = fft2(fftshift(u1)); % Fourier transform original function into frequency domain
    U2 = U1.*H; % add propagation phase
    u2 = ifftshift(ifft2(U2)); % Shift function back from frequency domain to spatial domain
end