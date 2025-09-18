function screen = FourierPhaseScreen(N,dx,r0,rs,numsub)
    % rs - Random seed for repeatability in phase screen generation 
    % r0 - Fried parameter define correlation length of field and thus strength of the turbulence
    % N  - number of pixels in real space
    % dx - pixel size in real space
    % numsub - number of subharmonics for subharmonic sampling    

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % High Frequency Screen %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Real Space Parameters %
    Lx    = N*dx               ; % length of our grid in the spatial domain (number of pixels*length of each pixel)
    x     = (-N/2 : N/2-1).*dx ; % Creates a vector of length N (One entry for each pixel) and centres it so that the middle entry acts as the origin/centre of our one coordinate axes
    [X,Y] = meshgrid(x)        ; % creates spatial coordinate grid
    
    % Fourier Space Parameters %
    df      = 1/(Lx)              ; % spacing in the frequency domain (inverse lengt of the spatial domain as defined by a Fourier transform)
    f       = (-N/2 : N/2-1).*df  ; % Creates a vector of length N (One entry for each pixel) and centres it so that the middle entry acts as the origin/centre of our one coordinate axes
    [Fx,Fy] =  meshgrid(f)        ; % creates frequency coordinate axes 
    F       = sqrt(Fx.^2 + Fy.^2) ; % define absolute frequency value for each point
    
    % Define frequency coefficients as Gaussian random variables
    rng(rs) % sets MATLABS's random number generation (RNG) seed
    PSD      = 0.023 .* r0.^(-5/3) .* F.^(-11/3)                            ; % defines the power spectral density (PSD)
    cnm_high = (randn(N) + 1i*randn(N)) .* sqrt(PSD./(Lx.*Lx)) ; % alter the variance of the randomly sampled coefficients according to the PSD
    cnm_high(floor(N/2)+1,floor(N/2)+1) = 0                                 ; % set zero frequency contributions to zero 
    
    % the randn() function pulls from the standard Gaussian/normal
    % distribution which has the form (1/sqrt(2*pi))*exp(0.5*((x-mu)/sigma)^2) 
    % we therefore multiply by 1/sqrt(2) so that the variance of the 
    % distribution is 1 before wer multiply by the PSD

    screen_high = ifftshift(ifft2(ifftshift(cnm_high)))*N^2;    
    screen = screen_high;  

    if numsub > 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Subharmonic Sampling %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        clc
        rng(rs) % sets random number generation (RNG) seed
        screen_low = zeros(N);
        
        for b =1:numsub
            dfs       = (1/(3^b))*df          ; % spacing in the frequency domain (inverse lengt of the spatial domain as defined by a Fourier transform)
            fs        = [-1 0 1].*dfs         ; % Creates a vector of length N (One entry for each pixel) and centres it so that the middle entry acts as the origin/centre of our one coordinate axes
            [Fxs,Fys] =  meshgrid(fs)         ; % creates frequency coordinate axes 
            Fs        = sqrt(Fxs.^2 + Fys.^2) ; % define absolute frequency value for each point
            
            PSDs         = 0.023 .* r0.^(-5/3) .* Fs.^(-11/3) ;
            cnm_low      = (randn(3) + 1i*randn(3)) .* sqrt(PSDs./(Lx.*Lx)).*(1/3^(b));
            cnm_low(2,2) = 0;
            
            sh = 0;
            for n = 1:9
                sh = sh + cnm_low(n) .* exp(1i *2*pi* (Fxs(n).*X + Fys(n).*Y));
            end
            screen_low = screen_low + sh;
        
        end
        screen = screen_high + screen_low;
    end
        
    screen = real(screen) - mean(mean(real(screen)));

end