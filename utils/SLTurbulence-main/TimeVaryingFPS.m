function TV_screens = TimeVaryingFPS(N_large, N, dx, time_steps, dt, v_meanx, v_varx, v_meany, v_vary, r0, rs, numsub)

    % rs - Random seed for repeatability in phase screen generation 
    % r0 - Fried parameter define correlation length of field and thus strength of the turbulence
    % N  - number of pixels in real space for large screen
    % dx - pixel size in real space
    % numsub - number of subharmonics for subharmonic sampling    
    % M - number of pixels in real space for output screens
    % time_steps - number of time events to sample
    % dt - time difference between time setps
    % v_meanx - mean wind velocity in the x direction
    % v_varx - wind velocity variance in the x direction
    % v_meany - mean wind velocity in the y direction 
    % v_vary - wind velocity variance in the y direction
    
    % Generate Large Phase Screen %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % rs - Random seed for repeatability in phase screen generation 
    % r0 - Fried parameter define correlation length of field and thus strength of the turbulence
    % N  - number of pixels in real space
    % dx - pixel size in real space
    % numsub - number of subharmonics for subharmonic sampling    
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % High Frequency Screen %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Real Space Parameters %
    Lx    = N_large*dx               ; % length of our grid in the spatial domain (number of pixels*length of each pixel)
    x     = (-N_large/2 : N_large/2-1).*dx ; % Creates a vector of length N (One entry for each pixel) and centres it so that the middle entry acts as the origin/centre of our one coordinate axes
    [X,Y] = meshgrid(x)        ; % creates spatial coordinate grid
    
    % Fourier Space Parameters %
    df      = 1/(Lx)              ; % spacing in the frequency domain (inverse lengt of the spatial domain as defined by a Fourier transform)
    f       = (-N_large/2 : N_large/2-1).*df  ; % Creates a vector of length N (One entry for each pixel) and centres it so that the middle entry acts as the origin/centre of our one coordinate axes
    [Fx,Fy] =  meshgrid(f)        ; % creates frequency coordinate axes 
    F       = sqrt(Fx.^2 + Fy.^2) ; % define absolute frequency value for each point
    
    % Define frequency coefficients as Gaussian random variables
    rng(rs) % sets MATLABS's random number generation (RNG) seed
    PSD      = 0.023 .* r0.^(-5/3) .* F.^(-11/3)                            ; % defines the power spectral density (PSD)
    cnm_high = 1/sqrt(2) .* (randn(N_large) + 1i*randn(N_large)) .* sqrt(PSD./(Lx.*Lx)) ; % alter the variance of the randomly sampled coefficients according to the PSD
    cnm_high(floor(N_large/2)+1,floor(N_large/2)+1) = 0                                 ; % set zero frequency contributions to zero 
    
    % the randn() function pulls from the standard Gaussian/normal
    % distribution which has the form (1/sqrt(2*pi))*exp(0.5*((x-mu)/sigma)^2) 
    % we therefore multiply by 1/sqrt(2) so that the variance of the 
    % distribution is 1 before wer multiply by the PSD
    
    screen_high = ifftshift(ifft2(ifftshift(cnm_high)))*N_large^2;    
    screen = screen_high;  
    
    if numsub > 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % Subharmonic Sampling %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        clc
        rng(rs) % sets random number generation (RNG) seed
        screen_low = zeros(N_large);
        
        for b =1:numsub
            dfs       = (1/(3^b))*df          ; % spacing in the frequency domain (inverse lengt of the spatial domain as defined by a Fourier transform)
            fs        = [-1 0 1].*dfs         ; % Creates a vector of length N (One entry for each pixel) and centres it so that the middle entry acts as the origin/centre of our one coordinate axes
            [Fxs,Fys] =  meshgrid(fs)         ; % creates frequency coordinate axes 
            Fs        = sqrt(Fxs.^2 + Fys.^2) ; % define absolute frequency value for each point
            
            PSDs         = 0.023 .* r0.^(-5/3) .* Fs.^(-11/3) ;
            cnm_low      = 1/sqrt(2) .* (randn(3) + 1i*randn(3)) .* sqrt(PSDs./(Lx.*Lx)).*(1/3^(b));
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Shift Phase Screen %%
    %%%%%%%%%%%%%%%%%%%%%%%%
    % clc
    % time_steps = 100; % number of time steps we will be sampling
    % dt = 1e-5; % time difference between each time step
    % v_meanx = 3; % mean velocity for windspeed
    % v_varx = 1; % mean variance in the wind's velocity
    % v_meany = 0; % mean velocity for windspeed
    % v_vary = 2; % mean variance in the wind's velocity
    % N = 2^9; % size of the time varyong phase screens
    
    TV_screens = zeros(N,N,time_steps);
    
    % TV_screens(:,:,1) = screen(end/2 - M/2+1:end/2 + M/2, end/2 - M/2+1:end/2 + M/2) ;
    % 
    % imagesc(TV_screens(:,:,1)); axis image off; colormap jet
    
    rng(rs) % sets random number generation (RNG) seed
    delta_x = 0;
    delta_y = 0;
    for time_ind = 1:time_steps
        
        velx = randn(1)*v_varx + v_meanx;
        vely = randn(1)*v_vary + v_meany;
    
        delta_x = delta_x + velx*dt;
        delta_y = delta_y + vely*dt;
    
        FT_screen = fftshift(fft2(fftshift(screen)));
        
        shift_FT_screen = FT_screen.*exp(1i.*2.*pi.*(Fx.*delta_x)).*exp(1i.*2.*pi.*(Fy.*delta_y));
    
        IFT_screen =  ifftshift(ifft2(ifftshift(shift_FT_screen))) ;
        
        TV_screens(:,:,time_ind) = real(IFT_screen(end/2 - N/2+1:end/2 + N/2, end/2 - N/2+1:end/2 + N/2) );
               
    end


end