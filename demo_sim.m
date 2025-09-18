% =========================================================================
% Introduction
% =========================================================================
% This code provides a simple simulation demonstration of wavefront 
% reconstruction of complex optical fields via spatial and Fourier-domain 
% regularized inversion (SAFARI).
%
% Author: Yunhui Gao (gyh21@mails.tsinghua.edu.cn)
% =========================================================================
%%
% =========================================================================
% Simulation settings
% =========================================================================
close all;
clc;clear;

% load functions
addpath(genpath('./utils'))

% key parameters
params.pxsize = 2.4e-3;         % pixel size (mm)
params.wavlen = 0.532e-3;       % wavelength (mm)

n = 1000;                       % wavefront dimension

% =========================================================================
% Wavefront generation (uncomment the corresponding block)
% =========================================================================

% speckle fields
rng(0)                          % random seed, for reproducibility
grain_size = 8;                 % speckle grain size (pixel)
m = round(n/grain_size);        % number of speckles in one dimension
m = round((m+1)/2)*2;           % make m an even number
u = exp(1i*rand(m,m)*2*pi);     % random phase uniformly sampled in [0,2pi)
wavefront = fftshift(fft2(fftshift(zeropad(u,(n-m)/2)))); % fft to obtain speckles

% % Laguerre Gaussian beams
% [X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
% z = 0;      % z position (mm)
% w0 = 0.2;   % beam waist (mm)
% l = 3;      % azimuthal index
% p = 3;      % radial index
% wavefront = genLaguerreGaussian(X,Y,z,params.wavlen,w0,l,p);

% % Hermite Gaussian beams
% [X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
% z = 0;      % z position (mm)
% w0 = 0.2;   % beam waist (mm)
% mi = 3;     % mode index
% ni = 3;     % mode index
% wavefront = genHermiteGaussian(X,Y,z,params.wavlen,w0,mi,ni);

% % Ince Gaussian beams
% z = 0;      % z position (mm)
% w0 = 0.2;   % beam waist (mm)
% p = 12;     % order 
% m = 8;      % degree
% e = 2;      % ellipticity parameter
% parity = 0; % parity of the beam (0: even, 1:odd)
% wavefront = genInceGaussian(params.pxsize*n/2,n+1,parity,p,m,e,w0,2*pi/params.wavlen,z);
% wavefront = wavefront(1:n,1:n);

% % Ariy beams
% [X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
% w0 = 0.2;                       % scaling factor
% x0 = -n*0.4*params.pxsize;      % x center location (mm)
% y0 = -n*0.4*params.pxsize;      % y center location (mm)
% a = 1e-3;                       % exponential truncation factor
% wavefront = genAiry(X,Y,w0,x0,y0,a);

% % Bessel beams
% [X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
% z = 0;              % z position (mm)
% n_charge = 3;       % topological charge
% theta = pi/3e3;     % axicon angle (rad)
% wavefront = genBessel(X,Y,z,params.wavlen,n_charge,theta);

% % parabolic phase
% [X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
% f = 200;    % focal length (mm)
% a = 0.2;    % amplitude attenuation
% k = 2*pi/params.wavlen;     % wave number
% wavefront = exp(-a*(X.^2+Y.^2)) .* exp(-1i*k*(X.^2+Y.^2)/2/f);

% % Zernike aberrations
% [X,Y] = meshgrid(linspace(-1,1,n));     % define 2D coordinate
% [theta,r] = cart2pol(X,Y);              % convert to polar coordinate
% idx = r <= 1;                           % define the circular aperture
% n_max = 5;                              % define maximum Zernike order
% s_fac = 4;                              % define scaling factor controlling the dynamic range
% n_modes = (n_max+2)*(n_max+1)/2;        % number of Zernike modes
% z_n = nan(n_modes,1);
% z_m = nan(n_modes,1);
% for i = 0:n_max
%     z_n(i*(i+1)/2+1:(i+1)*(i+2)/2) = i;
%     z_m(i*(i+1)/2+1:(i+1)*(i+2)/2) = -i:2:i;
% end
% rng(0)                                  % set random seed, for reproducibility
% coef = 2*rand(n_modes,1)-1;             % uniformly sample Zernike coefficients between [-1,1]
% zer = zeros(n,n);
% for i = 1:n_modes
%     bfun = zeros(n,n);
%     bfun(idx) = zernfun(z_n(i),z_m(i),r(idx),theta(idx));
%     zer(idx) = zer(idx) + coef(i)*bfun(idx);
% end
% phase = 2*pi*s_fac*zer;                 % scale the phase profile
% phase = imresize(phase(ceil(n/2-sqrt(2)/4*n+1):floor(n/2+sqrt(2)/4*n-1),...
%     ceil(n/2-sqrt(2)/4*n+1):floor(n/2+sqrt(2)/4*n-1)),[n,n]);   % crop the central rectangular region
% a = 1;    % amplitude attenuation
% wavefront = exp(-a*(X.^2 + Y.^2)).*exp(1i*phase);

% % turbulence (single phase screen)
% seed = 0;       % set random seed, for reproducibility
% numsub = 10;    % number of subharmonics for subharmonic sampling
% r0 = 0.1;       % Fried parameter (mm)
% phase = FourierPhaseScreen(n,params.pxsize,r0,seed,numsub);     % calculate phase screen
% [X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
% a = 0.1;        % amplitude attenuation
% wavefront = exp(-a*(X.^2 + Y.^2)).*exp(1i*phase);

% % turbulence (scintillated wavefront) [Note: Optimization Toolbox is requried]
% rng(0)          % set random seed, for reproducibility
% n_screen = 10;  % number of phase screens
% D1 = 5;         % length of one side of square phase screen (m)
% D2 = 1;         % diameter of the observation aperture (m)
% Dz = 2e5;       % propagation distance (m)
% wavefront = genTurbulence(n,n_screen,D1,D2,Dz,params.wavlen*1e-3);

% % amplitude pattern
% img = im2double(imread('data/thulogo.bmp'));
% img = imresize(img,[n/2,n/2]);
% img = padarray(1-img,[n/4,n/4],0);
% amp = img;
% pha = zeros(n,n);
% wavefront = amp.*exp(1i*pha);

% % phase pattern
% img = im2double(imread('data/cityulogo.bmp'));
% img = imresize(img,[n/2,n/2]);
% img = padarray(1-img,[n/4,n/4],0);
% amp = ones(n,n);
% pha = img*pi;
% wavefront = amp.*exp(1i*pha);

% % prism
% [X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
% kmax = 2*pi/params.pxsize/4/2;
% kx = kmax*0.5;
% ky = kmax*0.0;
% pha = kx*X + ky*Y;
% amp = ones(n,n);
% wavefront = amp.*exp(1i*pha);

% =========================================================================
% Wavefront tilt / propagation
% =========================================================================

% wavefront tilt
[X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);
k_alpha_x = 0.0;    % x tilt angle (degree)
k_alpha_y = 0.0;    % y tilt angle (degree)
kx = 2*pi/params.wavlen*sind(k_alpha_x);
ky = 2*pi/params.wavlen*sind(k_alpha_y);
pha = kx*X + ky*Y;
wavefront = wavefront.*exp(1i*pha);

% wavefront propagation
prop_dist = 0;      % propagation distance (mm)
wavefront = propagate(wavefront,prop_dist,params.pxsize,params.wavlen);

% =========================================================================
% Measurement simulation
% =========================================================================

% normalization of the amplitude
wavefront = wavefront ./ max(abs(wavefront(:)));

% physical parameters
params.dist = 3.6;      % diffuser-to-sensor distance (mm)
cropsize = 50;          % image cropping size

% define the diffuser profile (random binary phase modulation)
rng(0)
diff_feat_size = 5;
diffuser = imresize(rand(floor(n/diff_feat_size),floor(n/diff_feat_size)),[n,n],'nearest');
index_1 = diffuser <  0.5;
index_2 = diffuser >= 0.5;
diffuser(index_1) = exp(1i*0);
diffuser(index_2) = exp(1i*pi);

% calculate the transfer function for diffraction modeling
HQ = fftshift(transfunc_propagate(n,n, params.dist,params.pxsize,params.wavlen)); % forward propagation

% define function handles for the measurement operators
M  = @(x) x.*diffuser;                  % diffuser modulation
MH = @(x) x.*conj(diffuser);            % Hermitian operator of M
Q  = @(x) ifft2(fft2(x).*HQ);           % free-space propagation
QH = @(x) ifft2(fft2(x).*conj(HQ));     % Hermitian operator of Q
C  = @(x) imgcrop(x,cropsize);          % image cropping to avoid circular convolution artifacts
CH = @(x) zeropad(x,cropsize);          % Hermitian operator of C
A  = @(x) C(Q(M(x)));                   % overall measurement operator
AH = @(x) MH(QH(CH(x)));                % Hermitian operator of A

% generate measurement data
rng(0)                                  % random seed, for reproducibility
x = wavefront;
y = abs(A(x)).^2;
y = y./max(y(:));
y = max(y + 1e-3*randn(size(y)),0);     % simulate noisy measurement

% display the wavefront and measurement
figure,set(gcf,'Units','Normalized','Position',[0.1,0.25,0.8,0.5],'Color','w')
cmap_g = colormap('gray'); cmap_i = inferno; cmap_s = sinebow(256);
ax = subplot(1,3,1);imshow(abs(wavefront),[0,1]);colormap(ax,cmap_g);colorbar
title('Amplitude of the wavefront','fontsize',10)
ax = subplot(1,3,2);imshow(angle(wavefront),[-pi,pi]);colormap(ax,cmap_i);colorbar
title('Phase of the wavefront','fontsize',10)
ax = subplot(1,3,3);imshow(y,[0,1]);colormap(ax,cmap_g);colorbar
title('Encoded intensity image','fontsize',10)
drawnow;

%%
% =========================================================================
% Wavefront reconstruction algorithm
% =========================================================================

% set the regularization parameters
lam_c = [1e-2,1e-2];        % regularization parameter (complex amplitude) [start, end]
lam_a = [1e-2,1e-2];        % regularization parameter (amplitude) [start, end]
alpha = 1e-2;               % parameter tuning weight (controlling the varying speed between the starting and ending values)

% set the support region in the Fourier domain
R = n/8;    % radius of the support region
support = aperture(n,n,n/2,n/2,R);

% Lipschitz bound for the fidelity term
LF = 1*max(abs(diffuser(:)))^2;

% algorithm settings
gpu = 1;            % whether using GPU or not
display = false;    % whether display the intermediate results during iterations (disable to speed up)
verbose = true;     % whether print status in the command line

term_std = 1e-4;    % termination criteria
cache_iter = 10;    % number of previous loss values stored during the reconstruction
cache_loss = inf(cache_iter,1);

% optimization variables
seed = 0;
rng(seed);
x_est = ones(n,n).*exp(1i*2*pi*rand(n,n));
z_est = x_est;

% initialize GPU
if gpu
    device   = gpuDevice(gpu);
    reset(device)
    y        = gpuArray(y);
    HQ       = gpuArray(HQ);
    diffuser = gpuArray(diffuser);
    support  = gpuArray(support);
    x_est = gpuArray(x_est);
    z_est = gpuArray(z_est);
end

% main loop
if display; figure,set(gcf,'Units','Normalized','Position',[0.25,0.3,0.5,0.4],'Color','w'); end
timer = tic;
loss  = inf;
iter = 1;

while(true)

    % update regularization parameters
    lam_a_val = (lam_a(1)-lam_a(2)) * exp(-alpha*(iter-1)) + lam_a(2);
    lam_c_val = (lam_c(1)-lam_c(2)) * exp(-alpha*(iter-1)) + lam_c(2);

    % gradient calculation
    u = A(z_est);
    res = abs(u) - sqrt(y);
    u = u.*(res./abs(u));
    g = AH(u);
    g = g + lam_c_val * DTf(Df(z_est));
    g = g + lam_a_val * exp(1i*angle(z_est)) .* DTf(Df(abs(z_est)));
    
    % proximal gradient update
    gamma = 1/(LF + (lam_a_val + lam_c_val)*8);     % step size (see the paper for details)
    u = z_est - gamma*g;
    x_next = ifft2(fftshift(support.*fftshift(fft2(u))));
    
    z_est = x_next + iter/(iter+3)*(x_next - x_est);      % Nesterov extrapolation
    x_est = x_next;
    
    loss = norm(res(:),2)^2;
    criterion = std(cache_loss)/mean(cache_loss);

    % print and display status
    if verbose
        fprintf('iter: %4d | loss: %5.2e | loss std: %5.2e | runtime: %5.1f s\n', iter, loss, criterion, toc(timer));
    end
    if display
        ax = subplot(1,2,1);imshow(abs(C(x_est)),[]);colormap(ax,cmap_g);colorbar
        title('Reconstructed amplitude','fontsize',10)
        ax = subplot(1,2,2);imshow(angle(C(x_est)),[-pi,pi]);colormap(ax,cmap_i);colorbar
        title('Reconstructed phase','fontsize',10)
        drawnow
    end
    
    % terminate if criterion is met
    if criterion < term_std
        fprintf('Terminated.\n')
        break;
    end
    
    % update loss cache
    idx = rem(iter-1,cache_iter)+1;
    cache_loss(idx) = loss;

    iter = iter+1;

end

% wait for GPU
if gpu; wait(device); end
toc(timer)

% gather data from GPU
if gpu
    x_est    = gather(x_est);
    y        = gather(y);
    HQ       = gather(HQ);
    diffuser = gather(diffuser);
    support  = gather(support);
end

% =========================================================================
% Display results
% =========================================================================

a_max = prctile(abs(x_est(:)),99.99);

figure,set(gcf,'Units','Normalized','Position',[0.25,0.3,0.5,0.4],'Color','w')
ax = subplot(1,2,1);imshow(abs(C(x_est)),[]);colormap(ax,cmap_g);colorbar
title('Reconstructed amplitude','fontsize',10)
ax = subplot(1,2,2);imshow(angle(C(x_est)),[-pi,pi]);colormap(ax,cmap_i);colorbar
title('Reconstructed phase','fontsize',10)

figure,set(gcf,'Units','Normalized','Position',[0.2,0.3,0.6,0.4],'Color','w')
img = visualizeComplex(C(x),cmap_s,[0,1],'hsv',false);
subplot(1,3,1),imshow(img);
title('Ground-truth wavefront','fontsize',10)
img = visualizeComplex(C(x_est),cmap_s,[0,a_max],'hsv',false);
subplot(1,3,2),imshow(img);
title('Reconstructed wavefront','fontsize',10)
img = visualizeComplex(C(x_est.*exp(-1i*angle(x))),cmap_s,[0,a_max],'hsv',false);
subplot(1,3,3),imshow(img);
title('Phase-conjugated wavefront','fontsize',10)
drawnow;

%%
% =========================================================================
% Auxiliary functions
% =========================================================================

function u = imgcrop(x,cropsize)
% =========================================================================
% Crop the central part of the image.
% -------------------------------------------------------------------------
% Input:    - x        : Original image.
%           - cropsize : Cropping pixel number along each dimension.
% Output:   - u        : Cropped image.
% =========================================================================
u = x(cropsize+1:end-cropsize,cropsize+1:end-cropsize);
end


function u = zeropad(x,padsize)
% =========================================================================
% Zero-pad the image.
% -------------------------------------------------------------------------
% Input:    - x        : Original image.
%           - padsize  : Padding pixel number along each dimension.
% Output:   - u        : Zero-padded image.
% =========================================================================
u = padarray(x,[padsize,padsize],0);
end


function H = transfunc_propagate(n1, n2, dist, pxsize, wavlen)
% =========================================================================
% Calculate the transfer function of the free-space diffraction.
% -------------------------------------------------------------------------
% Input:    - n1, n2   : The image dimensions (pixel).
%           - dist     : Propagation distance.
%           - pxsize   : Pixel (sampling) size.
%           - wavlen   : Wavelength of the light.
% Output:   - H        : Transfer function.
% =========================================================================
% sampling in the spatial frequency domain
k1 = pi/pxsize*(-1:2/n1:1-2/n1);
k2 = pi/pxsize*(-1:2/n2:1-2/n2);
[K2,K1] = meshgrid(k2,k1);

k = 2*pi/wavlen;    % wave number

ind = (K1.^2 + K2.^2 >= k^2);  % remove evanescent orders
K1(ind) = 0; K2(ind) = 0;

H = exp(1i*dist*sqrt(k^2-K1.^2-K2.^2));
end


function w = Df(x)
% =========================================================================
% Calculate the 2D gradient (finite difference) of an input image.
% -------------------------------------------------------------------------
% Input:    - x  : The input 2D array.
% Output:   - w  : The gradient (3D array).
% =========================================================================
w = cat(3, x(1:end,:) - x([2:end,end],:), ...
           x(:,1:end) - x(:,[2:end,end]));
end


function u = DTf(w)
% =========================================================================
% Calculate the transpose of the gradient operator.
% -------------------------------------------------------------------------
% Input:    - w  : 3D array.
% Output:   - u  : 2D array.
% =========================================================================
u1 = w(:,:,1) - w([end,1:end-1],:,1);
u1(1,:) = w(1,:,1);
u1(end,:) = -w(end-1,:,1);

u2 = w(:,:,2) - w(:,[end,1:end-1],2);
u2(:,1) = w(:,1,2);
u2(:,end) = -w(:,end-1,2);

u = u1 + u2;
end
