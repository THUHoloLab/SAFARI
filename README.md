# **SAFARI**: **S**patial **A**nd **F**ourier-dom**A**in **R**egularized **I**nversion for complex wavefront characterization

Authors: **[Yunhui Gao](https://github.com/Yunhui-Gao)**, **[Liangcai Cao](https://scholar.google.com/citations?user=FYYb_-wAAAAJ&hl=en)**, and **[Din Ping Tsai](https://www.cityu.edu.hk/stfprofile/dptsai.htm)**

:school: *[**HoloLab**](http://www.holoddd.com/), Tsinghua University*, and *[**Meta-Devices Lab**](https://dinpingtsai.wixsite.com/mysite), City University of Hong Kong*



<p align="left">
<img src="imgs/fig1.png", width='800'>
</p>


**SAFARI** is a general inverse problem framework for wavefront recovery. It leverages the intrinsic physical properties of wavefronts and enables single-shot characterization of diverse classes of complex optical wavefronts, including **structured beams**, **aberrations**, and **speckle fields**.

## News

- 2025.06 :fire: MATLAB code released!


## Requirements
- The code has been implemented using Matlab 2022b. Older visions may be sufficient but have not been tested.

## Quick Start
Run [`demo_sim.m`](https://github.com/THUHoloLab/SAFARI/blob/master/demo_sim.m) with default parameters. Set `gpu = true;` in the code to enable GPU usage.

We have included the example code for generating different types of wavefront, which are listed below. To simulate the desired wavefront, simply uncomment the corresponding block in [`demo_sim.m`](https://github.com/THUHoloLab/SAFARI/blob/master/demo_sim.m).

<details>
<summary> Generation of <strong>speckel fields</strong></summary>

```matlab
rng(1)                          % random seed, for reproducibility
grain_size = 16;                % speckle grain size (pixel)
m = round(n/grain_size);        % number of speckles in one dimension
u = exp(1i*rand(m,m)*2*pi);     % random phase uniformly sampled in [0,2pi)
wavefront = fftshift(fft2(fftshift(zeropad(u,(n-m)/2)))); % fft to obtain speckles
```

</details>

<details>
<summary> Generation of <strong>Laguerre-Gaussian beams</strong></summary>

```matlab
[X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
z = 0;      % z position (mm)
w0 = 0.5;   % beam waist (mm)
l = 3;      % azimuthal index
p = 3;      % radial index
wavefront = genLaguerreGaussian(X,Y,z,params.wavlen,w0,l,p);
```

</details>

<details>
<summary> Generation of <strong>Hermite-Gaussian beams</strong></summary>

```matlab
[X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
z = 0;      % z position (mm)
w0 = 0.5;   % beam waist (mm)
m = 3;      % mode index
n = 3;      % mode index
wavefront = genHermiteGaussian(X,Y,z,params.wavlen,w0,m,n);
```

</details>

<details>
<summary> Generation of <strong>Ince-Gaussian beams</strong></summary>

```matlab
[X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
z = 0;      % z position (mm)
w0 = 0.5;   % beam waist (mm)
p = 12;     % order 
m = 8;      % degree
e = 2;      % ellipticity parameter
parity = 0; % parity of the beam (0: even, 1:odd)
wavefront = genInceGaussian(params.pxsize*n/2,n+1,parity,p,m,e,w0,2*pi/params.wavlen,z);
wavefront = wavefront(1:n,1:n);
```

</details>

<details>
<summary> Generation of <strong>Airy beams</strong></summary>

```matlab
[X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
w0 = 0.3;   % scaling factor
x0 = -n/2*params.pxsize + 1;    % x center location (mm)
y0 = -n/2*params.pxsize + 1;    % y center location (mm)
a = 1e-3;   % exponential truncation factor
wavefront = genAiry(X,Y,w0,x0,y0,a);
```

</details>

<details>
<summary> Generation of <strong>Bessel beams</strong></summary>

```matlab
[X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
z = 0;              % z position (mm)
n_charge = 3;       % topological charge
theta = pi/3e3;     % axicon angle (rad)
wavefront = genBessel(X,Y,z,params.wavlen,n_charge,theta);
```

</details>

<details>
<summary> Generation of <strong>parabolic phase</strong></summary>

```matlab
[X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
f = 300;    % focal length (mm)
a = 0.1;    % amplitude attenuation
k = 2*pi/params.wavlen;     % wave number
wavefront = exp(-a*(X.^2+Y.^2)) .* exp(-1i*k*(X.^2+Y.^2)/2/f);
```

</details>

<details>
<summary> Generation of <strong>Zernike aberrations</strong></summary>

```matlab
[X,Y] = meshgrid(linspace(-1,1,n));     % define 2D coordinate
[theta,r] = cart2pol(X,Y);              % convert to polar coordinate
idx = r <= 1;                           % define the circular aperture
n_max = 5;                              % define maximum Zernike order
s_fac = 4;                              % define scaling factor controlling the dynamic range
n_modes = (n_max+2)*(n_max+1)/2;        % number of Zernike modes
z_n = nan(n_modes,1);
z_m = nan(n_modes,1);
for i = 0:n_max
    z_n(i*(i+1)/2+1:(i+1)*(i+2)/2) = i;
    z_m(i*(i+1)/2+1:(i+1)*(i+2)/2) = -i:2:i;
end
rng(1)                                  % set random seed, for reproducibility
coef = 2*rand(n_modes,1)-1;             % uniformly sample Zernike coefficients between [-1,1]
zer = zeros(n,n);
for i = 1:n_modes
    bfun = zeros(n,n);
    bfun(idx) = zernfun(z_n(i),z_m(i),r(idx),theta(idx));
    zer(idx) = zer(idx) + coef(i)*bfun(idx);
end
phase = 2*pi*s_fac*zer;                 % scale the phase profile
phase = imresize(phase(ceil(n/2-sqrt(2)/4*n+1):floor(n/2+sqrt(2)/4*n-1),...
    ceil(n/2-sqrt(2)/4*n+1):floor(n/2+sqrt(2)/4*n-1)),[n,n]);
                                        % crop the central rectangular region
a = 1;    % amplitude attenuation
wavefront = exp(-a*(X.^2 + Y.^2)).*exp(1i*phase);
```

</details>

<details>
<summary> Generation of <strong>Fourier phase screens</strong></summary>

```matlab
seed = 1;       % set random seed, for reproducibility
numsub = 10;    % number of subharmonics for subharmonic sampling
r0 = 0.3;       % Fried parameter (mm)
phase = FourierPhaseScreen(n,params.pxsize,r0,seed,numsub);     % calculate phase screen
[X,Y] = meshgrid((-n/2:n/2-1)*params.pxsize);   % define 2D coordinate (mm)
a = 0.1;        % amplitude attenuation
wavefront = exp(-a*(X.^2 + Y.^2)).*exp(1i*phase);
```

</details>

<details>
<summary> Generation of <strong>turbulence</strong></summary>

```matlab
rng(1)          % set random seed, for reproducibility
n_screen = 10;  % number of phase screens
D1 = 5;         % length of one side of square phase screen (m)
D2 = 1;         % diameter of the observation aperture (m)
Dz = 2e5;       % propagation distance (m)
wavefront = genTurbulence(n,n_screen,D1,D2,Dz,params.wavlen*1e-3);
```

</details>


## Theory and principle

<p align="left">
<img src="imgs/fig2.png", width='800'>
</p>

In a typical wavefront sensor, the relationship between the incident wavefront profile $\bm{x}$ and the measured field amplitude $\bm{y}$ can be described by a generalized forward model of coherent optical imaging systems:
$$
    \bm{y} = \lvert \bm{A} \bm{x} \rvert,
$$
where $\bm{A}$ denotes the measurement matrix of the system. Recovering the complex wavefront from a single exposure is an ill-posed problem known as phase retrieval, and thus requires exploiting additional prior knowledge about the wavefront. To address this, SAFARI employs joint spatial and Fourier-domain regularization, leading to the following inverse problem:
$$
\bm{\hat{x}} = \mathop{\mathrm{argmin}}_{\bm{x}} \, \left\{ \underbrace{\left\| \lvert \bm{Ax} \rvert - \bm{y} \right\|_2^2}_{F(\bm{x})}  + \underbrace{\lambda_a \left\| \bm{D} \lvert \bm{x} \rvert \right\|_2^2}_{R_a(\bm{x})} + \underbrace{\lambda_c \left\| \bm{D} \bm{x} \right\|_2^2}_{R_c(\bm{x})} + \underbrace{I_C(\bm{x})}_{R_f(\bm{x})} \right\},
$$
where $F(\bm{x})$ is the data-fidelity loss term ensuring consistency with the forward model. $R_a(\bm{x})$ and $R_c(\bm{x})$ represent the spatial-domain regularization terms for the amplitude and complex amplitude, respectively, where $\bm{D}$ is the spatial finite difference operator and $\lambda_a$, $\lambda_c$ are the corresponding regularization weights. $R_f(\bm{x})$ is the Fourier-domain regularization term with $I_C$ denoting the $\{0,+\infty\}$-valued indicator function of set $C$, defined as
$$
    C \stackrel{\text{def}}{=} \left\{ \bm{x} = \bm{F}^{-1} \bm{u} \, \vert \, u(k_x,k_y) = 0, k_x^2 + k_y^2 > k_{\text{max}}^2  \right\},
$$
where $\bm{F}$ denotes the Fourier transform, $k_x$ and $k_y$ represent the coordinates in the spatial frequency domain, and $k_{\text{max}}$ defines the radius of the support region. These regularization terms are motivated by the empirical observation that most wavefronts exhibit smooth spatial distributions with energy predominantly concentrated in low-frequency components. The spatial and Fourier-domain regularizers can be interpreted as *soft* and *hard* low-pass filters and, when combined into a unified framework, offer flexibility in handling different types of wavefront with varying degrees of distortions. The resulting non-convex and non-smooth optimization problem can be solved via standard optimizers such as the proximal gradient algorithm.