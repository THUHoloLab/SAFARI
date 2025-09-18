clear all;
close all;
clc;

addpath(genpath('../../'))

D = 2;    % length of one side of square phase screen [m]
r0 = 0.1; % coherence diameter [m]
N = 512;  % number of grid points per side
L0 = 100; % outer scale [m]
l0 = 0.01;% inner scale [m]

delta = D/N;    % grid spacing [m]
% spatial grid
x = (-N/2 : N/2-1) * delta;
y = x;
% generate a random draw of an atmospheric phase screen
[phz_lo phz_hi] ...
    = ft_sh_phase_screen(r0, N, delta, L0, l0);
phz = phz_lo + phz_hi;

% example_pt_source_atmos_setup.m

% determine geometry
D2 = 0.5;   % diameter of the observation aperture [m]
wvl = 1e-6;  % optical wavelength [m]
k = 2*pi / wvl; % optical wavenumber [rad/m]
Dz = 200e3;     % propagation distance [m]

% use sinc to model pt source
DROI = 4 * D2;  % diam of obs-plane region of interest [m]
D1 = wvl*Dz / DROI;    % width of central lobe [m]
R = Dz; % wavefront radius of curvature [m]

% atmospheric properties
Cn2 = 1e-16;    % structure parameter [m^-2/3], constant
% SW and PW coherence diameters [m]
r0sw = (0.423 * k^2 * Cn2 * 3/8 * Dz)^(-3/5);
r0pw = (0.423 * k^2 * Cn2 * Dz)^(-3/5);
p = linspace(0, Dz, 1e3);
% log-amplitude variance
rytov = 0.563 * k^(7/6) * sum(Cn2 * (1-p/Dz).^(5/6) ...
    .* p.^(5/6) * (p(2)-p(1)));

% screen properties
nscr = 11; % number of screens
A = zeros(2, nscr); % matrix
alpha = (0:nscr-1) / (nscr-1);
A(1,:) = alpha.^(5/3);
A(2,:) = (1 - alpha).^(5/6) .* alpha.^(5/6);
b = [r0sw.^(-5/3); rytov/1.33*(k/Dz)^(5/6)];
% initial guess
x0 = (nscr/3*r0sw * ones(nscr, 1)).^(-5/3);
% objective function
fun = @(X) sum((A*X(:) - b).^2);
% constraints
x1 = zeros(nscr, 1);
rmax = 0.1; % maximum Rytov number per partial prop
x2 = rmax/1.33*(k/Dz)^(5/6) ./ A(2,:);
x2(A(2,:)==0) = 50^(-5/3);
[X,fval,exitflag,output] ...
    = fmincon(fun,x0,[],[],[],[],x1,x2);
% check screen r0s
r0scrn = X.^(-3/5);
r0scrn(isinf(r0scrn)) = 1e6;
% check resulting r0sw & rytov
bp = A*X(:); [bp(1)^(-3/5) bp(2)*1.33*(Dz/k)^(5/6)];
[r0sw rytov];


% analysis_pt_source_atmos_samp.m

c = 2;
D1p = D1 + c*wvl*Dz/r0sw;
D2p = D2 + c*wvl*Dz/r0sw;

delta1 = linspace(0, 1.1*wvl*Dz/D2p, 100);
deltan = linspace(0, 1.1*wvl*Dz/D1p, 100);
% constraint 1
deltan_max = -D2p/D1p*delta1 + wvl*Dz/D1p;
% constraint 3
d2min3 = (1+Dz/R)*delta1 - wvl*Dz/D1p;
d2max3 = (1+Dz/R)*delta1 + wvl*Dz/D1p;
[delta1 deltan] = meshgrid(delta1, deltan);
% constraint 2
N2 = (wvl * Dz + D1p*deltan + D2p*delta1) ...
    ./ (2 * delta1 .* deltan);

% constraint 4
d1 = 10e-3;
d2 = 10e-3;
% N = 512;
d1*d2 * N / wvl;
zmax = min([d1 d2])^2 * N / wvl;
nmin = ceil(Dz / zmax) + 1;

% example_pt_source_vac_prop.m

delta1 = d1;    % source-plane grid spacing [m]
deltan = d2;    % observation-plane grid spacing [m]
n = nscr;         % number of planes

% coordinates
[x1 y1] = meshgrid((-N/2 : N/2-1) * delta1);
[theta1 r1] = cart2pol(x1, y1);

% point source
pt = exp(-i*k/(2*R) * r1.^2) / D1^2 ...
    .* sinc(x1/D1) .* sinc(y1/D1) ...
    .* exp(-(r1/(4*D1)).^2);
% partial prop planes
z = (1 : n-1) * Dz / (n-1);

% simulate vacuum propagation
sg = exp(-(x1/(0.47*N*d1)).^16) ...
    .* exp(-(y1/(0.47*N*d1)).^16);
t = repmat(sg, [1 1 n]);
[xn yn Uvac] = ang_spec_multi_prop(pt, wvl, ...
    delta1, deltan, z, t);
% collimate the beam
Uvac = Uvac .* exp(-i*pi/(wvl*R)*(xn.^2+yn.^2));

% example_pt_source_turb_prop.m

l0 = 0;     % inner scale [m]
L0 = inf;   % outer scale [m]

zt = [0 z];  % propagation plane locations
Delta_z = zt(2:n) - zt(1:n-1);    % propagation distances
% grid spacings
alpha = zt / zt(n);
delta = (1-alpha) * delta1 + alpha * deltan;

% initialize array for phase screens
phz = zeros(N, N, n);
nreals = 1;    % number of random realizations
% initialize arrays for propagated fields,
% aperture mask, and MCF
Uout = zeros(N);
mask = circ(xn/D2, yn/D2, 1);
MCF2 = zeros(N);
sg = repmat(sg, [1 1 n]);
for idxreal = 1 : nreals     % loop over realizations
    if k>1; fprintf(repmat('\b',1,9)); end
    fprintf('%03d / %03d',idxreal,nreals);
    idxreal;
    % loop over screens
    for idxscr = 1 : 1 : n
        [phz_lo phz_hi] ...
            = ft_sh_phase_screen ...
            (r0scrn(idxscr), N, delta(idxscr), L0, l0);
        phz(:,:,idxscr) = phz_lo + phz_hi;
    end
    % simulate turbulent propagation
    [xn yn Uout] = ang_spec_multi_prop(pt, wvl, ....
        delta1, deltan, z, sg.*exp(i*phz));
    % collimate the beam
    Uout = Uout .* exp(-i*pi/(wvl*R)*(xn.^2+yn.^2));
    % accumulate realizations of the MCF
    MCF2 = MCF2 + corr2_ft(Uout, Uout, mask, deltan);
end
% modulus of the complex degree of coherence
MCDOC2 = abs(MCF2) / (MCF2(N/2+1,N/2+1));

figure,imshow(abs(Uout),[])
figure,imshow(angle(Uout),[])

figure,imshow(abs(Uvac),[])
figure,imshow(angle(Uvac),[])